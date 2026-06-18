"""Script principal del pipeline de priorización de genes candidatos
para un kit diagnóstico por PCR.

Uso:
    python priorizar_genes.py <proteinas.faa> <"Nombre especie"> [opciones]

Opciones:
    --interpro   <ruta.tsv>    Resultados de InterProScan ya generados.
    --primers    <ruta.tab>    Archivo primers.tab del pipeline previo.
    --restriccion<ruta.tab>    Archivo res_enzimas.tab del pipeline previo.
    --blast-out  <ruta.tsv>    Resultados de BLAST ya generados (omite ejecutar BLAST).
    --evalue     <float>       Umbral de e-value para BLAST (default: 1e-5).
"""

import sys
import argparse
import os

from carga_fasta  import cargar_proteinas, extraer_descripcion
from blast            import (
    ejecutar_blastp,
    parsear_blast_tabular,
    calcular_metricas_especificidad,
)
from dominios         import obtener_dominios
from clasificacion    import (
    clasificar_gen,
    puntuar_amplicon,
    puntuar_tm,
    puntuar_restriccion,
    calcular_score_integrado,
)
from pipeline_previo  import leer_primers, leer_restriccion
from salida           import escribir_tsv_priorizacion, escribir_reporte
from configuracion    import (
    EVALUE_UMBRAL,
    SALIDA_BLAST_DIR,
    ARCHIVO_PRIMERS_DEFAULT,
    ARCHIVO_RESTRICCION_DEFAULT,
)

#no comprendí del todo que realiza esta función
def construir_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Pipeline de priorización de genes para kit diagnóstico PCR.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Ejemplos:\n"
            '  python priorizar_genes.py proteinas.faa "Brucella suis"\n'
            '  python priorizar_genes.py proteinas.faa "Salmonella enterica" '
            '--blast-out blast_resultados/blast_resultados.tsv\n'
        )
    )
    parser.add_argument(
        "faa",
        help="Archivo FASTA de proteínas (ej: salida de Prokka .faa)"
    )
    parser.add_argument(
        "especie",
        help='Nombre científico completo del organismo objetivo (ej: "Brucella suis")'
    )
    parser.add_argument(
        "--interpro",
        metavar="RUTA_TSV",
        default=None,
        help="Archivo TSV de InterProScan ya generado (opcional)"
    )
    parser.add_argument(
        "--primers",
        metavar="RUTA_TAB",
        default=ARCHIVO_PRIMERS_DEFAULT,
        help=f"Archivo primers.tab del pipeline previo (default: {ARCHIVO_PRIMERS_DEFAULT})"
    )
    parser.add_argument(
        "--restriccion",
        metavar="RUTA_TAB",
        default=ARCHIVO_RESTRICCION_DEFAULT,
        help=f"Archivo res_enzimas.tab del pipeline previo (default: {ARCHIVO_RESTRICCION_DEFAULT})"
    )
    parser.add_argument(
        "--blast-out",
        metavar="RUTA_TSV",
        default=None,
        help="Resultados de BLAST ya generados (omite ejecutar blastp)"
    )
    parser.add_argument(
        "--evalue",
        type=float,
        default=EVALUE_UMBRAL,
        help=f"Umbral de e-value para BLAST (default: {EVALUE_UMBRAL})"
    )
    return parser


def obtener_blast(
    ruta_faa:   str,
    blast_out:  str | None,
    evalue:     float,
) -> dict:
    """Ejecuta BLAST o lee resultados ya existentes según los argumentos.
    Si se proporciona --blast-out, se leen esos resultados directamente
    sin ejecutar blastp (útil para reanálisis sin repetir el BLAST).
    devuelve
    dict[str, list[HitBlast]]
        Hits de BLAST filtrados, agrupados por ID de gen.
    """
    if blast_out and os.path.isfile(blast_out):
        print(f"[BLAST] Usando resultados existentes: {blast_out}")
        return parsear_blast_tabular(blast_out, evalue)

    ruta_resultado = ejecutar_blastp(ruta_faa, SALIDA_BLAST_DIR, evalue)
    return parsear_blast_tabular(ruta_resultado, evalue)


def ensamblar_info_gen(
    record,
    resultado_blast,
    resultado_dominios,
    info_primers,
    info_restriccion,
) -> dict:
    """
    Combina toda la información disponible de un gen en un solo diccionario. integra los resultados de
    BLAST, dominios, primers y restricción para calcular el score final.
    """
    gen_id = record.id
    descripcion = extraer_descripcion(record)

    #Datos de BLAST
    if resultado_blast:
        hits_validos       = len(resultado_blast.hits_validos)
        hits_especie       = resultado_blast.hits_especie
        hits_genero        = resultado_blast.hits_genero
        porcentaje_especie = resultado_blast.porcentaje_especie
        porcentaje_genero  = resultado_blast.porcentaje_genero
    else:
        hits_validos = hits_especie = hits_genero = 0
        porcentaje_especie = porcentaje_genero = 0.0

    #Datos de dominios
    texto_dominios = resultado_dominios.descripcion_combinada() if resultado_dominios else ""

    #Clasificación funcional
    categoria, puntuacion_funcional = clasificar_gen(descripcion, texto_dominios)

    #Datos de primers y PCR
    tm_fw = tm_rv = tamaño_amplicon = None
    if info_primers:
        tm_fw          = info_primers.tm_fw
        tm_rv          = info_primers.tm_rv
        tamaño_amplicon= info_primers.tamaño_pb

    puntuacion_pcr, calidad_pcr = puntuar_tm(tm_fw, tm_rv)
    puntuacion_amplicon          = puntuar_amplicon(tamaño_amplicon)
    # La puntuación de PCR combina calidad de Tm y tamaño del amplicón
    puntuacion_pcr_integrada = (puntuacion_pcr + puntuacion_amplicon) / 2

    #Datos de restricción y clonado 
    total_sitios = info_restriccion.total_sitios if info_restriccion else None
    puntuacion_clonacion, clonacion = puntuar_restriccion(total_sitios)

    #Score integrado
    score = calcular_score_integrado(
        porcentaje_genero    = porcentaje_genero,
        porcentaje_especie   = porcentaje_especie,
        puntuacion_funcional = puntuacion_funcional,
        puntuacion_pcr       = puntuacion_pcr_integrada,
        puntuacion_clonacion = puntuacion_clonacion,
    )

    return {
        "gen_id":               gen_id,
        "descripcion":          descripcion,
        "hits_validos":         hits_validos,
        "hits_especie":         hits_especie,
        "hits_genero":          hits_genero,
        "porcentaje_especie":   porcentaje_especie,
        "porcentaje_genero":    porcentaje_genero,
        "dominios":             texto_dominios or "ninguno",
        "categoria_funcional":  categoria,
        "puntuacion_funcional": puntuacion_funcional,
        "tamaño_amplicon":      tamaño_amplicon,
        "tm_fw":                tm_fw,
        "tm_rv":                tm_rv,
        "calidad_pcr":          calidad_pcr,
        "sitios_restriccion":   total_sitios,
        "clonacion":            clonacion,
        "score_total":          score,
    }


def main() -> None:
    """ todas las etapas del pipeline.
        1. Parsear argumentos de línea de comandos.
        2. Cargar proteínas desde el archivo .faa.
        3. Ejecutar / leer resultados de BLAST.
        4. Obtener / leer análisis de dominios (InterProScan).
        5. Leer archivos del pipeline previo (primers y restricción).
        6. Para cada gen: calcular métricas y score integrado.
        7. Escribir priorizacion_genes.tsv y reporte_priorizacion.txt.
    """
    parser = construir_parser()
    args   = parser.parse_args()

    print("\n" + "=" * 60)
    print(" PIPELINE DE PRIORIZACIÓN DE GENES - KIT DIAGNÓSTICO PCR")
    print("=" * 60)
    print(f"Organismo objetivo : {args.especie}")
    print(f"Archivo proteínas  : {args.faa}")
    print("=" * 60 + "\n")

    #Cargar proteína
    print("[1/5] Cargando proteínas...")
    proteinas = cargar_proteinas(args.faa)
    print(f"      {len(proteinas)} proteínas cargadas.\n")

    #BLAST
    print("[2/5] Análisis BLAST...")
    hits_por_gen = obtener_blast(args.faa, args.blast_out, args.evalue)
    print(f"      {len(hits_por_gen)} genes con hits BLAST válidos.\n")

    #Dominios
    print("[3/5] Análisis de dominios (InterProScan)...")
    dominios_por_gen = obtener_dominios(
        ruta_faa=args.faa,
        directorio_salida="interpro_resultados",
        ruta_tsv_existente=args.interpro,
    )
    print(f"      {len(dominios_por_gen)} genes con información de dominios.\n")

    #Archivos del pipeline previo
    print("[4/5] Leyendo archivos del pipeline de diseño de primers...")
    primers_por_gen     = leer_primers(args.primers)
    restriccion_por_gen = leer_restriccion(args.restriccion)
    print()

    #Ensambr info y calcular scores
    print("[5/5] Calculando métricas y scores integrados...")
    genes_info = []

    for record in proteinas:
        gen_id = record.id

        # Calcular métricas BLAST para este gen
        hits = hits_por_gen.get(gen_id, [])
        resultado_blast = calcular_metricas_especificidad(
            gen_id=gen_id,
            hits=hits,
            nombre_especie=args.especie,
        ) if hits else None

        info = ensamblar_info_gen(
            record             = record,
            resultado_blast    = resultado_blast,
            resultado_dominios = dominios_por_gen.get(gen_id),
            info_primers       = primers_por_gen.get(gen_id),
            info_restriccion   = restriccion_por_gen.get(gen_id),
        )
        genes_info.append(info)

    print(f"      {len(genes_info)} genes procesados.\n")

    #Escritura de salidas
    print("Escribiendo archivos de salida...")
    escribir_tsv_priorizacion(genes_info)
    escribir_reporte(genes_info, args.especie)

    # Mostrar top 5 en consola
    genes_ordenados = sorted(genes_info, key=lambda g: g["score_total"], reverse=True)
    print("\n" + "=" * 60)
    print("TOP 5 CANDIDATOS:")
    print("=" * 60)
    for i, gen in enumerate(genes_ordenados[:5], 1):
        print(
            f"  #{i} {gen['gen_id']:<20} "
            f"Score: {gen['score_total']:.3f}  "
            f"Cat: {gen['categoria_funcional']}"
        )
    print("\n¡Pipeline completado!\n")


if __name__ == "__main__":
    main()