"""
Generación de los archivos de salida del pipeline de priorización.

Este módulo escribe dos archivos:

    priorizacion_genes.tsv
        Tabla con todas las métricas de cada gen, ordenada por Score_total
        de mayor a menor. Es el archivo principal para selección del candidato.

    reporte_priorizacion.txt
        Reporte en texto plano con resumen estadístico, justificaciones
        automáticas para los mejores genes y explicación del sistema de score.

Estructura del diccionario `gen_info` que reciben estas funciones:
    {
        "gen_id":               str,
        "descripcion":          str,
        "hits_validos":         int,
        "hits_especie":         int,
        "hits_genero":          int,
        "porcentaje_especie":   float,
        "porcentaje_genero":    float,
        "dominios":             str,     # texto combinado de descripciones
        "categoria_funcional":  str,
        "puntuacion_funcional": float,
        "tamaño_amplicon":      int | None,
        "tm_fw":                float | None,
        "tm_rv":                float | None,
        "calidad_pcr":          str,
        "sitios_restriccion":   int | None,
        "clonacion":            str,
        "score_total":          float,
    }
"""

from datetime import datetime

from configuracion import (
    PESO_GENERO, PESO_ESPECIE, PESO_FUNCIONAL, PESO_PCR, PESO_CLONACION,
    TOP_N_GENES, SALIDA_PRIORIZACION, SALIDA_REPORTE,
)


def escribir_tsv_priorizacion(genes_info: list[dict], ruta_salida: str = SALIDA_PRIORIZACION) -> None:
    """
    Escribe priorizacion_genes.tsv ordenado por Score_total descendente.
    """
    # Ordenar de mayor a menor score
    genes_ordenados = sorted(genes_info, key=lambda g: g["score_total"], reverse=True)

    encabezado = "\t".join([
        "Gen", "Descripcion", "Hits_validos", "Hits_especie", "Hits_genero",
        "Porcentaje_especie", "Porcentaje_genero", "Dominios",
        "Categoria_funcional", "Tamano_amplicon", "Tm_fw", "Tm_rv",
        "Calidad_PCR", "Sitios_restriccion_totales", "Aptitud_clonacion",
        "Score_total",
    ])

    with open(ruta_salida, "w") as archivo:
        archivo.write(encabezado + "\n")
        for gen in genes_ordenados:
            fila = "\t".join([
                str(gen["gen_id"]),
                str(gen["descripcion"]),
                str(gen["hits_validos"]),
                str(gen["hits_especie"]),
                str(gen["hits_genero"]),
                str(gen["porcentaje_especie"]),
                str(gen["porcentaje_genero"]),
                str(gen["dominios"]),
                str(gen["categoria_funcional"]),
                str(gen["tamaño_amplicon"] if gen["tamaño_amplicon"] is not None else "N/D"),
                str(gen["tm_fw"] if gen["tm_fw"] is not None else "N/D"),
                str(gen["tm_rv"] if gen["tm_rv"] is not None else "N/D"),
                str(gen["calidad_pcr"]),
                str(gen["sitios_restriccion"] if gen["sitios_restriccion"] is not None else "N/D"),
                str(gen["clonacion"]),
                str(gen["score_total"]),
            ])
            archivo.write(fila + "\n")

    print(f"[salida] Tabla de priorización guardada en: {ruta_salida}")


def escribir_reporte(
    genes_info:      list[dict],
    nombre_especie:  str,
    ruta_salida:     str = SALIDA_REPORTE,
) -> None:
    """
    Escribe el reporte en texto plano con resumen y justificaciones automáticas.
    """
    genes_ordenados = sorted(genes_info, key=lambda g: g["score_total"], reverse=True)
    top_genes = genes_ordenados[:TOP_N_GENES]

    # Contar categorías funcionales
    conteo_categorias: dict[str, int] = {}
    for gen in genes_info:
        cat = gen["categoria_funcional"]
        conteo_categorias[cat] = conteo_categorias.get(cat, 0) + 1

    with open(ruta_salida, "w") as archivo:
        archivo.write(_seccion_encabezado(nombre_especie))
        archivo.write(_seccion_resumen(genes_info, conteo_categorias))
        archivo.write(_seccion_formula_score())
        archivo.write(_seccion_top_genes(top_genes, nombre_especie))

    print(f"[salida] Reporte guardado en: {ruta_salida}")

# Funciones auxiliares para las secciones del reporte

def _seccion_encabezado(nombre_especie: str) -> str:
    """Genera la sección de encabezado del reporte."""
    fecha = datetime.now().strftime("%Y-%m-%d %H:%M")
    return (
        "=" * 70 + "\n"
        "   REPORTE DE PRIORIZACIÓN DE GENES PARA KIT DIAGNÓSTICO\n"
        "=" * 70 + "\n"
        f"Organismo objetivo : {nombre_especie}\n"
        f"Fecha de análisis  : {fecha}\n"
        "=" * 70 + "\n\n"
    )


def _seccion_resumen(genes_info: list[dict], conteo_categorias: dict) -> str:
    """Genera el resumen estadístico."""
    total = len(genes_info)
    con_blast = sum(1 for g in genes_info if g["hits_validos"] > 0)
    lineas = [
        "RESUMEN GENERAL\n",
        "-" * 40 + "\n",
        f"Total de genes analizados     : {total}\n",
        f"Genes con hits BLAST válidos  : {con_blast}\n",
        f"Genes sin hits BLAST          : {total - con_blast}\n\n",
        "DISTRIBUCIÓN POR CATEGORÍA FUNCIONAL:\n",
    ]
    for cat, conteo in sorted(conteo_categorias.items(), key=lambda x: -x[1]):
        lineas.append(f"  {cat:<35} : {conteo}\n")
    lineas.append("\n")
    return "".join(lineas)


def _seccion_formula_score() -> str:
    """Genera la explicación de la fórmula del score."""
    return (
        "FÓRMULA DEL SCORE INTEGRADO\n"
        "-" * 40 + "\n"
        "Score = "
        f"{PESO_GENERO}  × %_género    "
        f"(conservación en el género objetivo)\n"
        f"      + {PESO_ESPECIE}  × %_especie   "
        f"(conservación en la especie objetivo)\n"
        f"      + {PESO_FUNCIONAL}  × func_norm   "
        f"(relevancia funcional normalizada 0–100)\n"
        f"      + {PESO_PCR}  × cal_PCR     "
        f"(calidad del par de primers)\n"
        f"      + {PESO_CLONACION}  × aptit_clon  "
        f"(aptitud para clonado)\n\n"
        "Todos los pesos son modificables en configuracion.py\n\n"
    )


def _seccion_top_genes(top_genes: list[dict], nombre_especie: str) -> str:
    """Genera la sección de justificaciones para los mejores genes."""
    genero = nombre_especie.split()[0]
    lineas = [
        f"TOP {len(top_genes)} GENES CANDIDATOS\n",
        "=" * 70 + "\n\n",
    ]
    for posicion, gen in enumerate(top_genes, start=1):
        lineas.append(f"#{posicion} — {gen['gen_id']}\n")
        lineas.append(f"    Descripción         : {gen['descripcion']}\n")
        lineas.append(f"    Categoría funcional : {gen['categoria_funcional']}\n")
        lineas.append(f"    % hits en especie   : {gen['porcentaje_especie']:.1f}%\n")
        lineas.append(f"    % hits en género    : {gen['porcentaje_genero']:.1f}%\n")
        lineas.append(f"    Tamaño amplicón     : {gen['tamaño_amplicon'] or 'N/D'} pb\n")
        lineas.append(f"    Calidad PCR         : {gen['calidad_pcr']}\n")
        lineas.append(f"    Sitios restricción  : {gen['sitios_restriccion'] if gen['sitios_restriccion'] is not None else 'N/D'}\n")
        lineas.append(f"    Score total         : {gen['score_total']:.3f}\n")
        lineas.append(f"    Justificación:\n")
        lineas.append(f"        {_generar_justificacion(gen, genero)}\n")
        lineas.append("\n" + "-" * 70 + "\n\n")
    return "".join(lineas)
