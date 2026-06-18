"""
Ejecución de BLASTP contra nr y parseo de resultados tabulares.

Este módulo se encarga de:
    1. Ejecutar blastp en la línea de comandos para cada proteína.
    2. Leer el archivo tabular generado.
    3. Filtrar hits por e-value.
    4. Contar hits por especie y género del organismo objetivo.

Formato tabular usado (outfmt 6 personalizado):
    qseqid  → ID de la proteína consultada
    sseqid  → ID del hit en la base de datos
    pident  → porcentaje de identidad
    length  → longitud del alineamiento
    evalue  → e-value del hit
    bitscore→ bitscore del hit
    stitle  → título/descripción completa del hit (incluye el organismo)
    stitle contiene el nombre del organismo en formato libre, por ejemplo:
    "putative virulence protein [Brucella suis 1330]"
    A partir de ahí extraemos la especie y el género.
"""

import os
import subprocess
from dataclasses import dataclass, field

from configuracion import EVALUE_UMBRAL, MAX_HITS_BLAST, BLAST_OUTFMT, BLAST_DB


@dataclass
class HitBlast:
    """
    Estructura que representa un hit individual de BLAST.

    Atributos
    
    qseqid   : ID de la proteína consultada.
    sseqid   : ID del hit en nr.
    pident   : Porcentaje de identidad (0–100).
    length   : Longitud del alineamiento en aminoácidos.
    evalue   : E-value del hit (cuanto más bajo, mejor).
    bitscore : Bitscore (cuanto más alto, mejor).
    stitle   : Descripción completa del hit, incluye nombre del organismo.
    """
    qseqid:   str
    sseqid:   str
    pident:   float
    length:   int
    evalue:   float
    bitscore: float
    stitle:   str


@dataclass
class ResultadoBlast:
    """
    agrupa todos los hits válidos de una proteína y
    las métricas de especificidad calculadas a partir de ellos.

    Atributos
    gen_id           : Identificador de la proteína analizada.
    hits_validos     : Lista de HitBlast que pasaron el filtro de e-value.
    hits_especie     : Cantidad de hits del organismo objetivo exacto.
    hits_genero      : Cantidad de hits del mismo género (incluye la especie).
    porcentaje_especie : (hits_especie / total_hits) × 100
    porcentaje_genero  : (hits_genero  / total_hits) × 100
    """
    gen_id:              str
    hits_validos:        list[HitBlast] = field(default_factory=list)
    hits_especie:        int = 0
    hits_genero:         int = 0
    porcentaje_especie:  float = 0.0
    porcentaje_genero:   float = 0.0


def ejecutar_blastp(
    ruta_faa: str,
    directorio_salida: str,
    evalue: float = EVALUE_UMBRAL,
    max_hits: int = MAX_HITS_BLAST,
) -> str:
    """
    Ejecuta blastp contra la base de datos nr para todas las proteínas del archivo.

    Guarda los resultados en un archivo tabular dentro de directorio_salida.

    Parámetros
    
    ruta_faa          : Ruta al archivo FASTA de proteínas.
    directorio_salida : Carpeta donde se guarda el resultado tabular.
    evalue            : Umbral de e-value máximo.
    max_hits          : Número máximo de hits a recuperar por consulta.
    """
    os.makedirs(directorio_salida, exist_ok=True)
    ruta_salida = os.path.join(directorio_salida, "blast_resultados.tsv")

    comando = [
        "blastp",
        "-query",    ruta_faa,
        "-db",       BLAST_DB,
        "-out",      ruta_salida,
        "-outfmt",   BLAST_OUTFMT,
        "-evalue",   str(evalue),
        "-max_target_seqs", str(max_hits),
        "-remote",          # usar BLAST remoto (NCBI); quitar si se tiene nr local
    ]

    print(f"[BLAST] Ejecutando blastp remoto. Esto puede tardar varios minutos...")
    print(f"[BLAST] Comando: {' '.join(comando)}")

    resultado = subprocess.run(comando, capture_output=True, text=True)

    if resultado.returncode != 0:
        raise RuntimeError(
            f"blastp falló con código {resultado.returncode}.\n"
            f"stderr: {resultado.stderr}"
        )

    print(f"[BLAST] Resultados guardados en: {ruta_salida}")
    return ruta_salida


def parsear_blast_tabular(ruta_tsv: str, evalue_umbral: float = EVALUE_UMBRAL) -> dict[str, list[HitBlast]]:
    """
    Lee el archivo tabular de BLAST y devuelve un diccionario de hits por gen.

    El formato esperado es el definido en BLAST_OUTFMT de configuracion.py:
        qseqid  sseqid  pident  length  evalue  bitscore  stitle

    Los campos están separados por tabulaciones. El campo stitle puede
    contener espacios, por eso usamos split con maxsplit=6.

    Parámetros
    ruta_tsv      : Ruta al archivo tabular generado por blastp.
    evalue_umbral : E-value máximo para incluir un hit.

    Retorna
    dict[str, list[HitBlast]]
        Clave: ID del gen consultado.
        Valor: lista de HitBlast válidos (ya filtrados por e-value).
    """
    hits_por_gen: dict[str, list[HitBlast]] = {}

    with open(ruta_tsv) as archivo:
        for numero_linea, linea in enumerate(archivo, start=1):
            linea = linea.strip()
            if not linea or linea.startswith("#"):
                continue

            # Separa en exactamente 7 campos; stitle puede tener espacios
            partes = linea.split("\t", maxsplit=6)
            if len(partes) < 7:
                print(f"[BLAST] Línea {numero_linea} ignorada (formato inesperado): {linea}")
                continue

            qseqid, sseqid, pident, length, evalue, bitscore, stitle = partes

            try:
                evalue_val = float(evalue)
            except ValueError:
                continue

            if evalue_val > evalue_umbral:
                continue

            hit = HitBlast(
                qseqid=qseqid,
                sseqid=sseqid,
                pident=float(pident),
                length=int(length),
                evalue=evalue_val,
                bitscore=float(bitscore),
                stitle=stitle,
            )

            hits_por_gen.setdefault(qseqid, []).append(hit)

    return hits_por_gen


def _extraer_organismo_del_stitle(stitle: str) -> str:
    """
    Extrae el nombre del organismo desde el campo stitle de BLAST. BLAST incluye el organismo entre corchetes al final del stitle.

    Si no hay corchetes, devuelve el stitle completo 
    """
    import re
    # Buscamos el ÚLTIMO par de corchetes en el stitle
    matches = re.findall(r"\[([^\]]+)\]", stitle)
    if matches:
        return matches[-1]  # El último [organismo] es el nombre taxonómico
    return stitle


def calcular_metricas_especificidad(
    gen_id: str,
    hits: list[HitBlast],
    nombre_especie: str,
) -> ResultadoBlast:
    """
    Calcula métricas de especificidad taxonómica a partir de los hits de BLAST.

    La especie objetivo se extrae del nombre científico completo (ej: "Brucella suis").
    El género es la primera palabra del nombre (ej: "Brucella").

    Lógica de conteo:
        - Un hit cuenta para 'especie' si el organismo en stitle contiene
          el nombre completo de la especie objetivo.
        - Un hit cuenta para 'género' si el organismo en stitle contiene
          el nombre del género (primera palabra del nombre científico).
        - Los hits de la especie también se cuentan en el género.

    Retorna ResultadoBlast con todos los campos calculados.
    """
    # El género es la primera palabra del nombre científico
    genero = nombre_especie.split()[0]

    total = len(hits)
    conteo_especie = 0
    conteo_genero  = 0

    for hit in hits:
        organismo = _extraer_organismo_del_stitle(hit.stitle).lower()
        if nombre_especie.lower() in organismo:
            conteo_especie += 1
            conteo_genero  += 1
        elif genero.lower() in organismo:
            conteo_genero += 1

    porcentaje_especie = (conteo_especie / total * 100) if total > 0 else 0.0
    porcentaje_genero  = (conteo_genero  / total * 100) if total > 0 else 0.0

    return ResultadoBlast(
        gen_id=gen_id,
        hits_validos=hits,
        hits_especie=conteo_especie,
        hits_genero=conteo_genero,
        porcentaje_especie=round(porcentaje_especie, 2),
        porcentaje_genero=round(porcentaje_genero, 2),
    )