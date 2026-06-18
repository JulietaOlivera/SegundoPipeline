"""Carga y validación del archivo FASTA de proteínas.

Este módulo lee el archivo .faa (o cualquier FASTA de proteínas) y devuelve
una lista de objetos SeqRecord de Biopython. Cada SeqRecord representa una
proteína y contiene su identificador, descripción y secuencia.

Estructura de un SeqRecord relevante para este pipeline:
    record.id          → identificador (ej: "PROKKA_00001")
    record.description → línea completa del encabezado FASTA
    record.seq         → secuencia de aminoácidos
"""

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def cargar_proteinas(ruta_faa: str) -> list[SeqRecord]:

    try:
        registros = list(SeqIO.parse(ruta_faa, "fasta"))
    except FileNotFoundError:
        raise FileNotFoundError(f"No se encontró el archivo: {ruta_faa}")

    if not registros:
        raise ValueError(
            f"El archivo {ruta_faa} está vacío o no tiene secuencias FASTA válidas."
        )

    return registros


def extraer_descripcion(record: SeqRecord) -> str:
    """
    Extrae solo la parte descriptiva del encabezado FASTA, sin el ID.

    Por ejemplo, si el encabezado es:
        >PROKKA_00001 hypothetical protein
    Esta función devuelve: "hypothetical protein"

    Parámetros
    ----------
    record : SeqRecord
        Objeto SeqRecord de Biopython.

    Retorna
    -------
    str
        Descripción del gen sin el identificador inicial.
    """
    # record.description incluye el ID al principio; se saca
    descripcion = record.description
    if descripcion.startswith(record.id):
        descripcion = descripcion[len(record.id):].strip()
    return descripcion if descripcion else "sin descripción"