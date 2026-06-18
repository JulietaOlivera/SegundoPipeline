"""
dominios.py
===========
Análisis de dominios proteicos mediante InterProScan.

Este módulo soporta dos modos de operación:

    Modo 1 InterProScan local disponible:
        Se ejecuta interproscan.sh directamente sobre el archivo .faa.
        Requiere que InterProScan esté instalado y en el PATH.

    Modo 2 Leer resultados preexistentes:
        Si el usuario ya corrió InterProScan por su cuenta, se lee el
        archivo de resultados en formato TSV que InterProScan genera.

Formato TSV de InterProScan (los campos):
    Columna 0  → ID de la proteína
    Columna 3  → Base de datos de dominio (Pfam, SMART, TIGRFAM, etc.)
    Columna 4  → Accesión del dominio (ej: PF00005)
    Columna 5  → Descripción del dominio (ej: "ABC transporter")
    Columna 11 → Descripción del término GO (puede estar vacía)

    Referencia: https://interproscan-docs.readthedocs.io/en/latest/OutputFormats.html
"""

import os
import subprocess
from dataclasses import dataclass, field


@dataclass
class Dominio:
    """
    Representa un dominio proteico detectado por InterProScan.

    Atributos
    base_datos   : Base de datos de origen (Pfam, SMART, TIGRFAM, etc.).
    accesion     : Accesión del dominio en esa base de datos.
    descripcion  : Descripción legible del dominio.
    """
    base_datos:  str
    accesion:    str
    descripcion: str


@dataclass
class ResultadoDominios:
    """
    Agrupa todos los dominios encontrados para una proteína.

    Atributos
    gen_id   : Identificador de la proteína.
    dominios : Lista de objetos Dominio detectados.
    """
    gen_id:   str
    dominios: list[Dominio] = field(default_factory=list)

    def descripcion_combinada(self) -> str:
        """
        Devuelve una cadena con todas las descripciones de dominios
        concatenadas, útil para búsqueda de palabras clave.
        """
        return " | ".join(d.descripcion for d in self.dominios)


def verificar_interproscan() -> bool:
    """
    Verifica si interproscan.sh está disponible en el sistema.
    True si el comando existe, False si no está instalado.
    """
    resultado = subprocess.run(
        ["which", "interproscan.sh"],
        capture_output=True,
        text=True,
    )
    return resultado.returncode == 0


def ejecutar_interproscan(ruta_faa: str, directorio_salida: str) -> str: #Ejecuta InterProScan local sobre el archivo de proteínas.

    os.makedirs(directorio_salida, exist_ok=True)
    prefijo_salida = os.path.join(directorio_salida, "interproscan_resultado")

    comando = [
        "interproscan.sh",
        "-i",  ruta_faa,
        "-o",  prefijo_salida,
        "-f",  "TSV",
        "--goterms",       # incluir términos GO
        "--pathways",      # incluir pathways
        "--cpu", "4",
    ]

    print("[InterPro] Ejecutando InterProScan local...")
    resultado = subprocess.run(comando, capture_output=True, text=True)

    if resultado.returncode != 0:
        raise RuntimeError(
            f"InterProScan falló.\nstderr: {resultado.stderr}"
        )

    ruta_tsv = prefijo_salida + ".tsv"
    print(f"[InterPro] Resultados guardados en: {ruta_tsv}")
    return ruta_tsv


def parsear_interproscan_tsv(ruta_tsv: str) -> dict[str, ResultadoDominios]:
    """
    Lee el archivo TSV de InterProScan y devuelve dominios por proteína.

    Formato esperado (columnas separadas por tabulaciones):
        Col 0  : ID proteína
        Col 3  : Base de datos (Pfam, SMART, etc.)
        Col 4  : Accesión del dominio
        Col 5  : Descripción del dominio
        (otras columnas se ignoran)

    Retorna
    dict[str, ResultadoDominios]
        Clave: ID de la proteína.
        Valor: ResultadoDominios con todos los dominios detectados.
    """
    resultados: dict[str, ResultadoDominios] = {}

    with open(ruta_tsv) as archivo:
        for linea in archivo:
            linea = linea.strip()
            if not linea or linea.startswith("#"):
                continue

            campos = linea.split("\t")
            if len(campos) < 6:
                continue

            gen_id     = campos[0]
            base_datos = campos[3]
            accesion   = campos[4]
            descripcion= campos[5]

            dominio = Dominio(
                base_datos=base_datos,
                accesion=accesion,
                descripcion=descripcion,
            )

            if gen_id not in resultados:
                resultados[gen_id] = ResultadoDominios(gen_id=gen_id)
            resultados[gen_id].dominios.append(dominio)

    return resultados


def obtener_dominios(
    ruta_faa: str,
    directorio_salida: str,
    ruta_tsv_existente: str | None = None,
) -> dict[str, ResultadoDominios]:
    """
    Punto de entrada unificado para el análisis de dominios.

    Intenta las siguientes opciones en orden:
        1. Si se pasa ruta_tsv_existente, lee ese archivo directamente.
        2. Si InterProScan está instalado, lo ejecuta.
        3. Si ninguna opción está disponible, devuelve un dict vacío
           y emite una advertencia.

    Parámetros
    ruta_faa            : Ruta al archivo .faa (necesario para modo 2).
    directorio_salida   : Carpeta para guardar resultados de InterProScan.
    ruta_tsv_existente  : Ruta a un TSV de InterProScan ya generado (opcional).
    """
    if ruta_tsv_existente and os.path.isfile(ruta_tsv_existente):
        print(f"[InterPro] Usando resultados existentes: {ruta_tsv_existente}")
        return parsear_interproscan_tsv(ruta_tsv_existente)

    if verificar_interproscan():
        ruta_tsv = ejecutar_interproscan(ruta_faa, directorio_salida)
        return parsear_interproscan_tsv(ruta_tsv)

    print(
        "[AVISO] InterProScan no está disponible y no se proveyó un archivo "
        "TSV existente. El análisis de dominios se omitirá."
    )
    return {}