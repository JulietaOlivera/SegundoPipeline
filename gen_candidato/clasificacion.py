"""
Clasificación funcional automática y sistema de puntuación de genes.
dos responsabilidades: (creo que se supone que deben tener una sola...)

    1. Clasificar genes según su función biológica usando palabras clave.
       La clasificación NO elimina genes, solo les asigna una etiqueta
       que refleja su potencial como candidato diagnóstico.

    2. Calcular el score integrado que combina todas las métricas del pipeline
       en un único número comparable entre genes.
"""

from configuracion import (
    PUNTUACIONES_FUNCIONALES,
    PESO_GENERO, PESO_ESPECIE, PESO_FUNCIONAL, PESO_PCR, PESO_CLONACION,
    AMPLICON_MUY_FAVORABLE_MIN, AMPLICON_MUY_FAVORABLE_MAX, AMPLICON_FAVORABLE_MAX,
    PUNTUACION_AMPLICON_MUY_FAVORABLE, PUNTUACION_AMPLICON_FAVORABLE, PUNTUACION_AMPLICON_DESFAVORABLE,
    TM_DIFERENCIA_OPTIMA, TM_DIFERENCIA_ACEPTABLE,
    PUNTUACION_TM_OPTIMA, PUNTUACION_TM_ACEPTABLE, PUNTUACION_TM_DESFAVORABLE,
    RESTRICCION_IDEAL, RESTRICCION_ACEPTABLE,
    PUNTUACION_RESTRICCION_IDEAL, PUNTUACION_RESTRICCION_ACEPTABLE, PUNTUACION_RESTRICCION_DESFAVORABLE,
)


def clasificar_gen(descripcion: str, texto_dominios: str) -> tuple[str, float]:
    """
    Clasifica un gen según palabras clave en su descripción y dominios.

    Estrategia:
        Se concatena la descripción del gen con el texto de sus dominios
        (ambos en minúsculas) y se buscan las palabras clave definidas en
        PUNTUACIONES_FUNCIONALES (configuracion.py).
        Se toma la categoría con la puntuación más alta encontrada.
        Si no se encuentra ninguna, la categoría es "desconocida" con score 0.
    """
    # Texto unificado para búsqueda (todo en minúsculas)
    texto_busqueda = (descripcion + " " + texto_dominios).lower()

    mejor_categoria  = "desconocida"
    mejor_puntuacion = 0.0

    for palabra_clave, puntuacion in PUNTUACIONES_FUNCIONALES.items():
        if palabra_clave in texto_busqueda:
            # Si hay varias, la de mayor puntuación
            if puntuacion > mejor_puntuacion:
                mejor_puntuacion = puntuacion
                mejor_categoria  = _nombre_legible(palabra_clave)

    return mejor_categoria, float(mejor_puntuacion)


def _nombre_legible(palabra_clave: str) -> str:
    """Convierte una palabra clave interna a una etiqueta legible para el reporte.
    palabra_clave : Clave del diccionario PUNTUACIONES_FUNCIONALES.
    """
    mapa = {
        "virulence":        "virulence-associated",
        "pathogenicity":    "pathogenicity-island",
        "secretion":        "secretion-system",
        "outer membrane":   "outer-membrane-protein",
        "membrane protein": "membrane-protein",
        "transporter":      "transporter",
        "host interaction": "host-interaction",
        "invasion":         "invasion-factor",
        "toxin":            "toxin",
        "adhesin":          "adhesin",
        "ribosomal":        "housekeeping-ribosomal",
        "atp synthase":     "housekeeping-ATPsynthase",
        "dna polymerase":   "housekeeping-DNApol",
        "rna polymerase":   "housekeeping-RNApol",
        "elongation factor":"housekeeping-EF",
        "gyrase":           "housekeeping-gyrase",
        "housekeeping":     "housekeeping",
        "chaperone":        "housekeeping-chaperone",
    }
    return mapa.get(palabra_clave, palabra_clave)



# Funciones de puntuación por componente

def puntuar_amplicon(tamaño_pb: int | None) -> float:
    """ Asigna una puntuación (0–100) según el tamaño del amplicón. Los rangos están definidos en configuracion.py.
    Un amplicón más corto es preferible en PCR diagnóstica porque:
        - Amplifica más eficientemente.
        - Es más robusto en muestras degradadas.
        - Da bandas más claras en gel de agarosa.
        Puntuación entre 0 y 100.
    """
    if tamaño_pb is None:
        return 50.0  # valor neutro cuando no hay información

    if AMPLICON_MUY_FAVORABLE_MIN <= tamaño_pb <= AMPLICON_MUY_FAVORABLE_MAX:
        return PUNTUACION_AMPLICON_MUY_FAVORABLE
    if tamaño_pb <= AMPLICON_FAVORABLE_MAX:
        return PUNTUACION_AMPLICON_FAVORABLE
    return PUNTUACION_AMPLICON_DESFAVORABLE


def puntuar_tm(tm_fw: float | None, tm_rv: float | None) -> tuple[float, str]:
    """Evalúa la diferencia de Tm entre los primers forward y reverse. 
    tuple[float, str]
        (puntuación, etiqueta_calidad)
        Etiqueta puede ser: "óptimo", "aceptable", "desfavorable".
    """
    if tm_fw is None or tm_rv is None:
        return 50.0, "sin datos"

    diferencia = abs(tm_fw - tm_rv)

    if diferencia <= TM_DIFERENCIA_OPTIMA:
        return PUNTUACION_TM_OPTIMA, "óptimo"
    if diferencia <= TM_DIFERENCIA_ACEPTABLE:
        return PUNTUACION_TM_ACEPTABLE, "aceptable"
    return PUNTUACION_TM_DESFAVORABLE, "desfavorable"


def puntuar_restriccion(total_sitios: int | None) -> tuple[float, str]:
    """Puntúa la aptitud del amplicón para clonado según los sitios de restricción.
Para clonar el amplicón en un vector, se usan enzimas de restricción que deben cortar SOLO en los extremos del inserto
    """
    if total_sitios is None:
        return 50.0, "sin datos"

    if total_sitios == RESTRICCION_IDEAL:
        return PUNTUACION_RESTRICCION_IDEAL, "ideal"
    if total_sitios <= RESTRICCION_ACEPTABLE:
        return PUNTUACION_RESTRICCION_ACEPTABLE, "aceptable"
    return PUNTUACION_RESTRICCION_DESFAVORABLE, "desfavorable"


def calcular_score_integrado(
    porcentaje_genero:   float,
    porcentaje_especie:  float,
    puntuacion_funcional:float,
    puntuacion_pcr:      float,
    puntuacion_clonacion:float,
) -> float:
    """Calcula el score integrado que combina todas las métricas en un solo número.

    Fórmula:
        Score = w_genero  × %genero
              + w_especie × %especie
              + w_func    × puntuacion_funcional_normalizada
              + w_pcr     × puntuacion_pcr
              + w_clon    × puntuacion_clonacion

    Todos los pesos se definen en configuracion.py y deben sumar 1.0.

    La puntuación funcional se normaliza de un rango [-10, 20] a [0, 100]
    para que sea comparable con los porcentajes.

    porcentaje_genero    : % de hits en el género objetivo (0–100).
    porcentaje_especie   : % de hits en la especie objetivo (0–100).
    puntuacion_funcional : Puntuación raw de la clasificación (-10 a 20).
    puntuacion_pcr       : Puntuación de calidad de primers (0–100).
    puntuacion_clonacion : Puntuación de aptitud para clonado (0–100).

   va a devolver
    float
        Score integrado, en rango 0–100.
    """
    # Normalizar puntuación funcional: de [-10, 20] a [0, 100]
    # Fórmula: (valor - min) / (max - min) * 100
    puntuacion_funcional_norm = ((puntuacion_funcional - (-10)) / (20 - (-10))) * 100
    puntuacion_funcional_norm = max(0.0, min(100.0, puntuacion_funcional_norm))

    score = (
        PESO_GENERO    * porcentaje_genero   +
        PESO_ESPECIE   * porcentaje_especie  +
        PESO_FUNCIONAL * puntuacion_funcional_norm +
        PESO_PCR       * puntuacion_pcr      +
        PESO_CLONACION * puntuacion_clonacion
    )
    return round(score, 3)