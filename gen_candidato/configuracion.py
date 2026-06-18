"""Todas las constantes que controlan el comportamiento del pipeline están
centralizadas aca para que sean fáciles de encontrar y modificar sin
tocar el resto del código.

Estructura:
    - Parámetros de BLAST
    - Pesos del score integrado
    - Puntuaciones funcionales
    - Rangos de tamaño de amplicón
    - Umbrales de diferencia de Tm
    - Umbrales de sitios de restricción
    - Palabras clave para clasificación funcional
"""


# PARÁMETROS DE BLAST
# E-value máximo para considerar un hit como válido.
# Hits con e-value mayor a este umbral son descartados.
EVALUE_UMBRAL: float = 1e-5

# Número máximo de hits a recuperar por consulta BLAST.
# Aumentar este valor da más información pero hace el análisis más lento.
MAX_HITS_BLAST: int = 100

# Formato de salida de BLAST (formato 6 = tabular).
# Cada campo que pedimos está separado por tabulaciones.
# Referencia: https://www.ncbi.nlm.nih.gov/books/NBK153387/
BLAST_OUTFMT: str = (
    "6 qseqid sseqid pident length evalue bitscore stitle"
)

# Base de datos de BLAST a consultar.
BLAST_DB: str = "nr"


# PESOS DEL SCORE INTEGRADO
# Deben sumar 1.0 para que el score esté en escala de 0 a 100.

PESO_GENERO:    float = 0.35   # Conservación en el género objetivo
PESO_ESPECIE:   float = 0.25   # Conservación en la especie objetivo
PESO_FUNCIONAL: float = 0.15   # Relevancia de la función biológica
PESO_PCR:       float = 0.15   # Calidad del par de primers
PESO_CLONACION: float = 0.10   # Aptitud para clonado (sitios de restricción)


# PUNTUACIONES FUNCIONALES
# Se suman o restan al componente funcional del score (escala 0–100).
# Palabras clave que se buscan en la descripción y en los dominios del gen.
# El valor asociado es la puntuación que se suma si se detecta esa categoría.
PUNTUACIONES_FUNCIONALES: dict = {
    # Categorías favorables para diagnóstico
    "virulence":       20,
    "pathogenicity":   18,
    "secretion":       15,
    "outer membrane":  15,
    "membrane protein":12,
    "transporter":     10,
    "host interaction":18,
    "invasion":        16,
    "toxin":           18,
    "adhesin":         14,
    # Categorías desfavorables (genes muy conservados)
    "ribosomal":      -10,
    "atp synthase":   -10,
    "dna polymerase": -10,
    "rna polymerase": -10,
    "elongation factor": -10,
    "gyrase":          -8,
    "housekeeping":    -10,
    "chaperone":        -5,
}

# RANGOS DE TAMAÑO DE AMPLICÓN (en pares de bases)
AMPLICON_MUY_FAVORABLE_MIN: int = 100
AMPLICON_MUY_FAVORABLE_MAX: int = 500
AMPLICON_FAVORABLE_MAX:     int = 1000
# Amplicones > AMPLICON_FAVORABLE_MAX se consideran "menos favorables"

# Puntuaciones correspondientes a cada rango
PUNTUACION_AMPLICON_MUY_FAVORABLE: float = 100.0
PUNTUACION_AMPLICON_FAVORABLE:     float = 60.0
PUNTUACION_AMPLICON_DESFAVORABLE:  float = 20.0


# UMBRALES DE DIFERENCIA DE Tm (en °C)
TM_DIFERENCIA_OPTIMA:      float = 2.0   # ≤ 2°C → óptimo
TM_DIFERENCIA_ACEPTABLE:   float = 5.0   # entre 2 y 5°C → aceptable
# > 5°C → desfavorable

PUNTUACION_TM_OPTIMA:      float = 100.0
PUNTUACION_TM_ACEPTABLE:   float = 60.0
PUNTUACION_TM_DESFAVORABLE:float = 20.0


# UMBRALES DE SITIOS DE RESTRICCIÓN
RESTRICCION_IDEAL:       int = 0    # 0 sitios → ideal
RESTRICCION_ACEPTABLE:   int = 2    # 1–2 sitios → aceptable
# > 2 sitios → desfavorable

PUNTUACION_RESTRICCION_IDEAL:       float = 100.0
PUNTUACION_RESTRICCION_ACEPTABLE:   float = 50.0
PUNTUACION_RESTRICCION_DESFAVORABLE:float = 10.0


ARCHIVO_PRIMERS_DEFAULT:    str = "primers.tab"
ARCHIVO_RESTRICCION_DEFAULT:str = "res_enzimas.tab"

# ARCHIVOS DE SALIDA
SALIDA_PRIORIZACION: str = "priorizacion_genes.tsv"
SALIDA_REPORTE:      str = "reporte_priorizacion.txt"
SALIDA_BLAST_DIR:    str = "blast_resultados"


# TOP N genes para justificación automática en el reporte
TOP_N_GENES: int = 10