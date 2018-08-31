# ARTIFICIAL READS GENERATOR DOCS

## ESTRUCTURA DE UN READ
 * LINEA 1: (IDENTIFICADOR) Siempre empieza con @
 * LINEA 2:	(SECUENCIA)		AGNTAGNTAGNT
 * LINEA 3:	(COMENTARIO)	Siempre empieza con +
 * LINEA 4:	(QUALITY SCORE)	Es un valor que entrega la máquina

_NOTA: SE DEBE RESPETAR ESTA NOTACIÓN_

## TABLAS DE VARIABLES

### CONSTANTES

1. READ_BIAS:   Número de elementos extras que se agregan al Read para tener margen de overflow
2. ERROR_PER:   Porcentaje de error válido en cada Read
3. NAMES_SIZE:  Tamaño del char que contiene los nombres	
4. READID_SIZE: Tamaño de los IDs
5. READQ_SIZE:  Tamaño de los QSs

### VARIABLES DE ENTRADA

_NOTA: SE INGRESAN POR LA CONSOLA DE COMANDO, LA SINTAXIS ES -NombreVariable valorVariable, -I Y -Q SON BANDERAS_

VARIABLE | DESCRIPCIÓN | TYPE
-------- | ----------- | ----
DATA | Nombre del archivo de extensión fasta | char*
I | Bandera del identificador de los reads, es un valor fijo, cambia el id | char*
Q | Bandera de los quality score de los reads, es un valor fijo, no cambia | char*
L | Longitud de la secuencia de los reads | uint16_t [0-1024]
C | Valor entre [10 - 100], determina la cantidad de Reads | uint8_t
B | Valor fijo en 200000, B*C = Cantidad de reads a generar | uint32_t
P0 | Es el valor de ajuste de la distribución exponencial | double

### VARIABLES DEL PROCESO

TYPE | VARIABLE | DESCRIPCIÓN
-------- | ----------- | ----
char* | Reference | Contiene todos las bases del archivo .fa, no contiene los comentarios
char* | Read | Segmento de la referencia que se va a modificar y se imprimirá como la secuencia
char* | RefName | Nombre del archivo de referencia, sin la extensión
char* | RefFastq | Nombre del archivo fastq
char* | RefFastqseq | Nombre del archivo fastqseq
char* | RefAlign | Nombre del archivo align
char* | RefMeta | Nombre del archivo meta
uint32_t | TotalReads | Total de reads que se van a generar, se calcula como B*C (Mirar Pruebas)
uint64_t | TotalChars | Total de caracteres en la referencia (Número de bases en la secuencia)
FILE | *FASTQ, *FASTQSEQ | Punteros a los archivos de salida
FILE | *ALIGN,	*META | Punteros a los archivos de salida

### CONSIDERACIONES DE PRUEBA

1. B es un valor fijo esn 200000
2. L es un valor fijo en 1024
3. Sea K la cantidad máxima de errores de errores, que depende de P0, se hacen pruebas para
    K= 230, 170, 110, 80, 60, 40, 28, 24, 20, Todos entre [0-19]

## SALIDA DEL PROGRAMA
 1. ARCHIVO .fastq:     CON BXC READS
 2. ARCHIVO .fastqseq:  CON LAS SECUENCIAS DE CADA READ,UNA POR LÍNEA
 3. ARCHIVO .aling:     CONTIENE POR READ LOS DATOS DEL MATCHING DEL READ Y TODOS LOS DATOS 
    DE LOS ERRORES
 4. ARCHIVO .meta:      CONTIENE LOS METADATOS DEL EXPERIMENTO,LOS QUE SE PASAN POR LA CONSOLA
                        Y ALGUNOS QUE SE GENERAN EN EL TRANSCURSO DE LA EJECUCACIÓN COMO LOS 
                        NOMBRES DE LOS ARCHIVOS ANTERIORES

## _ACLARACIÓN SOBRE DISTRIBUCIÓN_
En este momento la distribución de los offsets, es uniforme, unicamente, en el intervalo [0 - L-200]. Funciona solo para pruebas con L = 1024
Esta decisión viene dada por errores técnicos en los excedentes al principio o al final del read al aplicar las mutaciones en los distintos casos
de matchings.

### DEVELOPERS:
_Juan Camilo Peña Vahos_ @MiloVahos96,
_Aníbal Guerra_,
_Sebastian Isaza Ramírez_,
