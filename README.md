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
C | Valor entre [10 - 100], determina la cobertura en el muestreo de los Reads | uint8_t
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

### VARIABLES DEL PROCESO DE MUTACIÓN


TYPE | VARIABLE | DESCRIPCIÓN
-------- | ----------- | ----
uint32_t  | id | Identificador consecutivo para cada uno de los Reads
uint32_t  |	Pos | Posición de Matching respecto a la referencia
uint16_t  |	lendesc | Cantidad de errores total en el Read
char      | strand | Caractér con el sentido del matching
uint8_t*  | Oper | Arreglo con la operación por error
uint16_t* |	Cnt	| Arreglo con los contadores por cada uno de los tipos de mutación
uint32_t* |	Hist | Arreglo con el acumulador de contadores Cnt
uint16_t* | Offsets | Arreglo de offsets por cada error
uint16_t* |	OffRel | Arreglo de offsets pero relativos a la mutación
uint8_t*  | BaseRef | Arreglo con la base de la referencia (Read Referencia)
uint8_t*  | BaseRead | Arreglo con la base después de la mutación (Read Destino)

### DESCRIPCIÓN DEL ALGORITMO 

1. El archivo lee la referencia de extensión .fa, de ella elimina todos los comentarios y los saltos de línea y solo toma la secuencia que almacena en el arrelo **Reference**.
2. Usando una distribución uniforme, se selecciona un punto en el arreglo **Reference** a partir de dicho punto se toman **L+READ_BIAS** elementos.
3. Usando una distribución exponencial se cálcula la máxima cantida de errores **Kmax** para el **P0** que se determinó en la línea de comandos y usando una distribución uniforme, se determina la cantidad de errores que va a tener dicho read en el intervalo [0-Kmax]
4. La distribución uniforme se usa también para determinar cuál tipo de matching va a tener el Read, qué mutación se va a aplicar en cierto offset y finalmente,si la mutación lo necesita, que base va a reemplazar la base de la referencia.
5. El algoritmo son dos ciclos anidados, el primero genera Reads completos, y el segundo aplica las mutaciones en los offsets determinados a cada Read.
6. En caso de que **lendesc = 0** se dice que el matching fué perfecto y se imprime directamente el Read.


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

## COMO COMPILAR 
- [x] gcc -o ARF MainARF.c Matching.c Mutation.c  Stats.c FilesUtils.c -lm
- [x] ./ARF -DATA lambda_virus.fa -I -Q -L 1024 -B 200000 -C 10 -P0 0.02

_NOTA:_ Recuerde que los parámtros C y P0 se van a variar en las pruebas    

## _ACLARACIÓN SOBRE DISTRIBUCIÓN_
En este momento la distribución de los offsets, es uniforme, unicamente, en el intervalo [0 , (L-280)]. Funciona solo para pruebas con L = 1024
Esta decisión viene dada por errores técnicos en los excedentes al principio o al final del read al aplicar las mutaciones en los distintos casos de matchings.

_NOTA:_ No modifique el branch master, si desea agregar cambios o enviar comentarios porfavor crear un nuevo branch y deje un comentario con los cambios que va a hacer, si cree que es necesario modificar el master, enviar un correo a **milovahos@gmail.com**

### DEVELOPERS:
_Juan Camilo Peña Vahos_ @MiloVahos96,
_Aníbal Guerra_,
_Sebastian Isaza Ramírez_,