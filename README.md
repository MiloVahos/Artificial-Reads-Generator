#DOCUMENTACIÓN DEL GENERADOR ARITIFICIAL DE READS

##ESTRUCTURA DE UN READ
 * LINEA 1: (IDENTIFICADOR) Siempre empieza con @
 * LINEA 2:	(SECUENCIA)		AGNTAGNTAGNT
 * LINEA 3:	(COMENTARIO)	Siempre empieza con +
 * LINEA 4:	(QUALITY SCORE)	Es un valor que entrega la máquina

##ENTRADA DEL PROGRAMA
First Header | Second Header
------------ | -------------
Content from cell 1 | Content from cell 2
Content in the first column | Content in the second column

_NOTA: SE DEBE RESPETAR ESTA NOTACIÓN_
VARIABLE | DESCRIPCIÓN | TYPE
-------- | ----------- | ----
DATA | Nombre del archivo de extensión fasta | char*
I | Identificador de los reads, es un valor fijo, cambia el id | char*
Q | Quality score de los reads, es un valor fijo, no cambia | char*
L | Longitud de la secuencia de los reads | uint16_t [0-1024]
C | Valor entre [10 - 100], determina la cantidad de Reads | uint8_t
B | Valor fijo en 200000, B*C = Cantidad de reads a generar | uint32_t
P0 | Es el valor de ajuste de la distribución exponencial | double



##SALIDA DEL PROGRAMA
 1. ARCHIVO .fastq:     CON BXC READS
 2. ARCHIVO .fastqseq:  CON LAS SECUENCIAS DE CADA READ,UNA POR LÍNEA
 3. ARCHIVO .aling:     CONTIENE POR READ LOS DATOS DEL MATCHING DEL READ Y TODOS LOS DATOS 
    DE LOS ERRORES
 4. ARCHIVO .meta:      CONTIENE LOS METADATOS DEL EXPERIMENTO,LOS QUE SE PASAN POR LA CONSOLA
                        Y ALGUNOS QUE SE GENERAN EN EL TRANSCURSO DE LA EJECUCACIÓN COMO LOS 
                        NOMBRES DE LOS ARCHIVOS ANTERIORES