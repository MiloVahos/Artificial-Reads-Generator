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

## TABLA DE K vs LAMBDA
MinK    =   0
MaxK    =   256

K | LAMBDA | PROBABILIDAD
-- | ------ | ------------
0 | 1.000000 | 1.000000
1 | 0.657500 | 0.824813
2 | 0.507500 | 0.889723
3 | 0.422500 | 0.999846
4 | 0.362500 | 0.942308
5 | 0.320000 | 0.954363
6 | 0.287500 | 0.962613
7 | 0.262500 | 0.969900
8 | 0.242500 | 0.976179
9 | 0.225000 | 0.979740
10 | 0.210000 | 0.982297
11 | 0.197500 | 0.985096
12 | 0.185000 | 0.984468
13 | 0.175000 | 0.985645
14 | 0.167500 | 0.989227
15 | 0.160000 | 0.990909
16 | 0.152500 | 0.990965
17 | 0.147500 | 0.999952
18 | 0.140000 | 0.991868
19 | 0.135000 | 0.993112
20 | 0.130000 | 0.993440
21 | 0.125000 | 0.992944
22 | 0.122500 | 0.996735
23 | 0.117500 | 0.994894
24 | 0.115000 | 0.997501
25 | 0.110000 | 0.994528
26 | 0.107500 | 0.996189
27 | 0.105000 | 0.997426
28 | 0.102500 | 0.998283
29 | 0.100000 | 0.998794
30 | 0.097500 | 0.998990
31 | 0.095000 | 0.998892
32 | 0.092500 | 0.998521
33 | 0.090000 | 0.997890
34 | 0.087500 | 0.997009
35 | 0.087500 | 0.997009
36 | 0.085000 | 0.999948
37 | 0.082500 | 0.998478
38 | 0.080000 | 0.996783
39 | 0.080000 | 0.996783
40 | 0.077500 | 0.998370
41 | 0.077500 | 0.998370
42 | 0.075000 | 0.999379
43 | 0.075000 | 0.999379
44 | 0.072500 | 0.999886
45 | 0.072500 | 0.999886
46 | 0.070000 | 0.999948
47 | 0.070000 | 0.999948
48 | 0.067500 | 0.999608
49 | 0.067500 | 0.999608
50 | 0.065000 | 0.998895
51 | 0.065000 | 0.998895
52 | 0.065000 | 0.999900
53 | 0.065000 | 0.999900
54 | 0.065000 | 0.999900
55 | 0.060000 | 0.998538
56 | 0.060000 | 0.998538
57 | 0.060000 | 0.998538
58 | 0.057500 | 0.998667
59 | 0.057500 | 0.998667
60 | 0.057500 | 0.998667
61 | 0.057500 | 0.999908
62 | 0.057500 | 0.999908
63 | 0.057500 | 0.999908
64 | 0.057500 | 0.999908
65 | 0.052500 | 0.999042
66 | 0.052500 | 0.999042
67 | 0.052500 | 0.999042
68 | 0.052500 | 0.999042
69 | 0.050000 | 0.999129
70 | 0.050000 | 0.999129
71 | 0.050000 | 0.999129
72 | 0.050000 | 0.999129
73 | 0.050000 | 0.999862
74 | 0.050000 | 0.999862
75 | 0.050000 | 0.999862
76 | 0.050000 | 0.999862
77 | 0.050000 | 0.999862
78 | 0.047500 | 0.999918
79 | 0.047500 | 0.999918
80 | 0.047500 | 0.999918
81 | 0.047500 | 0.999918
82 | 0.047500 | 0.999918
83 | 0.047500 | 0.999918
84 | 0.042500 | 0.999374
85 | 0.042500 | 0.999374
86 | 0.042500 | 0.999374
87 | 0.042500 | 0.999374
88 | 0.042500 | 0.999374
89 | 0.042500 | 0.999374
90 | 0.040000 | 0.999154
91 | 0.040000 | 0.999154
92 | 0.040000 | 0.999154
93 | 0.040000 | 0.999154
94 | 0.040000 | 0.999154
95 | 0.040000 | 0.999154
96 | 0.040000 | 0.999154
97 | 0.040000 | 0.999893
98 | 0.040000 | 0.999893
99 | 0.040000 | 0.999893
100 | 0.040000 | 0.999893
101 | 0.040000 | 0.999893
102 | 0.040000 | 0.999893
103 | 0.040000 | 0.999893
104 | 0.040000 | 0.999893
105 | 0.040000 | 0.999893
106 | 0.035000 | 0.999760
107 | 0.035000 | 0.999760
108 | 0.035000 | 0.999760
109 | 0.035000 | 0.999760
110 | 0.035000 | 0.999760
111 | 0.035000 | 0.999760
112 | 0.035000 | 0.999760
113 | 0.035000 | 0.999760
114 | 0.035000 | 0.999760
115 | 0.032500 | 0.999446
116 | 0.032500 | 0.999446
117 | 0.032500 | 0.999446
118 | 0.032500 | 0.999446
119 | 0.032500 | 0.999446
120 | 0.032500 | 0.999446
121 | 0.032500 | 0.999446
122 | 0.032500 | 0.999446
123 | 0.032500 | 0.999446
124 | 0.032500 | 0.999446
125 | 0.032500 | 0.999446
126 | 0.032500 | 0.999952
127 | 0.032500 | 0.999952
128 | 0.032500 | 0.999952
129 | 0.032500 | 0.999952
130 | 0.032500 | 0.999952
131 | 0.032500 | 0.999952
132 | 0.032500 | 0.999952
133 | 0.032500 | 0.999952
134 | 0.032500 | 0.999952
135 | 0.032500 | 0.999952
136 | 0.032500 | 0.999952
137 | 0.032500 | 0.999952
138 | 0.032500 | 0.999952
139 | 0.030000 | 0.999854
140 | 0.030000 | 0.999854
141 | 0.030000 | 0.999854
142 | 0.030000 | 0.999854
143 | 0.030000 | 0.999854
144 | 0.030000 | 0.999854
145 | 0.030000 | 0.999854
146 | 0.030000 | 0.999854
147 | 0.030000 | 0.999854
148 | 0.030000 | 0.999854
149 | 0.030000 | 0.999854
150 | 0.030000 | 0.999854
151 | 0.030000 | 0.999854
152 | 0.030000 | 0.999854
153 | 0.030000 | 0.999854
154 | 0.030000 | 0.999854
155 | 0.027500 | 0.999919
156 | 0.027500 | 0.999919
157 | 0.027500 | 0.999919
158 | 0.027500 | 0.999919
159 | 0.027500 | 0.999919
160 | 0.027500 | 0.999919
161 | 0.027500 | 0.999919
162 | 0.027500 | 0.999919
163 | 0.027500 | 0.999919
164 | 0.027500 | 0.999919
165 | 0.027500 | 0.999919
166 | 0.027500 | 0.999919
167 | 0.027500 | 0.999919
168 | 0.027500 | 0.999919
169 | 0.027500 | 0.999919
170 | 0.027500 | 0.999919
171 | 0.027500 | 0.999919
172 | 0.027500 | 0.999919
173 | 0.027500 | 0.999919
174 | 0.025000 | 0.999807
175 | 0.025000 | 0.999807
176 | 0.025000 | 0.999807
177 | 0.025000 | 0.999807
178 | 0.025000 | 0.999807
179 | 0.025000 | 0.999807
180 | 0.025000 | 0.999807
181 | 0.025000 | 0.999807
182 | 0.025000 | 0.999807
183 | 0.025000 | 0.999807
184 | 0.025000 | 0.999807
185 | 0.025000 | 0.999807
186 | 0.025000 | 0.999807
187 | 0.025000 | 0.999807
188 | 0.025000 | 0.999807
189 | 0.025000 | 0.999807
190 | 0.025000 | 0.999807
191 | 0.025000 | 0.999807
192 | 0.025000 | 0.999807
193 | 0.025000 | 0.999807
194 | 0.025000 | 0.999807
195 | 0.025000 | 0.999807
196 | 0.025000 | 0.999807
197 | 0.025000 | 0.999807
198 | 0.022500 | 0.999803
199 | 0.022500 | 0.999803
200 | 0.022500 | 0.999803
201 | 0.022500 | 0.999803
202 | 0.022500 | 0.999803
203 | 0.022500 | 0.999803
204 | 0.022500 | 0.999803
205 | 0.022500 | 0.999803
206 | 0.022500 | 0.999803
207 | 0.022500 | 0.999803
208 | 0.022500 | 0.999803
209 | 0.022500 | 0.999803
210 | 0.022500 | 0.999803
211 | 0.022500 | 0.999803
212 | 0.022500 | 0.999803
213 | 0.022500 | 0.999803
214 | 0.022500 | 0.999803
215 | 0.022500 | 0.999803
216 | 0.022500 | 0.999803
217 | 0.022500 | 0.999803
218 | 0.022500 | 0.999803
219 | 0.022500 | 0.999803
220 | 0.022500 | 0.999803
221 | 0.022500 | 0.999803
222 | 0.022500 | 0.999803
223 | 0.022500 | 0.999803
224 | 0.022500 | 0.999803
225 | 0.022500 | 0.999803
226 | 0.022500 | 0.999803
227 | 0.022500 | 0.999803
228 | 0.022500 | 0.999803
229 | 0.020000 | 0.999881
230 | 0.020000 | 0.999881
231 | 0.020000 | 0.999881
232 | 0.020000 | 0.999881
233 | 0.020000 | 0.999881
234 | 0.020000 | 0.999881
235 | 0.020000 | 0.999881
236 | 0.020000 | 0.999881
237 | 0.020000 | 0.999881
238 | 0.020000 | 0.999881
239 | 0.020000 | 0.999881
240 | 0.020000 | 0.999881
241 | 0.020000 | 0.999881
242 | 0.020000 | 0.999881
243 | 0.020000 | 0.999881
244 | 0.020000 | 0.999881
245 | 0.020000 | 0.999881
246 | 0.020000 | 0.999881
247 | 0.020000 | 0.999881
248 | 0.020000 | 0.999881
249 | 0.020000 | 0.999881
250 | 0.020000 | 0.999881
251 | 0.020000 | 0.999881
252 | 0.020000 | 0.999881
253 | 0.020000 | 0.999881
254 | 0.020000 | 0.999881
255 | 0.020000 | 0.999881
256 | 0.020000 | 0.999881


### DEVELOPERS:
_Juan Camilo Peña Vahos_ @MiloVahos96,
_Aníbal Guerra_,
_Sebastian Isaza Ramírez_,