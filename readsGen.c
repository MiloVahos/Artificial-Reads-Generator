/*
 *  @Developer:     Juan Camilo Peña Vahos
 *  @Last Revised:  09/07/2018
 *  @Description:   Generador artificial de reads para simular casos de prueba del compresor.
 *
*/
/*
 *	Todo: 	     	Inicializar Read usando memoria dinámica
*/
/*								ESTRUCTURA DE UN READ
 * LINEA 1: (IDENTIFICADOR) Siempre empieza con @
 * LINEA 2:	(SECUENCIA)		AGNTAGNTAGNT
 * LINEA 3:	(COMENTARIO)	Siempre empieza con +
 * LINEA 4:	(QUALITY SCORE)	Es un valor que entrega la máquina
 * 
 * NOTAS:
 * 		->EL QS SE GENERA CON UN SOLO VALOR QUE SE REPETIRÁ HASTA ALCANZAR LA MISMA
 * 		  LONGITUD QUE LA SECUENCIA
*/
/*                      PASOS DEL GENERADOR ARTIFICIAL DE READS
 *	1. Abrir y leer los diferentes tipos de archivos fasta y los argumentos de entrada
 * 	
 * 	ENTRADA: 
 * 		NOTA: 		SE DEBE RESPETAR ESTA NOTACIÓN
 * 
 * 		DATA	->	NOMBRE ARCHIVO DE EXTENSIÓN FASTA
 * 		I	  	->	IDENTIFICADOR PARA TODOS LOS READS, LOS DOS ÚLTIMOS NÚMEROS CAMBIAN DE 
 * 					FORMA ALEATORIA,TAMBIÉN ES EL ID DEL ARCHIVO .FASTQ
 *		Q		->	QUALITY SCORE PARA TODOS LOS ARCHIVOS
 * 		L		->	LONGITUD DE LOS READS
 * 		C		->	COVERAGE: ES INFORMACIÓN DE IMPORTANCIA BIOLÓGICA, NO PARA NOSOTROS
 * 		B		->	BASE DE READS: EL TOTAL DE READS A GENERAR SERÁ BxC 
 * 		E		-> 	CANTIDAD DE ERRORES DE CADA READ, DEBE SER MENOR QUE LA LONGITUD
 * 		P		->	ES UNA BANDERA, SI ESTÁ ACTIVADA, ENTONCES SOLO SE CONSIDERAN DOS TIPOS
 * 					DE MATCHING, FORWARD AND REVERSE EQUIPROBABLES DE SUCEDER (OPCIONAL)
 * 
 * 	2.	PROCESO DE GENERACIÓN DE READS
 *		A. GENERAR ALEATORIAMENTE LA POSICIÓN DE MAPEO
		B. CALCULAR EL MATCHING
*/
/*
 * 	SALIDA DEL PROGRAMA:
 * 		1. ARCHIVO .fastq CON BXC READS
 * 		2. ARCHIVO .fastqseq CON LAS SECUENCIAS DE CADA READ,UNA POR LÍNEA
 * 		3. ARCHIVO .aling CONTIENE POR READ LOS DATOS DEL MATCHING DEL READ Y TODOS LOS DATOS
 * 						  DE LOS ERRORES
 * 		4. ARCHIVO .meta CONTIENE LOS METADATOS DEL EXPERIMENTO,LOS QUE SE PASAN POR LA CONSOLA
 * 					     Y ALGUNOS QUE SE GENERAN EN EL TRANSCURSO DE LA EJECUCACIÓN COMO LOS
 * 						 NOMBRES DE LOS ARCHIVOS ANTERIORES
*/ 

//LIBRERÍAS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <inttypes.h>

//DECLARACIÓN DE CONSTANTES
#define MATCHING_TYPES 4
#define EXP 2.71828

//DECLARACIÓN DE FUNCIONES
void 	ReverseRead(char*,long);
void 	ComplementRead(char*,long);
int  	BusqBin_Rul(double[], int, double);
double 	CalculaTotal(int,double[]);
double 	LanzarDado();
long int contChars(char*);
void 	getReference(char*,char*);
void 	getRead(char*,char*,long int,int);

//DECLARACIÓN DE VARIABLES GLOBALES

int main (int argc, char *argv[]) {	
    
	//ARGUMENTOS DE ENTRADA
	char 	*DATA;
	char 	*I;
	char 	*Q;
	int 	L	=	0;
	int		C	=	0;
	int 	B	=	0;
	int 	E	=	0;

	//VARIABLES DEL PROCESO PARA SALIDA
	float 	lambda = 0;									//
	char	*MT;										//MATCHING TYPE
	char	*OT;										//OPERATION (MUTATION) TYPE

	//VARIABLES DEL PROCESO PARA OPERAR
	char		*Reference	=	NULL;					//REFERENCIA OBTENIDA DEL ARCHIVO FASTA
	char 		*Name		=	NULL;					//NOMBRE DE LOS ARCHIVOS
	int			TotalReads;								//TotalReads	=	BxC
	long int 	TotalChars;								//Total de caracteres en la referencia

	//Obtener los datos suministrados en la linea de comando
	if(argc>1){
		for (int i = 1; i < argc; i++)	{
			if (strcmp(argv[i], "-DATA") == 0)	
				DATA	=	argv[i+1];
			if (strcmp(argv[i], "-I") == 0)	
				I	=	argv[i+1];
			if (strcmp(argv[i], "-Q") == 0)	
				Q	=	argv[i+1];	
			if (strcmp(argv[i], "-L") == 0)
				L	=	atoi(argv[i+1]);
			if (strcmp(argv[i], "-C") == 0)	
				C	=	atoi(argv[i+1]);
			if (strcmp(argv[i], "-B") == 0)	
				B	=	atoi(argv[i+1]);
			if (strcmp(argv[i], "-E") == 0)	
				E	=	atoi(argv[i+1]);
		}
	}

	//OBTENER EL NOMBRE DE LA SECUENCIA Y CREAR EL NOMBRE DE LOS ARCHIVOS DE SALIDA
	Name	=	malloc(50*sizeof(char));
	char *dot;
	dot		=	strrchr(DATA,'.');
	strncpy(Name,DATA,dot - DATA);

	//A. GENERAR ALEATORIAMENTE LA POSICIÓN DE MAPEO
    //	1.SABER LA CANTIDAD DE CARACTERES DE LA REFERENCIA
	TotalChars	=	contChars(DATA);
	//		2.OBTENER LA REFERENCIA
	Reference	=	malloc(TotalChars*sizeof(char*));
	getReference(DATA,Reference);
	TotalReads	=	B*C;				//Número total de Reads a generar

	for(int ReadsCicle = 0;	ReadsCicle<TotalReads; ReadsCicle++){
		//			3.GENERAR UNA POSICIÓN ALEATORIA EN EL RANGO [0,LengthRef-LengthRead]
		srand(time(0));
		long int position = (rand() %((TotalChars-L) - 0 + 1)) + 0;
		//				4.SACAR LA PORCIÓN DEL ARREGLO
		char *Read	=	NULL;
		Read	=	malloc(L*sizeof(char*));
		getRead(Reference,Read,position,L);
		printf("Read:	%s\n",Read);
		/*//B.CALCULAR EL MATCHING
		
		//ORDEN DE LOS ARREGLOS 
		//FORWARD(F), REVERSE(R), COMPLEMENT(C), REVERSE COMPLEMENT(E)
		double MatTypeStats[MATCHING_TYPES]   =   {0.4,0.4,0.1,0.1};  
		double Acum_fun[MATCHING_TYPES];

		double total_adapt = CalculaTotal(MATCHING_TYPES,MatTypeStats);
		Acum_fun[0]= MatTypeStats[0] / total_adapt;
		for (int i=1; i<MATCHING_TYPES; i++){
			Acum_fun[i]= (MatTypeStats[i] / total_adapt) + Acum_fun[i-1];
		}

		//LANZAR EL DADO
		double dado = LanzarDado();
		printf("Resultado del dado = %lf\n", dado);

		//GIRAR LA RULETA
		int MatTypeSel  =   BusqBin_Rul(Acum_fun,MATCHING_TYPES,dado);
		printf("Selección de la ruleta = %d\n", MatTypeSel);

		switch(MatTypeSel){
			case 0:                     //FORWARD MATCH
				MT  =   "F";
				printf("F: %s\n",READ);
			break;
			case 1:                     //REVERSE MATCH
				MT  =   "R";
				ReverseRead(READ,L);
				printf("R: %s\n",READ);
			break;
			case 2:                     //COMPLEMENT MATCH
				MT  =   "C";
				ComplementRead(READ,L);
				printf("C: %s\n",READ);
			break;
			case 3:                     //REVERSE COMPLEMENT MATCH
				MT  =   "E";
				ComplementRead(READ,L);
				ReverseRead(READ,L);
				printf("E: %s\n",READ);                
			break;
			default:    printf ("**Error in the matching selection, wrong input base**");
		}*/

		free(Read);
	}

    


	free(Name);
	free(Reference);
	return 0;
}

//FUNCIONES

//contChars: Cuenta el número de caractéres que tiene el archivo de referencia
//@param: char *FileName : Nombre del archivo
long int contChars(char *FileName){
	long int cont	=	0;							
    int estado;
	char c;
	FILE *pf;

	pf = fopen(FileName,"r"); 
	if(pf!=NULL){								
		while((c = getc(pf)) != EOF){
            if(c == '>')    estado  =   0;		//EL ESTADO 0 SON LOS COMENTARIOS           
            switch(estado){
                case 0: if(c    ==  '\n')   estado  =   1;  break;
                case 1: if(c    !=  '\n')   cont++;    break;
            }
		}
	}

	return cont;
}

//getReference: Obtiene la referencia de los archivos FASTA y la guarda en un arreglo
//@param: char *FileName : Nombre del archivo
//@param: char *Reference: Arreglo en el que se guarda la referencia
void getReference(char *FileName, char *Reference){
	int estado;
	char c;
	int j	=	0;
	FILE *pf;
	pf	=	fopen(FileName,"r");
	if(pf	!=	NULL){
		while((c = getc(pf)) != EOF){
            if(c == '>')    estado  =   0;                          
            switch(estado){
                case 0: if(c == '\n')	estado  =   1;	break;
                case 1: if(c    !=  '\n'){   
                    Reference[j]	=	c;
					j++;
                }   
                break;
            }
		}
	}
	fclose(pf);
}

//getRead:	Obtiene el read de la referencia
//@param:	char *Reference		:Arreglo con la referencia
//@param:	char *READ			:Arreglo donde se va a almacenar el Read
//@param:	long int position	:Posición desde la cual se va a tomar el Read
//@param:	int L				:Longitud del Read
void getRead(char* Reference,char* READ,long int position,int L){
	for(int i = 0;	i < L;	i++){
		switch(Reference[position+i]){
			case 'a':	READ[i] = 'A'; break;
			case 'A':	READ[i]	= 'A'; break;
			case 'c':	READ[i] = 'C'; break;
			case 'C':	READ[i]	= 'C'; break;
			case 'g':	READ[i] = 'G'; break;
			case 'G':	READ[i]	= 'G'; break;
			case 't':	READ[i] = 'T'; break;
			case 'T':	READ[i]	= 'T'; break;
			default:	READ[i]	= 'N';
		}
	}
}

//ReverseRead:	Implementa el inversor de reads
//@param:	char *Read	:	Arreglo con el Read
//@param:	long length	:	Longitud del Read
void ReverseRead(char *Read, long length){
	char aux;
    int resto = length%2;
	for (int i=0; i<length/2;i++){
        //printf("READI = %c , READLEN = %c, i =  %d, resto = %d\n",Read[i],Read[length-i-1], i, length-i-1);
		aux=Read[length-i-1];
		Read[length-i-1]=Read[i];
		Read[i]=aux;
	}
}

//ComplementRead:	Implementa el complementador de reads (T,A)(C,G)
//@param:	char *Read	:	Arreglo con el Read
//@param:	long length	:	Longitud del Read
void ComplementRead(char *Read, long length){
	char Compl;
	int i;
	for (i=0; i<length;i++){
		switch(Read[i]){
		    case 'A': Compl='T'; break;
		    case 'a': Compl='T'; break;
		    case 'C': Compl='G'; break;
		    case 'c': Compl='G'; break;
		    case 'G': Compl='C'; break;
		    case 'g': Compl='C'; break;
		    case 'T': Compl='A'; break;
		    case 't': Compl='A'; break;
		    case 'N': Compl='N'; break; 
		    case 'n': Compl='N' ;break;
		    default: printf ("**Error building the complement, wrong input base**");
		}

		Read[i]=Compl;
	}

}

//BusqBin_Rul: función que utiliza el procedimiento de la Busq. Binaria
//             para ubicar el elemento seleccionado a través del mecanismo
//	            de la ruleta
//@param: double prob[] : vector de probabilidades acumuladas
//@param: n :   número de elementos en prob[]
//@param: x :   resultado de lanzar los dados (es decir, usar rand para generar un número)
int BusqBin_Rul(double prob[], int n, double x){

	int primero,ultimo,central;
	short encontrado;
    primero     =   0;
    ultimo      =   n-1;
    encontrado  =   0;
	while ((primero <= ultimo) && (encontrado==0)){
   	    central = (primero + ultimo)/2;
        if      (x  ==   prob[central]) encontrado = 1;
        else if (x  >    prob[central]) primero = central + 1;
        else                            ultimo = central - 1;
        //printf("primer %d\n",primero);
        //printf("ultimo %d\n",ultimo);
   }
   return(primero);
}

//CalculaTotal: Calculo de la Sumatoria de las Func. Aadaptación
//              de cada Individuo de la Poblacion.
//@param: a[] : vector de probabilidades de cada individuo, sin acumular
//@param: n :   número de elementos en a[]
double CalculaTotal(int n, double a[]){

	int i;
	double aux=0;
    for (i=0;i<n;i++) {	aux = aux + a[i];}
	return aux;
}

//LanzarDado: devuelve un valor entre 0 y 1 aleatorio
double LanzarDado(){
    srand(time(0));
    double dado = 0+(1-0)*rand()/((double)RAND_MAX);
    return dado;
}