/*
 *  @Developer:     Juan Camilo Peña Vahos
 *  @Last Revised:  26/06/2018
 *  @Description:   Generador artificial de reads para simular casos de prueba del compresor.
 *
*/
/*
 *	Todo: 	        .
*/
/*
 *	ULTIMA LABOR:		 
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


//LIBRERÍAS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

//DECLARACIÓN DE FUNCIONES

//DECLARACIÓN DE VARIABLES GLOBALES


int main (int argc, char *argv[]) {	
    
	//ARGUMENTOS DE ENTRADA
	char* DATA;
	char* I;
	char* Q;
	int L	=	0;
	int	C	=	0;
	int B	=	0;
	int E	=	0;

	//VARIABLES DEL PROCESO PARA SALIDA
	float lambda = 0;									//
	char	*MT;										//MATCHING TYPE
	char	*OT;										//OPERATION (MUTATION) TYPE

	//VARIABLES DEL PROCESO PARA OPERAR
	char	*Reference	=	NULL;						//REFERENCIA OBTENIDA DEL ARCHIVO FASTA
	

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
	//C1

	//A.GENERAR ALEATORIAMENTE LA POSICIÓN DE MAPEO
	//	1. SE DEBE TENER EN UN ARREGLO LA REFERENCIA
	int contChars	=	0;							
    int estado;
	char c;
	FILE *pf;

	pf = fopen(DATA,"r"); 
	if(pf!=NULL){								
		while((c = getc(pf)) != EOF){
            if(c == '>')    estado  =   0;                          //EL ESTADO 0 SON LOS COMENTARIOS           
            switch(estado){
                case 0: if(c    ==  '\n')   estado  =   1;  break;
                case 1: if(c    !=  '\n')   contChars++;    break;
            }
		}
	}
	fclose(pf);
    //      *GENERAR UNA POSICIÓN ALEATORIA EN EL RANGO [0,LengthRef-LengthRead]
    srand(time(0));
    int position = (rand() %((contChars-L) - 0 + 1)) + 0;
    //      *SACAR LA PORCIÓN DEL ARCHIVO
    char REF[L];
    //REF = malloc(L*sizeof(char*));                                //REFERENCIA
    int i   =   0;                                                  //Contador para la referencia
    int contAux =   0;  
    int upper   =   position+L;                    
    pf = fopen(DATA,"r"); 
	if(pf!=NULL){								
		while((c = getc(pf)) != EOF){
            if(c == '>')    estado  =   0;                          //EL ESTADO 0 SON LOS COMENTARIOS           
            switch(estado){
                case 0: 
                    if(c    ==  '\n'){
                        estado  =   1;
                    }
                break;
                case 1: if(c    !=  '\n'){   
                    contAux++;
                }   
                break;
            }
            if((contAux  >= position)&&(contAux<upper)){
                REF[i]  =   c;
                i++;
            }
		}
	}
    printf("REFERENCIA: %s\n",REF);
	fclose(pf);
	return 0;
}
