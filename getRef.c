/*
 *  @Developer:     Juan Camilo Peña Vahos
 *  @Last Revised:  26/06/2018
 *  @Description:   Obtener la referencia de un archivo
 *
*/
	


//LIBRERÍAS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

int main (int argc, char *argv[]) {	
    
	//ARGUMENTOS DE ENTRADA
	char* DATA;
	int L	=	0;
	//VARIABLES DEL PROCESO PARA SALIDA
	float lambda = 0;									//
	char	*MT;										//MATCHING TYPE
	char	*OT;										//OPERATION (MUTATION) TYPE

	//Obtener los datos suministrados en la linea de comando
	if(argc>1){
		for (int i = 1; i < argc; i++)	{
			if (strcmp(argv[i], "-DATA") == 0)	DATA	=	argv[i+1];
			if (strcmp(argv[i], "-L") == 0)		L	=	atoi(argv[i+1]);
		}
	}
	//A.GENERAR ALEATORIAMENTE LA POSICIÓN DE MAPEO
	//	1. SE DEBE TENER EN UN ARREGLO LA REFERENCIA
    //      *SABER LA CANTIDAD DE CARACTERES DE LA REFERENCIA
	char *Reference;
	int contChars	=	0;							
    int estado;
	char c;
	FILE *pf;

	pf = fopen(DATA,"r"); 
	if(pf!=NULL){								
		while((c = getc(pf)) != EOF){
            if(c == '>')    estado  =   0;		//EL ESTADO 0 SON LOS COMENTARIOS           
            switch(estado){
                case 0: if(c    ==  '\n')   estado  =   1;  break;
                case 1: if(c    !=  '\n')   contChars++;    break;
            }
		}
	}
    //printf("Numero de caracteres = %d\n", contChars);
	fclose(pf);
	Reference	=	malloc(contChars*sizeof(char*));
	//		*OBTENER LA REFERENCIA
	int j	=	0;
	pf	=	fopen(DATA,"r");
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

    //      	*GENERAR UNA POSICIÓN ALEATORIA EN EL RANGO [0,LengthRef-LengthRead]
    srand(time(0));
    int position = (rand() %((contChars-L) - 0 + 1)) + 0;
    //printf("%d\n", position);

    //      		*SACAR LA PORCIÓN DEL ARREGLO
	printf("L	=	%d;  POS	=	%d\n",L,position);
	char READ[L];
	for(int i = 0;	i < L;	i++){
		READ[i]	=	Reference[position+i];
	}
	printf("Read:	%s\n",READ);






	//OBTENER LA REFERENCIA DIRECTAMENTE DEL ARCHIVO
    /*char REF[L];
    //REF = malloc(L*sizeof(char*));                                //REFERENCIA
    int i   =   0;                                                  //Contador para la referencia
    int contAux =   0;  
    int upper   =   position+L;
    //printf("Limite inferior = %d, Limite Superior = %d\n",position,upper);                        
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
	fclose(pf);*/
	return 0;
}