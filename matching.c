/*
 *  @Developer:     Juan Camilo Peña Vahos
 *  @Last Revised:  27/06/2018
 *  @Description:   Algoritmo que determina el tipo de matching del read, de acuerdo con
 *                  una probabilidad.
*/	


//LIBRERÍAS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

//DECLARACIÓN DE FUNCIONES
void ReverseRead(char*,long);
void ComplementRead(char*,long);

int main (int argc, char *argv[]) {	
    
	//ARGUMENTOS DE ENTRADA
	int L	=   31;

	//VARIABLES DEL PROCESO PARA SALIDA								
	char *MT;										//MATCHING TYPE

    //PARA EL PROBAR EL MATCHING, SE ASUME LA SIGUIENTE CADENA.
    int MatType;
    char Read[L];
    
    //Read = (char*) malloc(L*sizeof(char));
    if (Read==NULL) {
        printf("There is no space in memory\n");
        exit (1);
    }
    strcpy(Read,"GGGCGGCGACCTCGCGGGTTTTCGCTATTTA");

    for(int i = 0; i<4; i++){
        switch(i){
            case 0:                     //DIRECT MATCH
                MT  =   "F";
                printf("DM: %s\n",Read);
            break;
            case 1:                     //REVERSE MATCH
                MT  =   "R";
                ReverseRead(Read,L);
                printf("RM: %s\n",Read);
            break;
            case 2:                     //COMPLEMENT MATCH
                MT  =   "C";
                ComplementRead(Read,L);
                printf("CM: %s\n",Read);
            break;
            case 3:                     //REVERSE COMPLEMENT MATCH
                MT  =   "E";
                ComplementRead(Read,L);
                ReverseRead(Read,L);
                printf("RCM: %s\n",Read);                
            break;
            default:    printf ("**Error in the matching selection, wrong input base**");
        }
        strcpy(Read,"GGGCGGCGACCTCGCGGGTTTTCGCTATTTA");
    }
}

//Implementa el Inversor de reads 
void ReverseRead(char *Read, long length){
	char aux;
    int resto = length%2;
    printf("%d\n",resto);
	for (int i=0; i<length/2;i++){
        //printf("READI = %c , READLEN = %c, i =  %d, resto = %d\n",Read[i],Read[length-i-1], i, length-i-1);
		aux=Read[length-i-1];
		Read[length-i-1]=Read[i];
		Read[i]=aux;
	}
};

//Implementa el complementador de reads (T,A)(C,G)
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

};