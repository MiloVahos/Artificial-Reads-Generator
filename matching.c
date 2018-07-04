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
int BusqBin_Rul(double[], int, double);

int main (int argc, char *argv[]) {	
    
	//ARGUMENTOS DE ENTRADA
	int L	=   30;

	//VARIABLES DEL PROCESO PARA SALIDA								
	char *MT;										//MATCHING TYPE

    //PARA EL PROBAR EL MATCHING, SE ASUME LA SIGUIENTE CADENA.
    int MatType;
    char Read[L];
    double MatTypeStats[]   =   {0.4,0.8,0.9,1};  //FORWARD, REVERSE, COMPLEMENT, REVERSE COMPLEMENT
                                                  //PROBABILITIES:  0.4 - 0.4 - 0.1 - 0.1

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
                printf("F: %s\n",Read);
            break;
            case 1:                     //REVERSE MATCH
                MT  =   "R";
                ReverseRead(Read,L);
                printf("R: %s\n",Read);
            break;
            case 2:                     //COMPLEMENT MATCH
                MT  =   "C";
                ComplementRead(Read,L);
                printf("C: %s\n",Read);
            break;
            case 3:                     //REVERSE COMPLEMENT MATCH
                MT  =   "E";
                ComplementRead(Read,L);
                ReverseRead(Read,L);
                printf("E: %s\n",Read);                
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

//BusqBin_Rul: función que utiliza el procedimiento de la Busq. Binaria
//             para ubicar el elemento seleccionado a través del mecanismo
//	            de la ruleta a fin de seleccionar los Padres.
int BusqBin_Rul(double prob[], int n, double x){
	int primero,ultimo,central;
	short encontrado;

	primero= 1 ; ultimo= n ; encontrado=0;
	while ((primero <= ultimo) && (encontrado==0)){
   	central = (primero + ultimo)/2;
      if (x == prob[central]) encontrado = 1;
      else
      	if (x > prob[central]) primero = central + 1;
         else ultimo = central - 1;
   }
	return(central);
}

