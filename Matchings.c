/*
 *  @Developers:    Juan Camilo Peña Vahos	- Aníbal Guerra	- Sebastian Isaza Ramirez
 *  @Last Revised:  28/07/2018
 *  @Description:   Módulo para las funciones que aplican los matching
*/ 

#include <stdio.h>
#include "Matchings.h"

//ReverseRead:	Implementa el inversor de reads
//@param:	    char *Read	:	Apuntador al arreglo que contiene el Read
//@param:	    long length	:	Longitud del Read
void ReverseRead(char *Read, long length){

	char aux;

	for (int i=0;   i<length/2;   i++){

        aux                 =   Read[length-i-1];
		Read[length-i-1]    =   Read[i];
		Read[i]             =   aux;

	}

}

//ComplementRead:	Implementa el complementador de reads (T,A)(C,G)
//@param:	        char *Read	:	Apuntador al arreglo que contiene el Read
//@param:	        long length	:	Longitud del Read
void ComplementRead(char *Read, long length){

	char Compl;
	int i;
	for (i=0;   i<length;   i++){
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