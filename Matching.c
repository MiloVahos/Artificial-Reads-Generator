/*
 ============================================================================
 Name        : 	Matching.c
 Author      : 	Juan Camilo Peña Vahos - Aníbal Guerra - Sebastian Isaza Ramírez
 Version     :
 Copyright   : 	This project is totally opensource
 Description :	Cabecera del archivo que contiene las funciones para aplicar los matchings
 ============================================================================
 */

#include <stdio.h>
#include <string.h>
#include "Matching.h"
#include <inttypes.h>

//ReverseRead:	Implementa el inversor de reads
//@param:	    char *Read	:	Apuntador al arreglo que contiene el Read
//@param:	    long length	:	Longitud del Read
void ReverseRead(char *Read, uint16_t L){

	char aux;

	for (int i=0;   i<L/2;   i++){
        	aux	=   Read[L-i-1];
		Read[L-i-1]    =   Read[i];
		Read[i]	=   aux;
	}
}

//ComplementRead:	Implementa el complementador de reads (T,A)(C,G)
//@param:	        char *Read	:	Apuntador al arreglo que contiene el Read
//@param:	        long length	:	Longitud del Read
void ComplementRead(char *Read, uint16_t L){

	char Compl;
	int i;
	for (i=0;   i<L;   i++){
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

//selMatching:		De acuerdo con MatTypeSel, aplica el matching seleccionado
//@param:			int MatTypeSel	:	Valor que determina el tipo de matching a aplicar
//@param:			int L			:	Longitud del read
//@param:			char *READ		:   Read
//@param:			char *MT		:	Matching Type
void selMatching(int MatTypeSel, uint16_t L, char *READ, char *MT){
	switch(MatTypeSel){
		case 0:                     //FORWARD MATCH
			memcpy(MT,"F",1);
		break;
		case 1:                     //REVERSE MATCH
			memcpy(MT,"R",1);
			ReverseRead(READ,L);
		break;
		case 2:                     //COMPLEMENT MATCH
			memcpy(MT,"C",1);
			ComplementRead(READ,L);
		break;
		case 3:                     //REVERSE COMPLEMENT MATCH
			memcpy(MT,"E",1);
			ComplementRead(READ,L);
			ReverseRead(READ,L);
		break;
		default: printf ("**Error in the matching selection, wrong input base**");
	}
}
