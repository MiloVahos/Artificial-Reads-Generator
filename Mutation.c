/*
 ============================================================================
 Name        : 	Mutation.c
 Author      : 	Juan Camilo Peña Vahos - Aníbal Guerra - Sebastian Isaza Ramírez
 Version     :
 Copyright   : 	This project is totally opensource
 Description :	Funciones para aplicar mutaciones a los reads
 ============================================================================
 */
#include "Mutation.h"
#include <inttypes.h>

//selMutation:	Determina a que tipo de mutación corresponde el lanzamiento del dado
//@param:	int ErrorSel	:	Valor resultante de la búsqueda binaria
uint8_t selMutation(int ErrorSel){
	switch(ErrorSel){
		case 0: return 's';	break;
		case 1: return 'd'; break;
		case 2: return 'i'; break;
		case 3: return 'D'; break;
		case 4: return 'I'; break;
		case 5: return 'T'; break;
		case 6: return 'S'; break;
		case 7: return 'C'; break;
		default: return 's';
	}
}

//selBase:	Selecciona una base de acuerdo a sel y a la Base de entrada
//@param:	int sel		:	Valor de la selección
//@param:	int Base	:	Base de entrada
char selBase(int sel,char Base){
	switch(Base){
		case 'A':
			switch(sel){
				case 0: return 'C'; break;
				case 1: return 'G'; break;
				case 2: return 'T'; break;
				case 3: return 'N'; break;
			}
		break;
		case 'C':
			switch(sel){
				case 0: return 'A'; break;
				case 1: return 'G'; break;
				case 2: return 'T'; break;
				case 3: return 'N'; break;
			}
		break;
		case 'G':
			switch(sel){
				case 0: return 'A'; break;
				case 1: return 'C'; break;
				case 2: return 'T'; break;
				case 3: return 'N'; break;
			}
		break;
		case 'T':
			switch(sel){
				case 0: return 'A'; break;
				case 1: return 'C'; break;
				case 2: return 'G'; break;
				case 3: return 'N'; break;
			}
		break;
	}
}
