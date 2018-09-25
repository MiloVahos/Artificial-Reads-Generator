/*
 ============================================================================
 Name        : 	ArtificialReadsGenerator.c
 Author      : 	Juan Camilo Peña Vahos - Aníbal Guerra - Sebastian Isaza Ramírez
 Version     :
 Copyright   : 	This project is totally opensource
 Description :	Funciones de Matching
 ============================================================================
 */

#ifndef MATCHING_H
#define MATCHING_H

#define MATCHING_TYPES 	4		//NÚMERO DE MATCHINGS POSIBLES

#include <inttypes.h>

void    ReverseRead     (char*,uint16_t);
void 	ComplementRead  (char*,uint16_t);
void	selMatching	    (int,uint16_t,char*,char*);

#endif
