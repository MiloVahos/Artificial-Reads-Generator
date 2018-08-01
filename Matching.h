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

void    ReverseRead     (char*,long);
void 	ComplementRead  (char*,long);
void	selMatching		(int,int,char*,char*);

#endif
