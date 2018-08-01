/*
 ============================================================================
 Name        : 	ArtificialReadsGenerator.c
 Author      : 	Juan Camilo Peña Vahos - Aníbal Guerra - Sebastian Isaza Ramírez
 Version     :
 Copyright   : 	This project is totally opensource
 Description :	Archivo de cabecera para las funciones de mutaciones
 ============================================================================
 */


#ifndef MUTATION_H
#define MUTATION_H

#define MUTATION_TYPES	8		//NÚMERO DE MUTACIONES POSIBLES
#include <inttypes.h>

uint8_t	selMutation	(int);
char 	selBase		(int,char);

#endif
