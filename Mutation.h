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
#include <stdio.h>

uint8_t	selMutation	        (int);
char 	selBase		        (int,char);
void    FordwardMutation	(uint8_t,char*,uint16_t,uint8_t*,uint8_t*,int,int,FILE*);
void    ReverseMutation     (uint8_t,char*,uint16_t,uint8_t*,uint8_t*,int,int,FILE*);

#endif
