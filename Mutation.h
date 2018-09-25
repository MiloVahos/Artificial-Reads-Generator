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

#define MUTATION_TYPES	8	// NÚMERO DE MUTACIONES POSIBLES
#define READ_BIAS	1024	// SE OBTIENEN ALGUNAS BASES EXTRAS PARA OPERAR
#include <inttypes.h>
#include <stdio.h>

uint8_t	selMutation	       	(int);
char 	selBase		      	(int,char);
void	mutsVector		    (uint16_t,uint8_t*,uint16_t*,uint32_t*,double*);
void	offsetsGen		    (uint16_t,uint16_t*,uint16_t);
void    genRelOffsets       (uint16_t,uint16_t*,uint16_t*);
void    intercambiar        (uint16_t*, int, int);
void    ordenarOffsets      (uint16_t*, uint16_t);
// Oper - Read - Offsets - OffRel - BaseRef - BaseRead - BasesAcum - L - ALIGN
void    FordwardMutation	(   uint8_t*,char*,uint16_t*,uint16_t*,uint16_t,
                                uint8_t*,uint8_t*,double*,int,FILE*);
void    ReverseMutation     (   uint8_t*,char*,uint16_t*,uint16_t*,uint16_t,
                                uint8_t*,uint8_t*,double*,int,FILE*);


#endif
