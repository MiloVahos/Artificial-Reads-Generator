/*
 ============================================================================
 Name        : 	ArtificialReadsGenerator.c
 Author      : 	Juan Camilo Peña Vahos - Aníbal Guerra - Sebastian Isaza Ramírez
 Version     :
 Copyright   : 	This project is totally opensource
 Description :	Archivo de cabecera para la funciones que manejan archivos
 ============================================================================
 */

#ifndef FILESUTILS_H
#define FILESUTILS_H
#include <inttypes.h>
#include <stdio.h>

#define LINE_LENGTH 1024
#define CHRM_LENGTH 1024

uint64_t 	contChars	    (char*);
void 		getReference	(char*,char*);
void		generateRead	(char*,uint32_t,uint16_t,char*,char*,FILE*,FILE*);
void		printCounters	(FILE*,uint16_t*);
int		    validBase	    (char);

#endif
