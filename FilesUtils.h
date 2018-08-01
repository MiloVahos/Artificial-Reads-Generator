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

long int 	contChars	(char*);
void 		getReference(char*,char*);
void		generateRead(char*,uint32_t,int,char*,char*,FILE*,FILE*);

#endif
