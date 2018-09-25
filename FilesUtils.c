/*
 ============================================================================
 Name        : 	FileUtils.c
 Author      : 	Juan Camilo Peña Vahos
 Version     :
 Copyright   : 	This project is totally opensource
 Description :	Funciones para la manipulación de la referencia
 ============================================================================
 */

#include "FilesUtils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>

//contChars: Cuenta el número de caractéres que tiene el archivo de referencia
//@param: char *FileName : Nombre del archivo
uint64_t contChars(char *FileName){

	uint64_t cont	=	0;	// Contador de las bases
    int estado;	// Determina si se esta en un comentario o no
	char c;	// Caracter que se esta leyendo en la referencia
	FILE *pf;	// Puntero al archivo

	pf = fopen(FileName,"r");
	if(pf!=NULL){
		while( (c = getc(pf)) != EOF ) {
		    if( c == '>' )    estado  =   0;	// EL ESTADO 0 SON LOS COMENTARIOS
			switch( estado ) {
				case 0: if( c == '\n' )		estado  =   1;  break;
				case 1: if( validBase(c) )	cont++;    		break;
			}
		}
	}
	fclose(pf);
	return cont;
}

//getReference: Obtiene la referencia de los archivos FASTA y la guarda en un arreglo
//@param: char *FileName : Nombre del archivo
//@param: char *Reference: Arreglo en el que se guarda la referencia
void getReference(char *FileName, char *Reference) {

	int estado;			// Estado en el que se está del documento, 0 es en comentarios
	char c;				// Caracter que se esta evaluando
	uint64_t j	=	0;	// Contador de la posición
	FILE *pf;			// Puntero al archivo de la referencia
	pf	=	fopen(FileName,"r");
	if(pf	!=	NULL){
		while((c = getc(pf)) != EOF){
		    if(c == '>')    estado  =   0;
		    switch(estado){
		        case 0: if(c == '\n')	estado  =   1;	break;
		        case 1: 
					if(validBase(c)){
						switch(c){
							case 'a':	case 'A':	Reference[j] = 'A'; break;
							case 'c':	case 'C':	Reference[j] = 'C'; break;
							case 'g':	case 'G':	Reference[j] = 'G'; break;
							case 't':	case 'T':	Reference[j] = 'T'; break;
						}
						j++;
		        	}
		        break;
		    }
		}
	}
	fclose(pf);
}

//generateRead:	Crea el read a partir de los datos de entrada y los imprime en el archivo FASTQ
//@param: char *Read:		El read que se va a imprimir en el archivo
//@param: uint32_t *id:		El ID del read que se va a imprimir
//@param: uint16_t L:		La longitud que se desea del read
//@param: char *Q:			Quality score
//@param: char *I:			Formato del ID de los reads
//@param: FILE *FASTQ:		Puntero al archivo de extensión FASTQ
//@param: FILE *FASTQSEQ:	Puntero al archivo de extensión FASTQSEQ
void generateRead(char *Read, uint32_t id, uint16_t L, char *Q, char *I, FILE *FASTQ, FILE *FASTQSEQ){

	//SE IMPRIME EL READ
	fprintf(FASTQ,"@%s:%"PRIu32"\n",I,id);	//ID
	fwrite(Read,sizeof(char),L,FASTQ);
	fprintf(FASTQ,"\n");
	//fprintf(FASTQ,"%s\n",FinalRead);	//SECUENCE
	fprintf(FASTQ,"+\n");			//EMPTY COMMENT
	fprintf(FASTQ,"%s\n",Q);			//QUALITY SCORE

	//SE IMPRIME LA SECUENCIA SOLA EN EL ARCHIVO FASTQSEQ
	fwrite(Read,sizeof(char),L,FASTQSEQ);
	fprintf(FASTQSEQ,"\n\n");

	//free(FinalRead);

}

//printCounter: Imprime los contadores de mutaciones de cada Read
//@param: FILE *ALIGN:	Puntero al archivo de extensión ALIGN
//@param: uint16_t Cnt:	Vector con los contadores que se van a imprimir.
void printCounters(FILE* ALIGN,uint16_t* Cnt){
	//0->s 1->d 2->i 3->D 4->I 5->T 6->S 7->C
	fprintf(ALIGN,"Mutations per Read counter\n");
	fprintf(ALIGN,"s: %"PRIu16" - d: %"PRIu16" - i: %"PRIu16" - D: %"PRIu16" - I: %"PRIu16" - T: %"PRIu16" - S: %"PRIu16" - C: %"PRIu16" \n",Cnt[0],Cnt[1],Cnt[2],Cnt[3],Cnt[4],Cnt[5],Cnt[6],Cnt[7]);
}

//validBase:	Determina si un caracter pertenece o no al alfabeto deseado
//@return:		1 (true) si pertence, 0 (false) si no pertenece
//@param:		char c:	Caracter que se va a evaluar	
int validBase(char c){
	return ( (c == 'A') || (c == 'C') || (c == 'T') || (c == 'G') ||
		 	 (c == 'a') || (c == 'c') || (c == 't') || (c == 'g') );
}



