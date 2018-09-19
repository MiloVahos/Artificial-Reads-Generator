/*
 ============================================================================
 Name        : 	File.c
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
	long int cont	=	0;
    int estado;
	char c;
	FILE *pf;

	pf = fopen(FileName,"r");
	if(pf!=NULL){
		while((c = getc(pf)) != EOF){
		    if(c == '>')    estado  =   0;		//EL ESTADO 0 SON LOS COMENTARIOS
			    switch(estado){
				case 0: if(c    ==  '\n')   estado  =   1;  break;
				case 1: if(validBase(c))   cont++;    		break;
			    }
		}
	}

	return cont;
}

//getReference: Obtiene la referencia de los archivos FASTA y la guarda en un arreglo
//@param: char *FileName : Nombre del archivo
//@param: char *Reference: Arreglo en el que se guarda la referencia
void getReference(char *FileName, char *Reference){
	int estado;
	char c;
	uint64_t j	=	0;
	FILE *pf;
	pf	=	fopen(FileName,"r");
	if(pf	!=	NULL){
		while((c = getc(pf)) != EOF){
		    if(c == '>')    estado  =   0;
		    switch(estado){
		        case 0: if(c == '\n')	estado  =   1;	break;
		        case 1: if(validBase(c)){
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
void generateRead(char *Read, uint32_t id, uint16_t L, char *Q, char *I, FILE *FASTQ, FILE *FASTQSEQ){

	//SE CORTA LA SECUENCIA DE MODO QUE QUEDE DE LONGITUD L
	//char	*FinalRead	=	(char*) malloc(L*sizeof(char));
	//memcpy(FinalRead,Read,L);

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
void printCounters(FILE* ALIGN,uint16_t* Cnt){
	//0->s 1->d 2->i 3->D 4->I 5->T 6->S 7->C
	fprintf(ALIGN,"Mutations per Read counter\n");
	fprintf(ALIGN,"s: %"PRIu16" - d: %"PRIu16" - i: %"PRIu16" - D: %"PRIu16" - I: %"PRIu16" - T: %"PRIu16" - S: %"PRIu16" - C: %"PRIu16" \n",Cnt[0],Cnt[1],Cnt[2],Cnt[3],Cnt[4],Cnt[5],Cnt[6],Cnt[7]);
}


int validBase(char c){
	return ( (c == 'A') || (c == 'C') || (c == 'T') || (c == 'G') ||
		 (c == 'a') || (c == 'c') || (c == 't') || (c == 'g') );
}



