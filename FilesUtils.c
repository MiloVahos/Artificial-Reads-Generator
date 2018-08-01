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
long int contChars(char *FileName){
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
                case 1: if(c    !=  '\n')   cont++;    break;
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
	int j	=	0;
	FILE *pf;
	pf	=	fopen(FileName,"r");
	if(pf	!=	NULL){
		while((c = getc(pf)) != EOF){
            if(c == '>')    estado  =   0;
            switch(estado){
                case 0: if(c == '\n')	estado  =   1;	break;
                case 1: if(c    !=  '\n'){
					switch(c){
						case 'a':	case 'A':	Reference[j] = 'A'; break;
						case 'c':	case 'C':	Reference[j] = 'C'; break;
						case 'g':	case 'G':	Reference[j] = 'G'; break;
						case 't':	case 'T':	Reference[j] = 'T'; break;
						default:				Reference[j]	= 'N';
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
void generateRead(char *Read, uint32_t id, int L, char *Q, char *I, FILE *FASTQ, FILE *FASTQSEQ){

	//SE CORTA LA SECUENCIA DE MODO QUE QUEDE DE LONGITUD L
	char	*FinalRead	=	(char*) malloc(L*sizeof(char));
	memcpy(FinalRead,Read,L);

	//SE IMPRIME EL READ
	fprintf(FASTQ,"@%s %"PRIu32"\n",I,id);			//ID
	fprintf(FASTQ,"%s\n",FinalRead);				//SECUENCE
	fprintf(FASTQ,"+\n");							//EMPTY COMMENT
	for(int i =	0; i<L; i++) fprintf(FASTQ,"%s",Q);	//QUALITY SCORE

	//SE IMPRIME LA SECUENCIA SOLA EN EL ARCHIVO FASTQSEQ
	fprintf(FASTQSEQ,"Read: %"PRIu32", Secuencia: %s\n",id,FinalRead);

	free(FinalRead);

}

