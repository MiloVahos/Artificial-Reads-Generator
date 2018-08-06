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
                case 1: if(c    !=  '\n')   cont++;    		break;
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
                case 1: if(c != '\n'){
					switch(c){
						case 'a':	case 'A':	Reference[j] = 'A'; break;
						case 'c':	case 'C':	Reference[j] = 'C'; break;
						case 'g':	case 'G':	Reference[j] = 'G'; break;
						case 't':	case 'T':	Reference[j] = 'T'; break;
						default:				Reference[j] = 'N';
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
	fprintf(FASTQ,"%s",Q);	//QUALITY SCORE

	//SE IMPRIME LA SECUENCIA SOLA EN EL ARCHIVO FASTQSEQ
	fprintf(FASTQSEQ,"Read: %"PRIu32", Secuencia: %s\n",id,FinalRead);

	free(FinalRead);

}

void printCounters(FILE* ALIGN,uint16_t* Cnt){
	//0->s 1->d 2->i 3->D 4->I 5->T 6->S 7->C
	fprintf(ALIGN,"Cantidad de mutaciones tipo 's': %d\n",Cnt[0]);
	fprintf(ALIGN,"Cantidad de mutaciones tipo 'd': %d\n",Cnt[1]);
	fprintf(ALIGN,"Cantidad de mutaciones tipo 'i': %d\n",Cnt[2]);
	fprintf(ALIGN,"Cantidad de mutaciones tipo 'D': %d\n",Cnt[3]);
	fprintf(ALIGN,"Cantidad de mutaciones tipo 'I': %d\n",Cnt[4]);
	fprintf(ALIGN,"Cantidad de mutaciones tipo 'T': %d\n",Cnt[5]);
	fprintf(ALIGN,"Cantidad de mutaciones tipo 'S': %d\n",Cnt[6]);
	fprintf(ALIGN,"Cantidad de mutaciones tipo 'C': %d\n",Cnt[7]);
}
/*char *load_reference(char *fname, int len){

	FILE *fp;
	fp	=	fopen(fname,"r");
		if(fp==NULL){
		fprintf(stderr, "Cannot open reference file %s \n",fname);fflush(stdout);
		exit(1);
	}

	char *buf;
	buf	=	(char*)	malloc	(sizeof(char)*LINE_LENGTH);
	if (buf==NULL) printf("Not enough memory for buf");

	char *ret;
	ret	=	(char*)	malloc	(sizeof(char)*CHRM_LENGTH);
	if (ret==NULL) printf("Not enough memory for buf");

	int i=0,cnt=0;

	while(fgets(buf, 1024, fp)!=NULL){

		if(buf[0]=='>') continue; //ignore comment lines

		for(i=0;	i<strlen(buf);	i++){ //ignore newline

			if(buf[i]>='A' && buf[i]<='z'){
				if((buf[i]=='a') || (buf[i]=='A')||(buf[i]=='c') || (buf[i]=='C')||
				   (buf[i]=='g') || (buf[i]=='G')||(buf[i]=='T') || (buf[i]=='t')){
					ret[cnt++]=buf[i];
				}else{
					ret[cnt++]='N';
				}
			}
		}
	}

	free(buf);
	fclose(fp);

	return ret;
}*/
