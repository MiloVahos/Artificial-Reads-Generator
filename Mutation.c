/*
 ============================================================================
 Name        : 	Mutation.c
 Author      : 	Juan Camilo Peña Vahos - Aníbal Guerra - Sebastian Isaza Ramírez
 Version     :
 Copyright   : 	This project is totally opensource
 Description :	Funciones para aplicar mutaciones a los reads
 ============================================================================
 */
#include "Mutation.h"
#include <inttypes.h>
#include <stdio.h>
#include "Stats.h"
#include <stdlib.h>

//selMutation:	Determina a que tipo de mutación corresponde el lanzamiento del dado
//@param:	int ErrorSel	:	Valor resultante de la búsqueda binaria
uint8_t selMutation(int ErrorSel){
	switch(ErrorSel){
		case 0: return 's';	break;
		case 1: return 'd'; break;
		case 2: return 'i'; break;
		case 3: return 'D'; break;
		case 4: return 'I'; break;
		case 5: return 'T'; break;
		case 6: return 'S'; break;
		case 7: return 'C'; break;
		default: return 's';
	}
}

//selBase:	Selecciona una base de acuerdo a sel y a la Base de entrada
//@param:	int sel		:	Valor de la selección
//@param:	int Base	:	Base de entrada
char selBase(int sel,char Base){
	switch(Base){
		case 'A':
			switch(sel){
				case 0: return 'C'; break;
				case 1: return 'G'; break;
				case 2: return 'T'; break;
				case 3: return 'N'; break;
			}
		break;
		case 'C':
			switch(sel){
				case 0: return 'A'; break;
				case 1: return 'G'; break;
				case 2: return 'T'; break;
				case 3: return 'N'; break;
			}
		break;
		case 'G':
			switch(sel){
				case 0: return 'A'; break;
				case 1: return 'C'; break;
				case 2: return 'T'; break;
				case 3: return 'N'; break;
			}
		break;
		case 'T':
			switch(sel){
				case 0: return 'A'; break;
				case 1: return 'C'; break;
				case 2: return 'G'; break;
				case 3: return 'N'; break;
			}
		break;
	}
}

//mutsVector:	Genera un vector con todas las operaciones de mutación que se van a aplicar
void mutsVector (uint16_t lendesc,uint8_t *Oper,uint16_t *Cnt,double *MutTypeAcumF){
	for(int i=0;	i<lendesc;	i++){
		double dado		=	LanzarDado();
		int ErrorSel	=	BusqBin_Rul(MutTypeAcumF,MUTATION_TYPES,dado);
		Oper[i]			=	selMutation(ErrorSel);
		switch((char)Oper[i]){
		//0->s 1->d 2->i 3->D 4->I 5->T 6->S 7->C
			case 's': Cnt[0]++; break;
			case 'i': Cnt[2]++; break;
			case 'I': Cnt[4]++; break;
			case 'S': Cnt[6]++; break;
		}
	}
}

//offsetsGen:	Genera un vector con los offsets donde se van a aplicar las mutaciones
void offsetsGen	(uint16_t lendesc,uint16_t *Offsets, uint16_t L){

	//1. GENERAR TODOS LOS OFFSETS EN DESORDEN
	for ( int i = 0; i < lendesc; i++ ) {

		int flag = 0;
		int OffsetAux = 0;
		do {
			OffsetAux	=	rand() %(((L-1)-(lendesc - i)) + 1 - 0) + 0;
			if ( i != 0 ) {
				int contador = 0;
				for (int j = 0; j < i; j++ ) {
					if ( OffsetAux	== Offsets[j] ){
						contador++;
					} 
				}
				if( contador == 0){
					flag = 1;
				}
			} else {
				flag = 1;
			}
			
		} while ( flag != 1 );
		Offsets[i]	=	OffsetAux;
	}
	// 2. ORDENAR EL VECTOR
    ordenarOffsets(Offsets,lendesc);
}

void genRelOffsets (uint16_t lendesc,uint16_t* Offsets,uint16_t* OffRel) { 
	OffRel[0]	=	Offsets[0];
	for ( int i = 1; i < lendesc; i++ ) {
		OffRel[i]	=	Offsets[i]	-	Offsets[i-1];
	}
}


void intercambiar (uint16_t *Offsets, int i, int j) {
    int tmp = Offsets[i];
    Offsets[i] = Offsets[j];
    Offsets[j] = tmp;
}

void ordenarOffsets (uint16_t *Offsets, uint16_t N){

    int i, j, k;
    for (i = 0; i < N - 1; i++) {
        for (k = i, j = i + 1; j < N; j++)
            if (Offsets[j] < Offsets[k]) k = j;
            if (k != i) intercambiar (Offsets, i, k);
    }
}

void FordwardMutation (	uint8_t *Oper,char *Read,uint16_t *Offsets,uint16_t *OffRel,uint16_t lendesc,
						uint8_t *BaseRef,uint8_t *BaseRead,double* BasesAcum,int L,FILE *ALIGN) {

	int posRead = 0;
	for ( int i = 0; i < lendesc; i++ ) {
		fprintf(ALIGN,"Operación de mutación: %c\n",Oper[i]);
		fprintf(ALIGN,"Offset[%d] =	%d, \n",i,Offsets[i]);

		switch ( (char) Oper[i] ){

			case 's':

				//CALCULAR LA BASE DESTINO
				double dado		=	LanzarDado();
				int sel			=	BusqBin_Rul(BasesAcum,4,dado);
				char BaseDest	=	selBase(sel,Read[Offsets[i]]);

				
				posRead			= 	posRead + OffRel[i];
				Read[posRead]	=	BaseDest;
				BaseRef[i]		=	Read[posRead];
				BaseRead[i]		=	BaseDest;
				posRead++;

				//IMPRIMIR RESULTADOS
				fprintf(ALIGN,"Base reference %c\n",BaseRef[i]);
				fprintf(ALIGN,"Base read: %c\n",BaseRead[i]);
				fprintf(ALIGN,"SECUENCE: %s\n",Read);

			break;
			case 'd':

				//APLICAR LA MUTACIÓN Y GUARDAR LAS BASES
				posRead	= 	posRead + OffRel[i];
				for( int q = posRead; q < L-1; q++ ){
  					Read[q]	=	Read[q+1];
  				}
					
				BaseRef[i] 	=	Read[posRead];
				BaseRead[i]	=	'0';

				//IMPRIMIR RESULTADOS
				fprintf(ALIGN,"Base reference %c\n",BaseRef[i]);
				fprintf(ALIGN,"Base read: %c\n",BaseRead[i]);
				fprintf(ALIGN,"SECUENCE: %s\n",Read);

			break;
			case 'i':

				//APLICAR LA MUTACIÓN Y GUARDAR LAS BASES
				posRead	= 	posRead + OffRel[i];
				for( int q = (L-1); q > posRead; q--){
  					Read[q]	=	Read[q-1];
  				}
				BaseRef[i] 	=	'0';
				BaseRead[i]	=	Read[posRead];
				// Read[posRead]	=	aux; 	//PREGUNTAR ESTO
				posRead++;
				//IMPRIMIR RESULTADOS
				fprintf(ALIGN,"Base reference %c\n",BaseRef[i]);
				fprintf(ALIGN,"Base read: %c\n",BaseRead[i]);
				fprintf(ALIGN,"SECUENCE: %s\n",Read);

			break;
			case 'D':

				//APLICAR LA MUTACIÓN Y GUARDAR LAS BASES
				posRead	=	posRead	+ OffRel[i];
				for( int q = posRead; q< L-2; q++ ){
  			   		Read[q]	=	Read[q+2];
  			   	}
				BaseRef[i] 	=	Read[posRead];
				BaseRead[i]	=	'0';
				//IMPRIMIR RESULTADOS
				fprintf(ALIGN,"Base reference %c\n",BaseRef[i]);
				fprintf(ALIGN,"Base read: %c\n",BaseRead[i]);
				fprintf(ALIGN,"SECUENCE: %s\n",Read);

			break;
			case 'I':

				//APLICAR LA MUTACIÓN Y GUARDAR LAS BASES
				posRead	=	posRead + OffRel[i];
				for( int q = (L-1); q > posRead; q-- ) {
					Read[q]	=	Read[q-1];
				}
				BaseRef[i] 		=	Read[posRead];
				BaseRead[i]		=	'N';
				Read[posRead]	=	'N';
				posRead++;

				//IMPRIMIR RESULTADOS
				fprintf(ALIGN,"Base reference %c\n",BaseRef[i]);
				fprintf(ALIGN,"Base read: %c\n",BaseRead[i]);
				fprintf(ALIGN,"SECUENCE: %s\n",Read);

			break;
			case 'T':
				
				//APLICAR LA MUTACIÓN Y GUARDAR LAS BASES
				posRead	=	posRead	+ OffRel[i];
				for( int q = posRead; q< L-3; q++ ){
					Read[q]	=	Read[q+3];
				}
				BaseRef[i] 	=	Read[posRead];
				BaseRead[i]	=	'0';
				//IMPRIMIR RESULTADOS
				fprintf(ALIGN,"Base reference %c\n",BaseRef[i]);
				fprintf(ALIGN,"Base read: %c\n",BaseRead[i]);
				fprintf(ALIGN,"SECUENCE: %s\n",Read);
			break;
			case 'S':

				//CALCULAR LA BASE DESTINO
				int flag = 0;
				posRead	=	posRead	+ OffRel[i];
				do{
					double dado		=	LanzarDado();
					int sel			=	BusqBin_Rul(BasesAcum,4,dado);
					char BaseDest	=	selBase(sel,Read[posRead]);
					if(BaseDest != Read[posRead+1]){
						//GUARDAR LAS BASES Y APLICAR LA MUTACIÓN
						Read[posRead]	=	BaseDest;
						BaseRef[i]		=	Read[posRead];
						BaseRead[i]		=	BaseDest;
						posRead++;
						Read[posRead]	=	BaseDest;
						BaseRef[i]		=	Read[posRead];
						BaseRead[i]		=	BaseDest;
						posRead++;
						flag = 1;
					}

				}while(flag != 1);


				fprintf(ALIGN,"Bases reference %c\n",BaseRef[i]);
				fprintf(ALIGN,"Base read: %c\n",BaseRead[i]);
				fprintf(ALIGN,"SECUENCE: %s\n",Read);

			break;
			case 'C':
				//APLICAR LA MUTACIÓN Y GUARDAR LAS BASES
				posRead	=	posRead	+ OffRel[i];
				for( int q = posRead; q< L-4; q++ ){
  			   		Read[q]	=	Read[q+4];
  			   	}
				BaseRef[i] 	=	Read[posRead];
				BaseRead[i]	=	'0';
				//IMPRIMIR RESULTADOS
				fprintf(ALIGN,"Base reference %c\n",BaseRef[i]);
				fprintf(ALIGN,"Base read: %c\n",BaseRead[i]);
				fprintf(ALIGN,"SECUENCE: %s\n",Read);
			break;
		}
	}
 }


void ReverseMutation (	uint8_t *Oper,char *Read,uint16_t *Offsets,uint16_t *OffRel,uint16_t lendesc,
						uint8_t *BaseRef,uint8_t *BaseRead,double* BasesAcum,int L,FILE *ALIGN) {
	
	int posRead = L-1;
	for ( int i = 0; i < lendesc; i++ ) {

		fprintf(ALIGN,"Operación de mutación: %c\n",Oper[i]);
		fprintf(ALIGN,"Offset[%d] =	%d, \n",i,Offsets[i]);

		switch ( (char) Oper[i] ){

			case 's':

				//CALCULAR LA BASE DESTINO
				double dado		=	LanzarDado();
				int sel			=	BusqBin_Rul(BasesAcum,4,dado);
				char BaseDest	=	selBase(sel,Read[Offsets[i]]);

				
				posRead			= 	posRead - OffRel[i];
				Read[posRead]	=	BaseDest;
				BaseRef[i]		=	Read[posRead];
				BaseRead[i]		=	BaseDest;
				posRead--;

				//IMPRIMIR RESULTADOS
				fprintf(ALIGN,"Base reference %c\n",BaseRef[i]);
				fprintf(ALIGN,"Base read: %c\n",BaseRead[i]);
				fprintf(ALIGN,"SECUENCE: %s\n",Read);

			break;
			case 'd':

				//APLICAR LA MUTACIÓN Y GUARDAR LAS BASES
				posRead	= 	posRead - OffRel[i];
				for( int q = posRead; q < L-1; q++ ){
  					Read[q]	=	Read[q+1];
  				}
					
				BaseRef[i] 	=	Read[posRead];
				BaseRead[i]	=	'0';
				posRead--;

				//IMPRIMIR RESULTADOS
				fprintf(ALIGN,"Base reference %c\n",BaseRef[i]);
				fprintf(ALIGN,"Base read: %c\n",BaseRead[i]);
				fprintf(ALIGN,"SECUENCE: %s\n",Read);

			break;
			case 'i':

				//APLICAR LA MUTACIÓN Y GUARDAR LAS BASES
				posRead	= 	posRead - OffRel[i];
				for( int q = (L-1); q > posRead+1; q--){
  					Read[q]	=	Read[q-1];
  				}
				BaseRef[i] 	=	'0';
				BaseRead[i]	=	Read[posRead];
				// Read[posRead]	=	aux; 	//PREGUNTAR ESTO

				//IMPRIMIR RESULTADOS
				fprintf(ALIGN,"Base reference %c\n",BaseRef[i]);
				fprintf(ALIGN,"Base read: %c\n",BaseRead[i]);
				fprintf(ALIGN,"SECUENCE: %s\n",Read);

			break;
			case 'D':

				//APLICAR LA MUTACIÓN Y GUARDAR LAS BASES
				posRead	=	posRead	- OffRel[i];
				for( int q = posRead-1; q< L-2; q++ ){
  			   		Read[q]	=	Read[q+2];
  			   	}
				BaseRef[i] 	=	Read[posRead];
				BaseRead[i]	=	'0';
				posRead 	= 	posRead-2;
				//IMPRIMIR RESULTADOS
				fprintf(ALIGN,"Base reference %c\n",BaseRef[i]);
				fprintf(ALIGN,"Base read: %c\n",BaseRead[i]);
				fprintf(ALIGN,"SECUENCE: %s\n",Read);

			break;
			case 'I':

				//APLICAR LA MUTACIÓN Y GUARDAR LAS BASES
				posRead	=	posRead - OffRel[i];
				for( int q = (L-1); q > posRead+1; q-- ) {
					Read[q]	=	Read[q-1];
				}
				BaseRef[i] 		=	Read[posRead];
				BaseRead[i]		=	'N';
				Read[posRead+1]	=	'N';

				//IMPRIMIR RESULTADOS
				fprintf(ALIGN,"Base reference %c\n",BaseRef[i]);
				fprintf(ALIGN,"Base read: %c\n",BaseRead[i]);
				fprintf(ALIGN,"SECUENCE: %s\n",Read);

			break;
			case 'T':
				
				//APLICAR LA MUTACIÓN Y GUARDAR LAS BASES
				posRead	=	posRead	- OffRel[i];
				for( int q = posRead-2; q< L-3; q++ ){
					Read[q]	=	Read[q+3];
				}
				BaseRef[i] 	=	Read[posRead];
				BaseRead[i]	=	'0';
				posRead		=	posRead-3;
				//IMPRIMIR RESULTADOS
				fprintf(ALIGN,"Base reference %c\n",BaseRef[i]);
				fprintf(ALIGN,"Base read: %c\n",BaseRead[i]);
				fprintf(ALIGN,"SECUENCE: %s\n",Read);
			break;
			case 'S':

				//CALCULAR LA BASE DESTINO
				int flag = 0;
				posRead	=	posRead	- OffRel[i];
				do{
					double dado		=	LanzarDado();
					int sel			=	BusqBin_Rul(BasesAcum,4,dado);
					char BaseDest	=	selBase(sel,Read[posRead]);
					if(BaseDest != Read[posRead+1]){
						//GUARDAR LAS BASES Y APLICAR LA MUTACIÓN
						Read[posRead]	=	BaseDest;
						BaseRef[i]		=	Read[posRead];
						BaseRead[i]		=	BaseDest;
						posRead--;
						Read[posRead]	=	BaseDest;
						BaseRef[i]		=	Read[posRead];
						BaseRead[i]		=	BaseDest;
						posRead--;
						flag = 1;
					}

				}while(flag != 1);


				fprintf(ALIGN,"Bases reference %c\n",BaseRef[i]);
				fprintf(ALIGN,"Base read: %c\n",BaseRead[i]);
				fprintf(ALIGN,"SECUENCE: %s\n",Read);

			break;
			case 'C':
				//APLICAR LA MUTACIÓN Y GUARDAR LAS BASES
				posRead	=	posRead	- OffRel[i];
				for( int q = posRead-3; q< L-4; q++ ){
  			   		Read[q]	=	Read[q+4];
  			   	}
				BaseRef[i] 	=	Read[posRead];
				BaseRead[i]	=	'0';
				posRead		=	posRead - 4;
				//IMPRIMIR RESULTADOS
				fprintf(ALIGN,"Base reference %c\n",BaseRef[i]);
				fprintf(ALIGN,"Base read: %c\n",BaseRead[i]);
				fprintf(ALIGN,"SECUENCE: %s\n",Read);
			break;
		}
	}

}
