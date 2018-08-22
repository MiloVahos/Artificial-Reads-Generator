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


void mutsVector (uint16_t lendesc,uint8_t *Oper,uint16_t *Cnt,double *MutTypeAcumF){
	for(int i=0;	i<lendesc;	i++){
		double dado	=	LanzarDado();
		int ErrorSel	=	BusqBin_Rul(MutTypeAcumF,MUTATION_TYPES,dado);
		Oper[i]		=	selMutation(ErrorSel);
		switch((char)Oper[i]){
		//0->s 1->d 2->i 3->D 4->I 5->T 6->S 7->C
			case 's': Cnt[0]++; break;
			case 'i': Cnt[2]++; break;
			case 'I': Cnt[4]++; break;
			case 'S': Cnt[6]++; break;
		}
	}
}

void offsetsGen	(uint16_t lendesc,char strand,uint16_t *Offsets,uint8_t *Oper, uint16_t L){

	uint32_t OffsetAnt	=	0;	//ACUMULADOR DE OFFSETS
	uint32_t OffsetAux	=	0;	//OFFSET ACTUAL
	for (int i=0;	i<lendesc;	i++){
		if((strand=='f')||(strand=='c')){
			if(i==0){
				OffsetAux 	= 	rand() %(((L-1)-lendesc) + 1 - 0) + 0;
				Offsets[i]	=	OffsetAux;
			}else{
				switch(Oper[i-1]){
					case 's': case 'S': case 'i': case 'I':
						OffsetAux	=	rand() %(((L-1)-(lendesc-i)) + 1 - (OffsetAnt+1)) + (OffsetAnt+1);
						Offsets[i]	=	OffsetAux	-	OffsetAnt;
					break;
					default:
						OffsetAux	=	rand() %(((L-1)-(lendesc-i)) + 1 - (OffsetAnt)) + (OffsetAnt);
						Offsets[i]	=	OffsetAux	-	OffsetAnt;
				}
				if(OffsetAux>1024) printf("forward mayor 1024\n");
				if(OffsetAux<0) printf("forward menor 0\n");
			}
		}else{
			if(i==0){
				OffsetAux = rand() %((L-1) + 1 - (0+lendesc)) + (0+lendesc);
				Offsets[i]	=	OffsetAux;
			}else{
				switch(Oper[i-1]){
					case 's': case 'S': case 'd':
						OffsetAux = rand() %((OffsetAnt) + 1 - (0+(lendesc-i))) + (0+(lendesc-i));
						Offsets[i]	=	OffsetAnt - OffsetAux;
					break;
					case 'D':
						OffsetAux = rand() %((OffsetAnt) + 1 - (0+(lendesc-i))) + (0+(lendesc-i));
						Offsets[i]	=	OffsetAnt - OffsetAux;
					break;
					case 'T':
						OffsetAux = rand() %((OffsetAnt) + 1 - (0+(lendesc-i))) + (0+(lendesc-i));
						Offsets[i]	=	OffsetAnt - OffsetAux;
					break;
					case 'C':
						OffsetAux = rand() %((OffsetAnt) + 1 - (0+(lendesc-i))) + (0+(lendesc-i));
						Offsets[i]	=	OffsetAnt - OffsetAux;
					break;
					default:
						OffsetAux = rand() %((OffsetAnt) + 1 - (0+(lendesc-i))) + (0+(lendesc-i));
						Offsets[i]	=	OffsetAnt - OffsetAux;
				}
				printf("%"PRIu32"\n",OffsetAux);
				if(OffsetAux>1024) printf("Reverse mayor 1024\n");
				if(OffsetAux<0) printf("Reverse menor 0\n");
			}
		}
		printf("r = %d\n",((L-1)-(lendesc-i)) + 1 - (OffsetAnt));
		OffsetAnt	=	OffsetAux;
	}

}



void FordwardMutation(uint8_t Oper,char *Read,uint16_t Offset, uint8_t *BaseRef,uint8_t *BaseRead, int L, int BaseActual, FILE *ALIGN){
 	//VECTOR DE PROBABILIDADES DE LAS BASES PARA LAS SUBSTITUCIONES
	double	BasesStats[4]	=	{0.3,0.3,0.3,0.1};
	double 	BasesAcum[4];

 	switch((char)Oper){
		case 's':
			BasesAcum[0]		= 	BasesStats[0];
			for (int i=1; i<4; i++){
				BasesAcum[i]= BasesStats[i] + BasesAcum[i-1];
			}

			//CALCULAR LA BASE DESTINO
			double dado		=	LanzarDado();
			int sel	=	BusqBin_Rul(BasesAcum,4,dado);
			char BaseDest			=	selBase(sel,Read[Offset]);

			//APLICAR LA MUTACIÓN Y GUARDAR LAS BASES
			BaseRef[BaseActual] 	=	Read[Offset];
			BaseRead[BaseActual]	=	BaseDest;
			Read[Offset]	=	BaseDest;

			//IMPRIMIR RESULTADOS
			fprintf(ALIGN,"Base reference %c\n",BaseRef[BaseActual]);
			fprintf(ALIGN,"Base read: %c\n",BaseRead[BaseActual]);
			fprintf(ALIGN,"SECUENCE: %s\n",Read);

		break;
		case 'd':

			//APLICAR LA MUTACIÓN Y GUARDAR LAS BASES
			BaseRef[BaseActual] 	=	Read[Offset];
			BaseRead[BaseActual]	=	'0';
			for(int q=Offset; q<L-1;q++){
			  Read[q] = Read[q+1];
			}

			//IMPRIMIR RESULTADOS
			fprintf(ALIGN,"Base reference %c\n",BaseRef[BaseActual]);
			fprintf(ALIGN,"Base read: %c\n",BaseRead[BaseActual]);
			fprintf(ALIGN,"SECUENCE: %s\n",Read);

		break;
		case 'i':

			//APLICAR LA MUTACIÓN Y GUARDAR LAS BASES
			BaseRef[BaseActual] 	=	'0';
			BaseRead[BaseActual]	=	Read[Offset];
			for(int q=(L-1);q>Offset;q--){
				Read[q]	=	Read[q-1];
			}

			//IMPRIMIR RESULTADOS
			fprintf(ALIGN,"Base reference %c\n",BaseRef[BaseActual]);
			fprintf(ALIGN,"Base read: %c\n",BaseRead[BaseActual]);
			fprintf(ALIGN,"SECUENCE: %s\n",Read);

		break;
		case 'D':

			//APLICAR LA MUTACIÓN Y GUARDAR LAS BASES
			BaseRef[BaseActual] 	=	Read[Offset];
			BaseRead[BaseActual]	=	'0';
			for(int q=Offset;q<L-2;q++){
				Read[q] = Read[q+2];
			}

			//IMPRIMIR RESULTADOS
			fprintf(ALIGN,"Base reference %c\n",BaseRef[BaseActual]);
			fprintf(ALIGN,"Base read: %c\n",BaseRead[BaseActual]);
			fprintf(ALIGN,"SECUENCE: %s\n",Read);

		break;
		case 'I':

			//APLICAR LA MUTACIÓN Y GUARDAR LAS BASES
			BaseRef[BaseActual] 	=	'0';
			BaseRead[BaseActual]	=	Read[Offset];

			for(int q=(L-1);q>Offset;q--){
				Read[q]	=	Read[q-1];
			}
			Read[Offset]	=	'N';

			//IMPRIMIR RESULTADOS
			fprintf(ALIGN,"Base reference %c\n",BaseRef[BaseActual]);
			fprintf(ALIGN,"Base read: %c\n",BaseRead[BaseActual]);
			fprintf(ALIGN,"SECUENCE: %s\n",Read);

		break;
		case 'T':
			//APLICAR LA MUTACIÓN Y GUARDAR LAS BASES
			BaseRef[BaseActual] 	=	Read[Offset];
			BaseRead[BaseActual]	=	'0';
			for(int q=Offset;q<L-3;q++){
				Read[q] = Read[q+3];
			}

			//IMPRIMIR RESULTADOS
			fprintf(ALIGN,"Base reference %c\n",BaseRef[BaseActual]);
			fprintf(ALIGN,"Base read: %c\n",BaseRead[BaseActual]);
			fprintf(ALIGN,"SECUENCE: %s\n",Read);

		break;
		case 'S':

			BasesAcum[0]		= 	BasesStats[0];
			for (int i=1; i<4; i++){
				BasesAcum[i]= BasesStats[i] + BasesAcum[i-1];
			}

			//CALCULAR LA BASE DESTINO
			int flag = 0;
			do{
				double dado		=	LanzarDado();
				int sel			=	BusqBin_Rul(BasesAcum,4,dado);
				char BaseDest	=	selBase(sel,Read[Offset]);
				if(BaseDest != Read[Offset+1]){
					//GUARDAR LAS BASES Y APLICAR LA MUTACIÓN
					BaseRef[BaseActual]		=	Read[Offset];
					BaseRef[BaseActual+1]	=	Read[Offset+1];
					BaseRead[BaseActual]	=	BaseDest;
					BaseRead[BaseActual+1]	=	BaseDest;
					Read[Offset]	=	BaseDest;
					Read[Offset+1]	=	BaseDest;
					flag = 1;
				}

			}while(flag != 1);


			fprintf(ALIGN,"Bases reference %c\n",BaseRef[BaseActual]);
			fprintf(ALIGN,"Base read: %c\n",BaseRead[BaseActual]);
			fprintf(ALIGN,"SECUENCE: %s\n",Read);

		break;
		case 'C':
			//APLICAR LA MUTACIÓN Y GUARDAR LAS BASES
			BaseRef[BaseActual] 	=	Read[Offset];
			BaseRead[BaseActual]	=	'0';

			for(int q=Offset;q<L-4;q++){
				Read[q] = Read[q+4];
			}

			//IMPRIMIR RESULTADOS
			fprintf(ALIGN,"Base reference %c\n",BaseRef[BaseActual]);
			fprintf(ALIGN,"Base read: %c\n",BaseRead[BaseActual]);
			fprintf(ALIGN,"SECUENCE: %s\n",Read);
		break;
	}
 }

void ReverseMutation(uint8_t Oper,char *Read,uint16_t Offset, uint8_t *BaseRef,uint8_t *BaseRead, int L, int BaseActual, FILE *ALIGN){
 	//VECTOR DE PROBABILIDADES DE LAS BASES PARA LAS SUBSTITUCIONES
	double	BasesStats[4]	=	{0.3,0.3,0.3,0.1};
	double 	BasesAcum[4];
 	switch((char)Oper){

		case 's':
			BasesAcum[0]		= 	BasesStats[0];
			for (int i=1; i<4; i++){
				BasesAcum[i]= BasesStats[i] + BasesAcum[i-1];
			}

			//CALCULAR LA BASE DESTINO
			double dado		=	LanzarDado();
			int sel	=	BusqBin_Rul(BasesAcum,4,dado);
			char BaseDest			=	selBase(sel,Read[Offset]);

			//APLICAR LA MUTACIÓN Y GUARDAR LAS BASES
			BaseRef[BaseActual] 	=	Read[Offset];
			BaseRead[BaseActual]	=	BaseDest;
			Read[Offset]	=	BaseDest;

			//IMPRIMIR RESULTADOS
			fprintf(ALIGN,"Base reference %c\n",BaseRef[BaseActual]);
			fprintf(ALIGN,"Base read: %c\n",BaseRead[BaseActual]);
			fprintf(ALIGN,"SECUENCE: %s\n",Read);

		break;
		case 'd':

			//APLICAR LA MUTACIÓN Y GUARDAR LAS BASES
			BaseRef[BaseActual] 	=	Read[Offset];
			BaseRead[BaseActual]	=	'0';
			for(int q=Offset; q<L-1;q++){
			  Read[q] = Read[q+1];
			}

			//IMPRIMIR RESULTADOS

			fprintf(ALIGN,"Base reference %c\n",BaseRef[BaseActual]);
			fprintf(ALIGN,"Base read: %c\n",BaseRead[BaseActual]);
			fprintf(ALIGN,"SECUENCE: %s\n",Read);

		break;
		case 'i':

			//APLICAR LA MUTACIÓN Y GUARDAR LAS BASES
			BaseRef[BaseActual] 	=	'0';
			BaseRead[BaseActual]	=	Read[Offset];
			for(int q=(L-1);q>Offset+1;q--){
				Read[q]	=	Read[q-1];
			}

			//IMPRIMIR RESULTADOS
			fprintf(ALIGN,"Base reference %c\n",BaseRef[BaseActual]);
			fprintf(ALIGN,"Base read: %c\n",BaseRead[BaseActual]);
			fprintf(ALIGN,"SECUENCE: %s\n",Read);
		break;
		case 'D':

			//APLICAR LA MUTACIÓN Y GUARDAR LAS BASES
			BaseRef[BaseActual] 	=	Read[Offset];
			BaseRead[BaseActual]	=	'0';
			for(int q=Offset-1;q<L-2;q++){
				Read[q] = Read[q+2];
			}

			//IMPRIMIR RESULTADOS
			fprintf(ALIGN,"Base reference %c\n",BaseRef[BaseActual]);
			fprintf(ALIGN,"Base read: %c\n",BaseRead[BaseActual]);
			fprintf(ALIGN,"SECUENCE: %s\n",Read);
		break;
		case 'I':

			//APLICAR LA MUTACIÓN Y GUARDAR LAS BASES
			BaseRef[BaseActual] 	=	'0';
			BaseRead[BaseActual]	=	Read[Offset];

			for(int q=(L-1);q>Offset+1;q--){
				Read[q]	=	Read[q-1];
			}
			Read[Offset]	=	'N';

			//IMPRIMIR RESULTADOS
			fprintf(ALIGN,"Base reference %c\n",BaseRef[BaseActual]);
			fprintf(ALIGN,"Base read: %c\n",BaseRead[BaseActual]);
			fprintf(ALIGN,"SECUENCE: %s\n",Read);
		break;
		case 'T':
			//APLICAR LA MUTACIÓN Y GUARDAR LAS BASES
			BaseRef[BaseActual] 	=	Read[Offset];
			BaseRead[BaseActual]	=	'0';
			for(int q=Offset-2;q<L-3;q++){
				Read[q] = Read[q+3];
			}

			//IMPRIMIR RESULTADOS
			fprintf(ALIGN,"Base reference %c\n",BaseRef[BaseActual]);
			fprintf(ALIGN,"Base read: %c\n",BaseRead[BaseActual]);
			fprintf(ALIGN,"SECUENCE: %s\n",Read);
		break;
		case 'S':

			BasesAcum[0]		= 	BasesStats[0];
			for (int i=1; i<4; i++){
				BasesAcum[i]= BasesStats[i] + BasesAcum[i-1];
			}

			//CALCULAR LA BASE DESTINO
			int flag = 0;
			do{
				double dado		=	LanzarDado();
				int sel			=	BusqBin_Rul(BasesAcum,4,dado);
				char BaseDest	=	selBase(sel,Read[Offset]);
				if(BaseDest != Read[Offset-1]){
					//GUARDAR LAS BASES Y APLICAR LA MUTACIÓN
					BaseRef[BaseActual]		=	Read[Offset];
					BaseRead[BaseActual]	=	BaseDest;
					Read[Offset]	=	BaseDest;
					Read[Offset+1]	=	BaseDest;
					flag = 1;
				}
			}while(flag != 1);

			fprintf(ALIGN,"Bases reference %c\n",BaseRef[BaseActual]);
			fprintf(ALIGN,"Base read: %c\n",BaseRead[BaseActual]);
			fprintf(ALIGN,"SECUENCE: %s\n",Read);

		break;
		case 'C':
			//APLICAR LA MUTACIÓN Y GUARDAR LAS BASES
			BaseRef[BaseActual] 	=	Read[Offset];
			BaseRead[BaseActual]	=	'0';
			for(int q=Offset-3;q<L-4;q++){
				Read[q] = Read[q+4];
			}

			//IMPRIMIR RESULTADOS
			fprintf(ALIGN,"Base reference %c\n",BaseRef[BaseActual]);
			fprintf(ALIGN,"Base read: %c\n",BaseRead[BaseActual]);
			fprintf(ALIGN,"SECUENCE: %s\n",Read);
		break;
	}
 }
