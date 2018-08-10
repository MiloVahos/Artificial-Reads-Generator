/*
 ============================================================================
 Name        : 	ArtificialReadsGenerator.c
 Author      : 	Juan Camilo Peña Vahos - Aníbal Guerra - Sebastian Isaza Ramírez
 Version     :
 Copyright   : 	This project is totally opensource
 Description :	Generador de reads artificiales, partiendo de un archivo fasta
 	 	 	 	referencia para generar reads
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <inttypes.h>

#include "Stats.h"
#include "Matching.h"
#include "Mutation.h"
#include "FilesUtils.h"

#define READ_BIAS		100		//SE OBTIENEN ALGUNAS BASES EXTRAS PARA OPERAR
#define ERROR_PER   	0.25    //PORCENTAJE DE ERRORES EN EL READ
#define NAMES_SIZE		40		//ESPACIO PARA EL NOMBRE DE LOS ARCHIVOS

void FordwardMutation	(uint8_t,char*,uint16_t,uint8_t*,uint8_t*,int,int,FILE*);
int ReverseMutation		(uint8_t,char*,uint16_t,uint8_t*,uint8_t*,int,int,FILE*);

int main(int argc, char *argv[]) {

	//ARGUMENTOS DE ENTRADA
	char 	*DATA;					//NOMBRE DEL ARCHIVO
	char 	*I;						//IDENTIFICADOR DE LOS READS
	char	*Q;						//VALOR DEL QUALITY SCORE
	int 	L	=	0;				//LONGITUD DEL READ
	int		C	=	0;				//COVERAGE
	int 	B	=	0;				//BASE
	double 	lambda 	= 0;			//VALOR AJUSTABLE DE ENTRADA

	//VARIABLES DEL PROCESO PARA SALIDA
	char *MT;						//MATCHING TYPE

	//VARIABLES DEL PROCESO PARA OPERAR
	char		*Reference;			//REFERENCIA OBTENIDA DEL ARCHIVO FASTA
	char		*Read;				//READ
	char 		*RefName;			//NOMBRE DEL ARCHIVO DE REFERENCIA
	char 		*RefFastq;			//NOMBRE ARCHIVO FASTQ
	char 		*RefFastqseq;		//NOMBRE ARCHIVO FASTQSEQ
	char 		*RefAlign;			//NOMBRE ARCHIVO ALIGN
	char 		*RefMeta;			//NOMBRE ARCHIVO META
	uint32_t	TotalReads;			//TotalReads	=	BxC
	uint64_t 	TotalChars;			//Total de caracteres en la referencia
	FILE 		*FASTQ, *FASTQSEQ;	//PUNTEROS A LOS ARCHIVOS
	FILE		*ALIGN,	*META;		//PUNTEROS A LOS ARCHIVOS
	int 		MaxK;				//TOTAL DE ERRORES L*ERROR_PER
	int			t;					//NÚMERO DE ELEMENTOS DE LA RULETA
	double		dado;				//DADO

	//VARIABLES CON RELACIÓN AL PROCESO DE MUTACIÓN
	uint32_t  	id;             	 	//Identificador consecutivo para cada uno de los Reads
	uint32_t 	Pos;            		//Posición de Matching respecto a la referencia
	uint16_t  	lendesc;         		//Cantidad de errores total en el Read
	char      	strand;         		//Caractér con el sentido del matching
	uint8_t   	*Oper;          		//Arreglo con la operación por error
	uint16_t	*Cnt;					//Arreglo con los contadores por cada uno de los tipos de mutación
	uint16_t  	*Offsets;       		//Arreglo de offsets por cada error
	uint8_t   	*BaseRef;       		//Arreglo con la base de la referencia (Read Referencia)
	uint8_t   	*BaseRead;      		//Arreglo con la base después de la mutación (Read Destino)

	//Obtener los datos suministrados en la linea de comando
	if(argc>1){
		for (int i = 1; i < argc; i++)	{
			if (strcmp(argv[i], "-DATA")	== 0)	DATA	=	argv[i+1];
			if (strcmp(argv[i], "-I")		== 0)	I		=	"@EAS139:136:FC706VJ:2:2104:15343:197393 1:N:18:1";
			if (strcmp(argv[i], "-Q")		== 0)	Q		=	"Valor temporal de QQQQQQQ";
			if (strcmp(argv[i], "-L") 		== 0)	L		=	atoi(argv[i+1]);
			if (strcmp(argv[i], "-C") 		== 0)	C		=	atoi(argv[i+1]);
			if (strcmp(argv[i], "-B") 		== 0)	B		=	atoi(argv[i+1]);
			if (strcmp(argv[i], "-lambda") 	== 0)	lambda	=	atof(argv[i+1]);
		}
	}

	//OBTENER EL NOMBRE DE LA SECUENCIA Y CREAR EL NOMBRE DE LOS ARCHIVOS DE SALIDA
	RefName		=	(char*)	malloc(NAMES_SIZE*sizeof(char));
	if (RefName	==	NULL) printf ("Not enough memory for RefName");
	RefFastq	=	(char*) malloc(NAMES_SIZE*sizeof(char));
	if (RefFastq	==	NULL) printf ("Not enough memory for RefFastq");
	RefFastqseq	=	(char*) malloc(NAMES_SIZE*sizeof(char));
	if (RefFastqseq	==	NULL) printf ("Not enough memory for RefFastqseq");
	RefAlign	=	(char*) malloc(NAMES_SIZE*sizeof(char));
	if (RefAlign	==	NULL) printf ("Not enough memory for RefAlign");
	RefMeta		=	(char*) malloc(NAMES_SIZE*sizeof(char));
	if (RefMeta	==	NULL) printf ("Not enough memory for RefMeta");

	char *dot	=	strrchr(DATA,'.');
	strncpy(RefName,DATA,dot - DATA);

	//CREAR Y ABRIR LOS ARCHIVOS DE SALIDA
	strcpy(RefFastq,RefName);	strcpy(RefFastqseq,RefName);
	strcpy(RefAlign,RefName);	strcpy(RefMeta,RefName);
	FASTQ		=	fopen(strcat(RefFastq	,".fastq")		,"w");
	FASTQSEQ	=	fopen(strcat(RefFastqseq,".fastqseq")	,"w");
	ALIGN		=	fopen(strcat(RefAlign	,".align")		,"w");
	META		=	fopen(strcat(RefName	,".meta")		,"w");

	//SE LLENA EL ARCHIVO META
	fprintf(META,"Nombre del archivo FASTA:	%s\n",DATA);
	fprintf(META,"Longitud de los reads: %d\n",L);
	fprintf(META,"Valor del coverage: %d\n",C);
	fprintf(META,"Valor de la base:	%d\n",B);
	fprintf(META,"Valor de lambda:	%f\n",lambda);

	//ANTES DE COMENZAR, SE HALLA LOS VECTOR DE PROBABILIDADES

	//VECTOR PARA LOS MATCHING	{FORWARD, REVERSE, COMPLEMENT, REVERSE COMPLEMENT}
	double	MatTypeStats[MATCHING_TYPES]   =   {0.4,0.4,0.1,0.1};
	double 	MatTypeAcumF[MATCHING_TYPES];
	MatTypeAcumF[0]		= MatTypeStats[0];
	for (int i=1; i<MATCHING_TYPES; i++){
		MatTypeAcumF[i]	= MatTypeStats[i] + MatTypeAcumF[i-1];
	}

	//VECTOR DE PROBABILIDADES DE ERROR EXPONENCIALES
	MaxK	=   (int) floor(ERROR_PER * L);
	double	*ErrorStat	=	(double*)	malloc((MaxK+1)*sizeof(double));
	t	=   AdjustKExp(TINICIAL,lambda,MaxK);
	generarRuletaExp(ErrorStat,t,lambda);

	//VECTOR DE PROBABILIDADES DE LAS MUTACIONES
	//Simple Substitution, Single deletion, Insertion, Contiguos Deletion
	//N Insertions, Triple Contiguos deletion, Contiguos Repeated Substitution
	//Quadruple Contiguos deletion
	double		MutTypeStats[MUTATION_TYPES]   =   {0.63,0.15,0.071,0.065,0.049,0.006,0.0009,0.0001};
	double 		MutTypeAcumF[MUTATION_TYPES];
	double		total_adaptMut	= CalculaTotal(MUTATION_TYPES,MutTypeStats);
	MutTypeAcumF[0]		= MutTypeStats[0] / total_adaptMut;
	for (int i=1; i<MUTATION_TYPES; i++){
		MutTypeAcumF[i]= (MutTypeStats[i] / total_adaptMut) + MutTypeAcumF[i-1];
	}

	TotalChars	=	contChars(DATA);				//CANTIDAD DE CARACTERES DE LA REFERENCIA
	Reference	=	(char*) malloc(TotalChars*sizeof(char));
	if (Reference	==	NULL){
		printf ("Not enough memory for Reference, Mission Abort!!!");
		exit(1);
	}
	getReference(DATA,Reference);					//SE OBTIENE LA REFERENCIA
	TotalReads	=	B*C;							//NUMERO TOTAL DE READS A GENERAR

	srand(time(NULL));								//SEMILLA DE LOS ALEATOREOS

	for(uint32_t ReadsCicle = 0;	ReadsCicle<TotalReads; ReadsCicle++){

		MT	=	(char*) malloc(sizeof(char));

		//IDENTIFICADOR DEL READ
		id	=	ReadsCicle+1;
		fprintf(ALIGN,"Read ID:	%"PRIu32"\n",id);

		//GENERAR UNA POSICIÓN ALEATORIA EN EL RANGO [0,LengthRef-LengthRead]
		Pos		=	(rand() %((TotalChars-L-READ_BIAS) - 0 + 1)) + 0;
		fprintf(ALIGN,"Offset de la referencia: %"PRIu32"\n",Pos);

		//OBTENER EL READ DE REFERENCIA DESDE LA POSICIÓN DE MAPEO
		Read	=	(char*) malloc((L+READ_BIAS)*sizeof(char));
		memcpy(Read,Reference+Pos,L+READ_BIAS);
		fprintf(ALIGN,"Referencia sin matching: %s\n",Read);

		//CALCULAR EL MATCHING
		//ORDEN DE LOS ARREGLOS
		//FORWARD(F), REVERSE(R), COMPLEMENT(C), REVERSE COMPLEMENT(E)
		dado 	= 	LanzarDado();
		int MatTypeSel  =   BusqBin_Rul(MatTypeAcumF,MATCHING_TYPES,dado);
		selMatching(MatTypeSel,L,Read,MT);
		fprintf(ALIGN,"Referencia: %s\n",Read);

		//CALCULAR LA CANTIDAD DE ERROES
		dado	=	LanzarDado();
		lendesc	=	BusqBin_Rul(ErrorStat,t,dado);
		fprintf(ALIGN,"Cantidad de errores:	%"PRIu16"\n",lendesc);

		//AQUÍ HAY DOS POSIBILIDADES, QUE EXISTAN ERRORES Y QUE NO
		if(lendesc	!=	0){
			//EN ESTE CASO TENEMOS MUTACIONES EN EL READ
			Oper        =   (uint8_t*)  malloc(lendesc*sizeof(uint8_t));
			if (Oper	==	NULL) printf ("Not enough memory for Oper");
			Cnt			=	(uint16_t*) malloc(MUTATION_TYPES*sizeof(uint16_t));
			if (Cnt	==	NULL) printf ("Not enough memory for Cnt");
			Offsets		=   (uint16_t*) malloc(lendesc*sizeof(uint16_t));
			if (Offsets	==	NULL) printf ("Not enough memory for Offsets");

			for(int i=0; i<MUTATION_TYPES;i++) Cnt[i]	=	0;
			//Este acumulador determina cuanto correr el offset de acuerdo con
			//la cantidad de deleciones que vayan a ocurrir
			int operShift	=	0;

			//EL MATCHING SE PONE EN MINÚSCULA
			strand	=	(char) tolower(*MT);
			fprintf(ALIGN,"Matching type %c\n",strand);

			//SE VA A DETERMINAR EL VECTOR DE OPERACIONES PARA CALCULAR
			//LOS OFFSETS DE MANERA EFICIENTE
			//GENERAR EL VECTOR DE MUTACIONES A APLICAR
			for(int i=0;	i<lendesc;	i++){
				dado			=	LanzarDado();
				int ErrorSel	=	BusqBin_Rul(MutTypeAcumF,MUTATION_TYPES,dado);
				Oper[i]			=	selMutation(ErrorSel);
				switch((char)Oper[i]){
					//0->s 1->d 2->i 3->D 4->I 5->T 6->S 7->C
					case 's': Cnt[0]++; break;
					case 'd': Cnt[1]++; operShift++;  break;
					case 'i': Cnt[2]++; break;
					case 'D': Cnt[3]++; operShift+=2; break;
					case 'I': Cnt[4]++; break;
					case 'T': Cnt[5]++; operShift+=3; break;
					case 'S': Cnt[6]++; break;
					case 'C': Cnt[7]++; operShift+=4; break;
				}
			}

			//GENERAR EL VECTOR DE OFFSETS
			//LOS OFFSETS DEPENDE DEL TIPO DE MATCHING, DE LOS ERRORES QUE
			//SE VAYAN A APLICAR
			int OffsetAnt	=	0;	//ACUMULADOR DE OFFSETS
			int OffsetAux	=	0;	//OFFSET ACTUAL
			for (int i=0;	i<lendesc;	i++){
				if((strand=='f')||(strand=='c')){
					if(i==0){
						OffsetAux 	= 	rand() %(((L-1)-lendesc-operShift) + 1 - 0) + 0;
						Offsets[i]	=	OffsetAux;
					}else{
						switch(Oper[i-1]){
							case 's': case 'S': case 'i': case 'I':
								OffsetAux	=	rand() %(((L-1)-(lendesc-i)-operShift) + 1 - (OffsetAnt+1)) + (OffsetAnt+1);
								Offsets[i]	=	OffsetAux	-	OffsetAnt;
							break;
							default:
								OffsetAux	=	rand() %(((L-1)-(lendesc-i)-operShift) + 1 - (OffsetAnt)) + (OffsetAnt);
								Offsets[i]	=	OffsetAux	-	OffsetAnt;
						}
					}
				}else{
					if(i==0){
						OffsetAux = rand() %((L-1) + 1 - (0+lendesc+operShift)) + (0+lendesc+operShift);
						Offsets[i]	=	OffsetAux;
					}else{
						switch(Oper[i-1]){
							case 's': case 'S': case 'd':
								OffsetAux = rand() %((OffsetAnt-1) + 1 - (0+(lendesc-i)+operShift)) + (0+(lendesc-i)+operShift);
								Offsets[i]	=	OffsetAnt - OffsetAux;
							break;
							case 'D':
								OffsetAux = rand() %((OffsetAnt-2) + 1 - (0+(lendesc-i)+operShift)) + (0+(lendesc-i)+operShift);
								Offsets[i]	=	OffsetAnt - OffsetAux;
							break;
							case 'T':
								OffsetAux = rand() %((OffsetAnt-3) + 1 - (0+(lendesc-i)+operShift)) + (0+(lendesc-i)+operShift);
								Offsets[i]	=	OffsetAnt - OffsetAux;
							break;
							case 'C':
								OffsetAux = rand() %((OffsetAnt-4) + 1 - (0+(lendesc-i)+operShift)) + (0+(lendesc-i)+operShift);
								Offsets[i]	=	OffsetAnt - OffsetAux;
							break;
							default:
								OffsetAux = rand() %((OffsetAnt) + 1 - (0+(lendesc-i)+operShift)) + (0+(lendesc-i)+operShift);
								Offsets[i]	=	OffsetAnt - OffsetAux;
						}
					}
					switch(Oper[i]){
						case 'd': operShift--;  break;
						case 'D': operShift-=2; break;
						case 'T': operShift-=3; break;
						case 'C': operShift-=4; break;
					}
				}
				printf("%d ",OffsetAux);
				OffsetAnt	=	OffsetAux;
			}
			printf("\n");
			for(int i=0; i<lendesc; i++){
				printf("%c ",Oper[i]);
			}
			printf("matching : %c =>",strand);
			for(int i=0; i<lendesc; i++){
				printf("%d ",Offsets[i]);
			}
			printf("\n");
			int BaseActual	=	0;

			/*for(int i=0; i<lendesc;	i++){

				//SE INICIALIZAN LOS ARREGLOS

				BaseRef		=	(uint8_t*)  malloc((lendesc+operShift)*sizeof(uint8_t));
				if (BaseRef	==	NULL) printf ("Not enough memory for BaseRef");
				BaseRead	=	(uint8_t*)  malloc((lendesc+operShift)*sizeof(uint8_t));
				if (BaseRead ==	NULL) printf ("Not enough memory for BaseRead");

				fprintf(ALIGN,"Operación de mutación: %c\n",Oper[i]);
				//EL OFFSET SE GENERA DE ACUERDO AL TIPO DE MATCHING, FORWARD OR REVERSE
				if ((strand	==	'f')||(strand	==	'c')) {
					//FORWARD MATCHINGS
					if(i==0){
						OffsetAux 	= 	rand() %(((L-1)-lendesc-operShift) + 1 - 0) + 0;
						Offsets[i]	=	OffsetAux;
					}else{
						OffsetAux	=	rand() %(((L-1)-(lendesc-i)-operShift) + 1 - (OffsetAnt+1)) + (OffsetAnt+1);
						Offsets[i]	=	OffsetAux	-	OffsetAnt;
					}
					fprintf(ALIGN,"Offset[%d]	=	%d\n",i,Offsets[i]);
					OffsetAnt	=	OffsetAux;
					//APLICAR LA MUTACIÓN
					FordwardMutation(Oper[i],Read,OffsetAux,BaseRef,BaseRead,L,BaseActual,ALIGN);
					BaseActual++;

				}else{
					//REVERSE MATCHINGS
					if(i==0){
						OffsetAux2 = rand() %((L-1) + 1 - (0+lendesc+operShift)) + (0+lendesc+operShift);
						Offsets[i]	=	OffsetAux2;
					}else{
						OffsetAux2 = rand() %((OffsetAnt2-1) + 1 - (0+(lendesc-i)+operShift)) + (0+(lendesc-i)+operShift);
						Offsets[i]	=	OffsetAnt2 - OffsetAux2;
					}
					fprintf(ALIGN,"Offset[%d]	=	%"PRIu16"\n",i,Offsets[i]);
					OffsetAnt2	=	OffsetAux2;

					//APLICAR LA MUTACIÓN
					//BaseActual = ReverseMutation(Oper[i],Read,OffsetAux2,BaseRef,BaseRead,L,BaseActual,ALIGN);

				}

				//ACTUALIZAR LOS PARÁMETRO DEL CÁLCULO DEL OFFSET

				if(BaseRef)		free(BaseRef);
				if(BaseRead)	free(BaseRead);
			}*/

			//IMPRIMIR LOS CONTADORES
			printCounters(ALIGN,Cnt);
			//GENERAR EL READ
			generateRead(Read,id,L,Q,I,FASTQ,FASTQSEQ);
			printf("\n");
			if(Oper)	free(Oper);
			if(Cnt)		free(Cnt);
			if(Offsets)		free(Offsets);
		}else{
			//EN ESTE CASO NO HAY ERRORES
			//COMO NO HAY ERRORES EL MATCHING SE REPRESENTA EN MAYÚSCULA Y ES PERFECTO
			strand	=	*MT;
			//EN EL ARCHIVO DE ALINEACIÓN SE MUESTRA EL TIPO DE MATCH NADA MÁS
			fprintf(ALIGN,"Match Perfecto del tipo %c\n",strand);
			//SE IMPRE EL READ Y LA SECUENCIA
			generateRead(Read,id,L,Q,I,FASTQ,FASTQSEQ);
		}

		fprintf(ALIGN,"\n\n");
		fprintf(FASTQ,"\n\n");

		if(Read)	free(Read);
		if(MT)		free(MT);
	}


	fclose(META);
	fclose(ALIGN);
	fclose(FASTQ);
	fclose(FASTQSEQ);

	if(Reference)	free(Reference);
	if(RefName)		free(RefName);
	if(RefMeta)		free(RefMeta);
	if(RefAlign)	free(RefAlign);
	if(RefFastq)	free(RefFastq);
	if(RefFastqseq) free(RefFastqseq);
	if(ErrorStat)   free(ErrorStat);
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

			printf("Offset	=	%d, BaseActual	=	%d\n",Offset,BaseActual);
			//APLICAR LA MUTACIÓN Y GUARDAR LAS BASES
			BaseRef[BaseActual] 	=	Read[Offset];
			BaseRead[BaseActual]	=	BaseDest;
			Read[Offset]	=	BaseDest;

			//IMPRIMIR RESULTADOS
			fprintf(ALIGN,"Base reference %c\n",BaseRef[BaseActual]);
			fprintf(ALIGN,"Base read: %c\n",BaseRead[BaseActual]);
			fprintf(ALIGN,"SECUENCE: %s\n",Read);

		break;
		/*case 'd':

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
			BaseRef[BaseActual+1] 	=	Read[Offset+1];
			BaseRead[BaseActual+1]	=	'0';
			for(int q=Offset;q<L-2;q++){
				Read[q] = Read[q+2];
			}

			//IMPRIMIR RESULTADOS
			fprintf(ALIGN,"Base reference %c%c\n",BaseRef[BaseActual],BaseRef[BaseActual+1]);
			fprintf(ALIGN,"Base read: %c%c\n",BaseRead[BaseActual],BaseRead[BaseActual+1]);
			fprintf(ALIGN,"SECUENCE: %s\n",Read);

		break;
		case 'I':

			for(int q=(L-1);q>Offset;q--){
				Read[q]	=	Read[q-1];
			}
			Read[Offset]	=	'N';
			//APLICAR LA MUTACIÓN Y GUARDAR LAS BASES
			BaseRef[BaseActual] 	=	'0';
			BaseRead[BaseActual]	=	Read[Offset];

			//IMPRIMIR RESULTADOS
			fprintf(ALIGN,"Base reference %c\n",BaseRef[BaseActual]);
			fprintf(ALIGN,"Base read: %c\n",BaseRead[BaseActual]);
			fprintf(ALIGN,"SECUENCE: %s\n",Read);

		break;
		case 'T':
			//APLICAR LA MUTACIÓN Y GUARDAR LAS BASES
			BaseRef[BaseActual] 	=	Read[Offset];
			BaseRead[BaseActual]	=	'0';
			BaseRef[BaseActual+1] 	=	Read[Offset+1];
			BaseRead[BaseActual+1]	=	'0';
			BaseRef[BaseActual+2] 	=	Read[Offset+2];
			BaseRead[BaseActual+2]	=	'0';
			for(int q=Offset;q<L-3;q++){
				Read[q] = Read[q+3];
			}

			//IMPRIMIR RESULTADOS
			fprintf(ALIGN,"Base reference %c%c%c\n",BaseRef[BaseActual],BaseRef[BaseActual+1],BaseRef[BaseActual+2]);
			fprintf(ALIGN,"Base read: %c%c%c\n",BaseRead[BaseActual],BaseRead[BaseActual+1],BaseRead[BaseActual+2]);
			fprintf(ALIGN,"SECUENCE: %s\n",Read);

		break;*/
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


			fprintf(ALIGN,"Bases reference %c%c\n",BaseRef[BaseActual],BaseRef[BaseActual+1]);
			fprintf(ALIGN,"Base read: %c%c\n",BaseRead[BaseActual],BaseRead[BaseActual+1]);
			fprintf(ALIGN,"SECUENCE: %s\n",Read);

		break;
		/*case 'C':
			//APLICAR LA MUTACIÓN Y GUARDAR LAS BASES
			BaseRef[BaseActual] 	=	Read[Offset];
			BaseRead[BaseActual]	=	'0';
			BaseRef[BaseActual+1] 	=	Read[Offset+1];
			BaseRead[BaseActual+1]	=	'0';
			BaseRef[BaseActual+2] 	=	Read[Offset+2];
			BaseRead[BaseActual+2]	=	'0';
			BaseRef[BaseActual+3] 	=	Read[Offset+3];
			BaseRead[BaseActual+3]	=	'0';
			for(int q=Offset;q<L-4;q++){
				Read[q] = Read[q+4];
			}

			//IMPRIMIR RESULTADOS
			fprintf(ALIGN,"Base reference %c%c%c%c\n",BaseRef[BaseActual],BaseRef[BaseActual+1],BaseRef[BaseActual+2],BaseRef[BaseActual+3]);
			fprintf(ALIGN,"Base read: %c%c%c%c\n",BaseRead[BaseActual],BaseRead[BaseActual+1],BaseRead[BaseActual+2],BaseRead[BaseActual+3]);
			fprintf(ALIGN,"SECUENCE: %s\n",Read);

		break;*/
	}
 }

int ReverseMutation(uint8_t Oper,char *Read,uint16_t Offset, uint8_t *BaseRef,uint8_t *BaseRead, int L, int BaseActual, FILE *ALIGN){
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

			return (BaseActual+1);
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

			return (BaseActual+1);
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

			return (BaseActual+1);
		break;
		case 'D':

			//APLICAR LA MUTACIÓN Y GUARDAR LAS BASES
			BaseRef[BaseActual] 	=	Read[Offset];
			BaseRead[BaseActual]	=	'0';
			BaseRef[BaseActual+1] 	=	Read[Offset-1];
			BaseRead[BaseActual+1]	=	'0';
			for(int q=Offset-1;q<L-2;q++){
				Read[q] = Read[q+2];
			}

			//IMPRIMIR RESULTADOS
			fprintf(ALIGN,"Base reference %c%c\n",BaseRef[BaseActual],BaseRef[BaseActual+1]);
			fprintf(ALIGN,"Base read: %c%c\n",BaseRead[BaseActual],BaseRead[BaseActual+1]);
			fprintf(ALIGN,"SECUENCE: %s\n",Read);

			return (BaseActual+2);
		break;
		case 'I':

			for(int q=(L-1);q>Offset;q--){
				Read[q]	=	Read[q-1];
			}
			Read[Offset]	=	'N';
			//APLICAR LA MUTACIÓN Y GUARDAR LAS BASES
			BaseRef[BaseActual] 	=	'0';
			BaseRead[BaseActual]	=	Read[Offset];

			//IMPRIMIR RESULTADOS

			fprintf(ALIGN,"Base reference %c\n",BaseRef[BaseActual]);
			fprintf(ALIGN,"Base read: %c\n",BaseRead[BaseActual]);
			fprintf(ALIGN,"SECUENCE: %s\n",Read);

			return (BaseActual+1);
		break;
		case 'T':
			//APLICAR LA MUTACIÓN Y GUARDAR LAS BASES
			BaseRef[BaseActual] 	=	Read[Offset];
			BaseRead[BaseActual]	=	'0';
			BaseRef[BaseActual+1] 	=	Read[Offset-1];
			BaseRead[BaseActual+1]	=	'0';
			BaseRef[BaseActual+2] 	=	Read[Offset-2];
			BaseRead[BaseActual+2]	=	'0';
			for(int q=Offset-2;q<L-3;q++){
				Read[q] = Read[q+3];
			}

			//IMPRIMIR RESULTADOS
			fprintf(ALIGN,"Base reference %c%c%c\n",BaseRef[BaseActual],BaseRef[BaseActual+1],BaseRef[BaseActual+2]);
			fprintf(ALIGN,"Base read: %c%c%c\n",BaseRead[BaseActual],BaseRead[BaseActual+1],BaseRead[BaseActual+2]);
			fprintf(ALIGN,"SECUENCE: %s\n",Read);

			return (BaseActual+2);
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
					BaseRef[BaseActual+1]	=	Read[Offset-1];
					BaseRead[BaseActual]	=	BaseDest;
					BaseRead[BaseActual+1]	=	BaseDest;
					Read[Offset]	=	BaseDest;
					Read[Offset+1]	=	BaseDest;
					flag = 1;
				}
			}while(flag != 1);

			fprintf(ALIGN,"Bases reference %c%c\n",BaseRef[BaseActual],BaseRef[BaseActual+1]);
			fprintf(ALIGN,"Base read: %c%c\n",BaseRead[BaseActual],BaseRead[BaseActual+1]);
			fprintf(ALIGN,"SECUENCE: %s\n",Read);

			return (BaseActual+2);
		break;
		case 'C':
			//APLICAR LA MUTACIÓN Y GUARDAR LAS BASES
			BaseRef[BaseActual] 	=	Read[Offset];
			BaseRead[BaseActual]	=	'0';
			BaseRef[BaseActual+1] 	=	Read[Offset-1];
			BaseRead[BaseActual+1]	=	'0';
			BaseRef[BaseActual+2] 	=	Read[Offset-2];
			BaseRead[BaseActual+2]	=	'0';
			BaseRef[BaseActual+3] 	=	Read[Offset-3];
			BaseRead[BaseActual+3]	=	'0';
			for(int q=Offset-3;q<L-4;q++){
				Read[q] = Read[q+4];
			}

			//IMPRIMIR RESULTADOS
			fprintf(ALIGN,"Base reference %c%c%c%c\n",BaseRef[BaseActual],BaseRef[BaseActual+1],BaseRef[BaseActual+2],BaseRef[BaseActual+3]);
			fprintf(ALIGN,"Base read: %c%c%c%c\n",BaseRead[BaseActual],BaseRead[BaseActual+1],BaseRead[BaseActual+2],BaseRead[BaseActual+3]);
			fprintf(ALIGN,"SECUENCE: %s\n",Read);

			return (BaseActual+2);
		break;
	}
 }
