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

//LIBRERÍAS
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

#define READ_BIAS	1024	//SE OBTIENEN ALGUNAS BASES EXTRAS PARA OPERAR
#define ERROR_PER   0.25    //PORCENTAJE DE ERRORES EN EL READ
#define NAMES_SIZE	40		//ESPACIO PARA EL NOMBRE DE LOS ARCHIVOS
#define READID_SIZE 50		//TAMAÑO DE I
#define READQ_SIZE 	1024	//TAMAÑO DE Q

int main(int argc, char *argv[]) {

	//ARGUMENTOS DE ENTRADA
	char 	*DATA;		//NOMBRE DEL ARCHIVO
	char 	*I;			//IDENTIFICADOR DE LOS READS
	char	*Q;			//VALOR DEL QUALITY SCORE
	uint16_t L	=	0;	//LONGITUD DEL READ
	uint8_t	 C	=	0;	//COVERAGE
	uint32_t B	=	0;	//BASE
	double 	 P0 = 	0;	//VALOR AJUSTABLE DE ENTRADA

	//VARIABLES DEL PROCESO PARA SALIDA
	char *MT;	//MATCHING TYPE

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
	uint32_t  	id; 		//Identificador consecutivo para cada uno de los Reads
	uint32_t 	Pos;        //Posición de Matching respecto a la referencia
	uint16_t  	lendesc;    //Cantidad de errores total en el Read
	char      	strand;     //Caractér con el sentido del matching
	uint8_t   	*Oper;      //Arreglo con la operación por error
	uint16_t	*Cnt;		//Arreglo con los contadores por cada uno de los tipos de mutación
	uint16_t  	*Offsets;   //Arreglo de offsets por cada error
	uint8_t   	*BaseRef;   //Arreglo con la base de la referencia (Read Referencia)
	uint8_t   	*BaseRead;  //Arreglo con la base después de la mutación (Read Destino)

	//Obtener los datos suministrados en la linea de comando
	if(argc>1){
		for (int i = 1; i < argc; i++)	{
			if (strcmp(argv[i], "-DATA")	== 0){
				DATA	=	(char*) malloc(NAMES_SIZE*sizeof(char));
				if (DATA ==	NULL) printf ("Not enough memory for DATA"); 
				strcpy(DATA,argv[i+1]);
			}	
			if (strcmp(argv[i], "-I")		== 0){
				I	=	(char*) malloc(READID_SIZE*sizeof(char));
				if (I	==	NULL) printf ("Not enough memory for I"); 
				strcpy(I,"@EAS139:136:FC706VJ:2:2104:15343:197393 1:N:18:190");
			}	
			if (strcmp(argv[i], "-Q")		== 0){
				Q	= 	(char*) malloc(READQ_SIZE*sizeof(char));
				if (Q	==	NULL) printf ("Not enough memory for Q"); 
				strcpy(Q,"3eEKj]-6LgGAE84dgkck14J0Gh[LK0d9jjC5Lk]JbD]cfe0I-6G6KcL2-I7k-kE98EK3hJhjLa^Bf4A9i3C6i6IaJc3a1L78Ai31a738Acj372FDAIedE98Ff*EB3B2iE2L5hajfcfAAFgL57l3-53]039HJ5JGEBc1E]HGJi[dK[lf2*_4B6Kc6K4D8DGLk_h]8jekHBiJ0H816BJ[fCe1]DLBI2]dB]1aHa]3kbK2CF^eifeKa3E*LH-_CJibb[^aBffKFB7dIH3Jgi6hLGjlJa5**l4bJclJaFf065f0be8heIG9j_lABK_LC568KhDFgBK]3Fc61bk_5HF5A7*2-l9-gCfA351l5B67ICh2gEJDG65eJa*LEK1Bk_[EgGHL1b0[]LJ222Id_hDi]K-Ea^l-BAl4ecbEH9k0h-bB[hL8^02e05B7fj]C1lK^ebaBcG3[clC]hKf]797F36C___i2igcJfcKh]4-5_5k[Cflc6j2[b11[E^-9k]5c-jDjED^H6fgbBhEH07B3BJD[d0cAj0^AEeklI*3[j1lg-JK0j6-5d-EjicGk3Lk6CBfLGL7]g80d97l97*bhbA2b]9H-fJCk8[d_G3^e*89kCD_cdjhkH^lg]LBK1-G0c_c*2KaHk7KDEfgkeD6aB61cEcIja[IeGbJ[0g-ghhI*A82[cFkJ40Hh]-7Hi0Fh^K47CfCbCI9La74bK3g_[bBH4^BfgkF8-^1dfj1FCA15af95B*8k^]5AJEDf27H6Kb^i1KHBAHic1bkIKal]B58-^HfA2kgI9_Hg^^4HB^LiBh]1F_8K*C^KJHgcI^3JhAeki[0_0AJec469h]5A2gCl9hA56f8a3c31JfKg9lhC33lKH91C2j4bC^f3Id2kffFB*6K3EH_fGH157lKKJkG8Ll6]gL91^3i_ik-EB4d1BB_8-fIG*kIgk^-c1H34gIa9B[eCHD[*H90A2AkFjli3d2dKJD1^FjAEH3CBh1g[h168F88bbdEH8-_hK]jj[");
			}
			if (strcmp(argv[i], "-L") 		== 0)	L		=	(uint16_t) atoi(argv[i+1]);
			if (strcmp(argv[i], "-C") 		== 0)	C		=	(uint8_t)  atoi(argv[i+1]);
			if (strcmp(argv[i], "-B") 		== 0)	B		=	(uint32_t) atoi(argv[i+1]);
			if (strcmp(argv[i], "-P0") 		== 0)	P0		=	atof(argv[i+1]);
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
	fprintf(META,"Longitud de los reads: %"PRIu16"\n",L);
	fprintf(META,"Valor del coverage: %"PRIu8"\n",C);
	fprintf(META,"Valor de la base:	%"PRIu32"\n",B);
	fprintf(META,"Valor de lambda:	%f\n",P0);

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
	double	*ErrorStat;
	ErrorStat	=	(double*) malloc((MaxK+1)*sizeof(double));
	if (ErrorStat == NULL) printf ("Not enough memory for ErrorStat"); 
	t	=   AdjustKExp(TINICIAL,P0,MaxK);
	generarRuletaExp(ErrorStat,t,P0);

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

	TotalChars	=	contChars(DATA);	//CANTIDAD DE CARACTERES DE LA REFERENCIA
	Reference	=	(char*) malloc(TotalChars*sizeof(char));
	if (Reference	==	NULL){
		printf ("Not enough memory for Reference, Mission Abort!!!");
		exit(1);
	}
	getReference(DATA,Reference);					//SE OBTIENE LA REFERENCIA
	TotalReads	=	B*C;							//NUMERO TOTAL DE READS A GENERAR
	srand(time(NULL));								//SEMILLA DE LOS ALEATOREOS

	for(uint32_t ReadsCicle = 0; ReadsCicle<TotalReads; ReadsCicle++){

		MT	=	(char*) malloc(sizeof(char));
		if (MT == NULL) printf ("Not enough memory for MT"); 

		//IDENTIFICADOR DEL READ
		id	=	ReadsCicle+1;
		fprintf(ALIGN,"ID: %"PRIu32"\n",id);

		//GENERAR UNA POSICIÓN ALEATORIA EN EL RANGO [0,LengthRef-LengthRead]
		Pos		=	(rand() %((TotalChars-L-READ_BIAS) - 0 + 1)) + 0;
		fprintf(ALIGN,"Offset Ref: %"PRIu32"\n",Pos);

		//OBTENER EL READ DE REFERENCIA DESDE LA POSICIÓN DE MAPEO
		Read	=	(char*) malloc((L+READ_BIAS)*sizeof(char));
		if (Read == NULL) printf ("Not enough memory for Read"); 		
		memcpy(Read,Reference+Pos,L+READ_BIAS);
		//fprintf(ALIGN,"Referencia sin matching: %s\n",Read);

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
		/*if(lendesc	!=	0){
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
				OffsetAnt	=	OffsetAux;
			}

			//SE INICIALIZAN LOS ARREGLOS
			BaseRef		=	(uint8_t*)  malloc((lendesc)*sizeof(uint8_t));
			if (BaseRef	==	NULL) printf ("Not enough memory for BaseRef");
			BaseRead	=	(uint8_t*)  malloc((lendesc)*sizeof(uint8_t));
			if (BaseRead ==	NULL) printf ("Not enough memory for BaseRead");

			int BaseActual	=	0;
			OffsetAnt	=	0;

			for(int i=0; i<lendesc;	i++){

				fprintf(ALIGN,"Operación de mutación: %c\n",Oper[i]);
				//EL OFFSET SE GENERA DE ACUERDO AL TIPO DE MATCHING, FORWARD OR REVERSE
				if ((strand	==	'f')||(strand	==	'c')) {
					//FORWARD MATCHINGS
					fprintf(ALIGN,"Offset[%d]	=	%d\n",i,Offsets[i]);
					if(i==0){
						//APLICAR LA MUTACIÓN
						//printf("Oper = %c, Offset = %d, L = %d, Base = %d\n",Oper[i],Offsets[0],L,BaseActual);
						FordwardMutation(Oper[i],Read,Offsets[0],BaseRef,BaseRead,L,BaseActual,ALIGN);
						OffsetAnt	=	Offsets[0];
					}else{
						//APLICAR LA MUTACIÓN
						//printf("Oper = %c, Offset = %d, L = %d, Base = %d\n",Oper[i],OffsetAnt+Offsets[i],L,BaseActual);
						FordwardMutation(Oper[i],Read,(OffsetAnt+Offsets[i]),BaseRef,BaseRead,L,BaseActual,ALIGN);
						OffsetAnt	= OffsetAnt+Offsets[i];
					}
				}else{
					//REVERSE MATCHINGS
					fprintf(ALIGN,"Offset[%d]	=	%d\n",i,Offsets[i]);
					if(i==0){
						//APLICAR LA MUTACIÓN
						ReverseMutation(Oper[i],Read,Offsets[0],BaseRef,BaseRead,L,BaseActual,ALIGN);
						OffsetAnt	=	Offsets[0];
					}else{
						//APLICAR LA MUTACIÓN
						ReverseMutation(Oper[i],Read,(OffsetAnt-Offsets[i]),BaseRef,BaseRead,L,BaseActual,ALIGN);
						OffsetAnt	= OffsetAnt-Offsets[i];
					}
				}
				BaseActual++;
			}
			if(BaseRef)		free(BaseRef);
			if(BaseRead)	free(BaseRead);

			//IMPRIMIR LOS CONTADORES
			printCounters(ALIGN,Cnt);
			//GENERAR EL READ
			generateRead(Read,id,L,Q,I,FASTQ,FASTQSEQ);
			if(Oper)	free(Oper);
			if(Cnt)		free(Cnt);
			if(Offsets)	free(Offsets);
		}else{
			//EN ESTE CASO NO HAY ERRORES
			//COMO NO HAY ERRORES EL MATCHING SE REPRESENTA EN MAYÚSCULA Y ES PERFECTO
			strand	=	*MT;
			//EN EL ARCHIVO DE ALINEACIÓN SE MUESTRA EL TIPO DE MATCH NADA MÁS
			fprintf(ALIGN,"Match Perfecto del tipo %c\n",strand);
			//SE IMPRE EL READ Y LA SECUENCIA
			generateRead(Read,id,L,Q,I,FASTQ,FASTQSEQ);
		}*/

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
	if(ErrorStat)   free(ErrorStat);
	if(RefName)		free(RefName);
	if(RefMeta)		free(RefMeta);
	if(RefAlign)	free(RefAlign);
	if(RefFastq)	free(RefFastq);
	if(RefFastqseq) free(RefFastqseq);
	if(DATA)		free(DATA);
	if(I)			free(I);
	if(Q)			free(Q);
}