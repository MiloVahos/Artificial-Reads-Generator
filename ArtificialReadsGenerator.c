/*
 *  @Developers:    Juan Camilo Peña Vahos	- Aníbal Guerra	- Sebastian Isaza Ramirez
 *  @Last Revised:  17/07/2018
 *  @Description:   Generador artificial de reads para simular casos de prueba del compresor.
 *
*/

//LIBRERÍAS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <inttypes.h>

#include "Mutations.h"
#include "Matchings.h"
#include "Stats.h"

//DECLARACIÓN DE CONSTANTES

#define MUTATION_TYPES	8		//NÚMERO DE MUTACIONES POSIBLES
#define READ_BIAS		50		//SE OBTIENEN ALGUNAS BASES EXTRAS PARA OPERAR 
#define ERROR_PER   	0.25    //PORCENTAJE DE ERRORES EN EL READ

//DECLARACIÓN DE FUNCIONES
long int 	contChars(char*);
void 		getReference(char*,char*);
void		selMatching(int,int,char*,char*);
uint8_t		selMutation(int);
char 		selBase(int,char);
void 		FordwardMutation(uint8_t,char*,uint16_t,uint8_t*,uint8_t*,uint32_t);

int main (int argc, char *argv[]) {	
    
	//ARGUMENTOS DE ENTRADA
	char 	*DATA;					//NOMBRE DEL ARCHIVO
	char 	*I;						//IDENTIFICADOR DE LOS READS
	char 	*Q;						//VALOR DEL QUALITY SCORE
	int 	L	=	0;				//LONGITUD DEL READ
	int		C	=	0;				//COVERAGE
	int 	B	=	0;				//BASE
	int 	E	=	0;				//CANTIDAD EXACTA DE ERRORES
	int 	P	=	0;				//BANDERA DE EVENTO SOLO FORWARD AND REVERSE
	double 	lambda 	= 0;			//VALOR AJUSTABLE DE ENTRADA


	//VARIABLES DEL PROCESO PARA SALIDA
	char	*MT		=	(char*) malloc(sizeof(char)); //MATCHING TYPE
	char	*OTy	=	(char*) malloc(sizeof(char)); //OPERATION (MUTATION) TYPE

	//VARIABLES DEL PROCESO PARA OPERAR
	char		*Reference		=	NULL;				//REFERENCIA OBTENIDA DEL ARCHIVO FASTA
	char		*Read			=	NULL;				//READ
	char 		*RefName		=	NULL;				//NOMBRE DEL ARCHIVO DE REFERENCIA
	char 		*RefFastq		=	NULL;				//NOMBRE ARCHIVO FASTQ
	char 		*RefFastqseq	=	NULL;				//NOMBRE ARCHIVO FASTQSEQ
	char 		*RefAlign		=	NULL;				//NOMBRE ARCHIVO ALIGN
	char 		*RefMeta		=	NULL;				//NOMBRE ARCHIVO META
	int			TotalReads;								//TotalReads	=	BxC
	long int 	TotalChars;								//Total de caracteres en la referencia
	FILE 		*FASTQ, *FASTQSEQ,	*ALIGN,	*META;		//PUNTEROS A LOS ARCHIVOS
	int 		MaxK;									//TOTAL DE ERRORES L*ERROR_PER
	int			t;										//NÚMERO DE ELEMENTOS DE LA RULETA
	double		dado;									//DADO

	//VARIABLES CON RELACIÓN AL PROCESO DE MUTACIÓN
	uint32_t 	Pos;            		//Posición de Matching respecto a la referencia
    char      	strand;         		//Caractér con el sentido del matching
    uint16_t  	*Offsets;       		//Arreglo de offsets por cada error 
    uint8_t   	*Oper;          		//Arreglo con la operación por error   
    uint8_t   	*BaseRef;       		//Arreglo con la base de la referencia (Read Referencia)
    uint8_t   	*BaseRead;      		//Arreglo con la base después de la mutación (Read Destino) 
    uint32_t  	id;             	 	//Identificador consecutivo para cada uno de los Reads
    uint16_t  	lendesc;         		//Cantidad de errores total en el Read
    uint16_t  	Cnts,CntS,Cnti,CntI;
	uint16_t  	Cntd,CntD,CntT,CntC;    //Cantidad de veces que ocurre una mutación

	//Obtener los datos suministrados en la linea de comando
	if(argc>1){
		for (int i = 1; i < argc; i++)	{
			if (strcmp(argv[i], "-DATA")	== 0)	DATA	=	argv[i+1];
			if (strcmp(argv[i], "-I")    	== 0)	I		=	argv[i+1];
			if (strcmp(argv[i], "-Q") 		== 0)	Q		=	argv[i+1];	
			if (strcmp(argv[i], "-L") 		== 0)	L		=	atoi(argv[i+1]);
			if (strcmp(argv[i], "-C") 		== 0)	C		=	atoi(argv[i+1]);
			if (strcmp(argv[i], "-B") 		== 0)	B		=	atoi(argv[i+1]);
			if (strcmp(argv[i], "-E") 		== 0)	E		=	atoi(argv[i+1]);
			if (strcmp(argv[i], "-P") 		== 0)	P		=	atoi(argv[i+1]);
			if (strcmp(argv[i], "-lambda") 	== 0)	lambda	=	atof(argv[i+1]);	
		}
	}

	//OBTENER EL NOMBRE DE LA SECUENCIA Y CREAR EL NOMBRE DE LOS ARCHIVOS DE SALIDA
	RefName		=	(char*)	malloc(40*sizeof(char));
	RefFastq	=	(char*) malloc(40*sizeof(char));
	RefFastqseq	=	(char*) malloc(40*sizeof(char));
	RefAlign	=	(char*) malloc(40*sizeof(char));
	RefMeta		=	(char*) malloc(40*sizeof(char));

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
	fprintf(META,"Identificador de los reads: %s\n",I);
	fprintf(META,"Valor de los puntajes de calidad:	%s\n",Q);
	fprintf(META,"Longitud de los reads: %d\n",L);
	fprintf(META,"Valor del coverage: %d\n",C);
	fprintf(META,"Valor de la base:	%d\n",B);
	fprintf(META,"Valor de lambda:	%f\n",lambda);

	if(E!=0){	fprintf(META,"Número predeterminado de errores:	%d\n",E);	}
	else{		fprintf(META,"Sin número predeterminado de errores\n");		}
	if(P==1){	fprintf(META,"Prueba especial, solo se consideran matchings tipo Forward and Reverse\n");	}
	else{		fprintf(META,"Se consideran los cuatro tipos de matching\n");		}

	/**
	 * POR EL MOMENTO ASUMIREMOS QUE E Y P NUNCA OCURREN
	*/
	 
	//ANTES DE COMENZAR, SE HALLA EL VECTOR DE PROBABILIDADES EXPONENCIALES Y DEMÁS VECTORES
	//VECTOR PARA LOS MATCHING	{FORWARD, REVERSE, COMPLEMENT, REVERSE COMPLEMENT}
	double		MatTypeStats[MATCHING_TYPES]   =   {0.4,0.4,0.1,0.1}; 
	double 		MatTypeAcumF[MATCHING_TYPES];
	double 		total_adaptMat	=	CalculaTotal(MATCHING_TYPES,MatTypeStats);
	MatTypeAcumF[0]		= MatTypeStats[0] / total_adaptMat;
	for (int i=1; i<MATCHING_TYPES; i++){
		MatTypeAcumF[i]= (MatTypeStats[i] / total_adaptMat) + MatTypeAcumF[i-1];
	}

	//VECTOR DE PROBABILIDADES DE ERROR
	MaxK    =   (int)floor(ERROR_PER	*	L);
	double		*ErrorStat	=	(double*)	malloc(MaxK*sizeof(double));
    t		=   AdjustKExp(TINICIAL,lambda,MaxK);
	generarRuletaExp(ErrorStat,t,lambda);


	//VECTOR DE PROBABILIDADES DE LAS MUTACIONES
	//Simple Substitution, Single deletion, Insertion, Contiguos Deletion
	//N Insertions, Triple Contiguos deletion, Contiguos Repeated Substitution
	//Quadruple Contiguos deletion
	//double		MutTypeStats[MUTATION_TYPES]   =   {0.63,0.15,0.071,0.065,0.049,0.006,0.0009,0.0001};
	double		MutTypeStats[MUTATION_TYPES]	=	{1,0,0,0,0,0,0,0} ;
	double 		MutTypeAcumF[MUTATION_TYPES];
	double		total_adaptMut	= CalculaTotal(MUTATION_TYPES,MutTypeStats);
	MutTypeAcumF[0]		= MutTypeStats[0] / total_adaptMut;
	for (int i=1; i<MUTATION_TYPES; i++){
		MutTypeAcumF[i]= (MutTypeStats[i] / total_adaptMut) + MutTypeAcumF[i-1];
	}
	
	TotalChars	=	contChars(DATA);				//CANTIDAD DE CARACTERES DE LA REFERENCIA
	Reference	=	(char*) malloc(TotalChars*sizeof(char));
	getReference(DATA,Reference);					//SE OBTIENE LA REFERENCIA
	TotalReads	=	B*C;							//NUMERO TOTAL DE READS A GENERAR

	srand(time(NULL));								//SEMILLA DE LOS ALEATOREOS

	for(int ReadsCicle = 0;	ReadsCicle<TotalReads; ReadsCicle++){

		id	=	ReadsCicle+1;
		fprintf(ALIGN,"Read ID:	%"PRIu32"\n",id);

		//B.GENERAR UNA POSICIÓN ALEATORIA EN EL RANGO [0,LengthRef-LengthRead]
		Pos		=	(rand() %((TotalChars-L) - 0 + 1)) + 0;
		fprintf(ALIGN,"Offset de la referencia: %"PRIu32"\n",Pos);
		Read	=	(char*) malloc((L+READ_BIAS)*sizeof(char));
		memcpy(Read,Reference+Pos,L+READ_BIAS);		//SE OBTIENE EL READ DE LA REFERENCIA

		//B.CALCULAR EL MATCHING
		//ORDEN DE LOS ARREGLOS 
		//FORWARD(F), REVERSE(R), COMPLEMENT(C), REVERSE COMPLEMENT(E)
		dado 	= 	LanzarDado();											//LANZAR DADO
		int MatTypeSel  =   BusqBin_Rul(MatTypeAcumF,MATCHING_TYPES,dado); 	//GIRAR LA RULETA
		selMatching(MatTypeSel,L,Read,MT);									//APLICAR EL MATCHING
		fprintf(ALIGN,"Referencia: %s\n",Read);

		//C.CALCULAR LA CANTIDAD DE ERROES
		dado	=	LanzarDado();
		lendesc	=	BusqBin_Rul(ErrorStat,t,dado);
		fprintf(ALIGN,"Cantidad de errores:	%"PRIu16"\n",lendesc);
		if(lendesc	!=	0){
			strand	=	(char) tolower(*MT);
			fprintf(ALIGN,"Matching type %c\n",strand);
			int OffsetAcum	=	0;
			int OffsetAux	=	0;
			for(int i	=	0;	i<lendesc;	i++){	

				Offsets		=   (uint16_t*) malloc(lendesc*sizeof(uint16_t));
   				Oper        =   (uint8_t*)  malloc(lendesc*sizeof(uint8_t));
				BaseRef		=	(uint8_t*)  malloc(lendesc*sizeof(uint8_t));
				BaseRead	=	(uint8_t*)  malloc(lendesc*sizeof(uint8_t));
				Cnts	=	0;	CntS	=	0;	Cntd	=	0;	CntD	=	0;
				Cnti	=	0;	CntI	=	0;	CntT	=	0;	CntC	=	0;

				
				
				//Determinar el tipo de error
				dado	=	LanzarDado();
				int ErrorSel	=	BusqBin_Rul(MutTypeAcumF,MUTATION_TYPES,dado);
				Oper[i]	=	selMutation(ErrorSel);
				fprintf(ALIGN,"Operación de mutación: %c\n",Oper[i]);

				//Aplicar la mutación
				if ((strand=='f')||(strand=='c')) { 
					//Calcular el offset del error para forward
					OffsetAux	=	(rand() %((L-lendesc+i) - (OffsetAcum+1) + 1)) + (OffsetAcum+1);
					Offsets[i]	=	OffsetAux	-	OffsetAcum;
					fprintf(ALIGN,"Offset[%d]	=	%d\n",i,Offsets[i]);
					OffsetAcum	=	OffsetAux;
					FordwardMutation(Oper[i],Read,OffsetAux,BaseRef,BaseRead,id);
					fprintf(ALIGN,"SECUENCE: %s\n",Read);
					fprintf(ALIGN,"Base referece %c",BaseRef[id-1]);
					fprintf(ALIGN,"Base read: %c",BaseRead[id-1]);
				}else{
					//RevereseMutation();
				}

				//Salida
			}
		}else{
			//EN ESTE CASO NO HAY ERRORES
			strand	=	*MT;
			fprintf(ALIGN,"Match Perfecto del tipo %c\n",strand);
			fprintf(FASTQ,"@%s %d\n",I,id);
			char	*FinalRead	=	(char*) malloc(L*sizeof(char));
			memcpy(FinalRead,Read,L);
			fprintf(FASTQ,"%s\n",FinalRead);
			fprintf(FASTQ,"+\n");
			for(int i	=	0;	i<L;	i++){
				fprintf(FASTQ,"%s",Q);
			}
			fprintf(FASTQ,"\n\n");
			fprintf(FASTQSEQ,"Read: %d, Secuencia:	%s\n",id,FinalRead);
		}

		fprintf(ALIGN,"\n");
		free(Read);

	}
	

	//LIBERAR ESPACIO DE MEMORIA
	free(RefName);
	free(RefFastq);
	free(RefFastqseq);
	free(RefAlign);
	free(RefMeta);
	free(Reference);

	//CERRAR LOS ARCHIVOS
	fclose(FASTQ);
	fclose(FASTQSEQ);
	fclose(ALIGN);
	fclose(META);
	return 0;
}

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


void selMatching(int MatTypeSel, int L, char *READ, char *MT){
	switch(MatTypeSel){
		case 0:                     //FORWARD MATCH
			memcpy(MT,"F",1);
		break;
		case 1:                     //REVERSE MATCH
			memcpy(MT,"R",1);
			ReverseRead(READ,L);
		break;
		case 2:                     //COMPLEMENT MATCH
			memcpy(MT,"C",1);
			ComplementRead(READ,L);
		break;
		case 3:                     //REVERSE COMPLEMENT MATCH
			memcpy(MT,"E",1);
			ComplementRead(READ,L);
			ReverseRead(READ,L);              
		break;
		default:    printf ("**Error in the matching selection, wrong input base**");
	}
}



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
	}
}

char selBase(int sel,char Base){
	switch(Base){
		case 'A':
			switch(sel){
				case 1: return 'C'; break;
				case 2: return 'G'; break;
				case 3: return 'T'; break;
				case 4: return 'N'; break;
			}
		break;
		case 'C':
			switch(sel){
				case 1: return 'A'; break;
				case 2: return 'G'; break;
				case 3: return 'T'; break;
				case 4: return 'N'; break;
			}
		break;
		case 'G':
			switch(sel){
				case 1: return 'A'; break;
				case 2: return 'C'; break;
				case 3: return 'T'; break;
				case 4: return 'N'; break;
			}
		break;
		case 'T':
			switch(sel){
				case 1: return 'A'; break;
				case 2: return 'C'; break;
				case 3: return 'G'; break;
				case 4: return 'N'; break;
			}
		break;
	}
}

void FordwardMutation(uint8_t Oper,char *Read,uint16_t Offset, uint8_t *BaseRef,uint8_t *BaseRead,uint32_t id){

	//VECTOR DE PROBABILIDADES DE LAS BASES PARA LAS SUBSTITUCIONES
	double	BasesStats[4]	=	{0.3,0.3,0.3,0.1};
	double 	BasesAcum[4];
	double	total_adaptBases	= CalculaTotal(4,BasesStats);

	switch((char)Oper){
		
		case 's':
			
			BasesAcum[0]		= 	BasesStats[0] / total_adaptBases;
			for (int i=1; i<4; i++){
				BasesAcum[i]= (BasesStats[i] / total_adaptBases) + BasesAcum[i-1];
			}
			BaseRef[id-1]	=	Read[Offset];
			//Calcular la base destino
			double dado		=	LanzarDado();
			int sel	=	BusqBin_Rul(BasesAcum,4,dado);
			char BaseDest			=	selBase(sel,Read[Offset]);
			BaseRead[id-1]	=	BaseDest;
			Read[Offset]	=	BaseDest;
		break;
	}

}

