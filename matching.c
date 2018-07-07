/*
 *  @Developer:     Juan Camilo Peña Vahos
 *  @Last Revised:  27/06/2018
 *  @Description:   Algoritmo que determina el tipo de matching del read, de acuerdo con
 *                  una probabilidad.
*/	


//LIBRERÍAS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#define MATCHING_TYPES 4

//DECLARACIÓN DE FUNCIONES
void ReverseRead(char*,long);
void ComplementRead(char*,long);
int  BusqBin_Rul(double[], int, double);
double CalculaTotal(int,double[]);
double LanzarDado();

int main (int argc, char *argv[]) {	
    
	//ARGUMENTOS DE ENTRADA
	int L	=   30;

	//VARIABLES DEL PROCESO PARA SALIDA								
	char *MT;										//MATCHING TYPE

    //PARA EL PROBAR EL MATCHING, SE ASUME LA SIGUIENTE CADENA.
    int MatType;                                   
    char Read[L];
    //ORDEN DE LOS ARREGLOS 
    //FORWARD(F), REVERSE(R), COMPLEMENT(C), REVERSE COMPLEMENT(E)
    double MatTypeStats[MATCHING_TYPES]   =   {0.4,0.4,0.1,0.1};  
    double Acum_fun[MATCHING_TYPES];

    strcpy(Read,"GGGCGGCGACCTCGCGGGTTTTCGCTATTTA");
    printf("READ ORIGINAL: %s\n",Read);

    //CALCULAR EL ARREGLO CON LA FUNCIÓN DE PROBABILIDADES ACUMULADAS
	double total_adapt = CalculaTotal(MATCHING_TYPES,MatTypeStats);
	Acum_fun[0]= MatTypeStats[0] / total_adapt;
	for (int i=1; i<MATCHING_TYPES; i++){
        Acum_fun[i]= (MatTypeStats[i] / total_adapt) + Acum_fun[i-1];
    }
    printf("PROBABILIDADES POR INDIVIDUO: ");
    for(int i = 0; i<MATCHING_TYPES; i++){
        printf("%lf, ",MatTypeStats[i]);
    }
    printf("\n");
    printf("PROBABILIDADES ACUMULADAS: ");
    for(int i = 0; i<MATCHING_TYPES; i++){
        printf("%lf, ",Acum_fun[i]);
    }
    printf("\n");

    //LANZAR EL DADO
    double dado = LanzarDado();
    printf("Resultado del dado = %lf\n", dado);

    //GIRAR LA RULETA
    int MatTypeSel  =   BusqBin_Rul(Acum_fun,MATCHING_TYPES,dado);
    printf("Selección de la ruleta = %d\n", MatTypeSel);

    switch(MatTypeSel){
        case 1:                     //FORWARD MATCH
            MT  =   "F";
            printf("F: %s\n",Read);
        break;
        case 2:                     //REVERSE MATCH
            MT  =   "R";
            ReverseRead(Read,L);
            printf("R: %s\n",Read);
        break;
        case 3:                     //COMPLEMENT MATCH
            MT  =   "C";
            ComplementRead(Read,L);
            printf("C: %s\n",Read);
        break;
        case 4:                     //REVERSE COMPLEMENT MATCH
            MT  =   "E";
            ComplementRead(Read,L);
            ReverseRead(Read,L);
            printf("E: %s\n",Read);                
        break;
        default:    printf ("**Error in the matching selection, wrong input base**");
    }
    
}

//Implementa el Inversor de reads 
void ReverseRead(char *Read, long length){
	char aux;
    int resto = length%2;
	for (int i=0; i<length/2;i++){
        //printf("READI = %c , READLEN = %c, i =  %d, resto = %d\n",Read[i],Read[length-i-1], i, length-i-1);
		aux=Read[length-i-1];
		Read[length-i-1]=Read[i];
		Read[i]=aux;
	}
}

//Implementa el complementador de reads (T,A)(C,G)
void ComplementRead(char *Read, long length){
	char Compl;
	int i;
	for (i=0; i<length;i++){
		switch(Read[i]){
		    case 'A': Compl='T'; break;
		    case 'a': Compl='T'; break;
		    case 'C': Compl='G'; break;
		    case 'c': Compl='G'; break;
		    case 'G': Compl='C'; break;
		    case 'g': Compl='C'; break;
		    case 'T': Compl='A'; break;
		    case 't': Compl='A'; break;
		    case 'N': Compl='N'; break; 
		    case 'n': Compl='N' ;break;
		    default: printf ("**Error building the complement, wrong input base**");
		}

		Read[i]=Compl;
	}

}

//BusqBin_Rul: función que utiliza el procedimiento de la Busq. Binaria
//             para ubicar el elemento seleccionado a través del mecanismo
//	            de la ruleta a fin de seleccionar los Padres.
//@param: double prob[] : vector de probabilidades acumuladas
//@param: n :   número de elementos en prob[]
//@param: x :   resultado de lanzar los dados (es decir, usar rand para generar un número)
int BusqBin_Rul(double prob[], int n, double x){

	int primero,ultimo,central;
	short encontrado;
    primero     =   1;
    ultimo      =   n;
    encontrado  =   0;
	while ((primero <= ultimo) && (encontrado==0)){
   	    central = (primero + ultimo)/2;
        if      (x  ==   prob[central]) encontrado = 1;
        else if (x  >    prob[central]) primero = central + 1;
        else                            ultimo = central - 1;
   }
   return(central);
}

//CalculaTotal: Calculo de la Sumatoria de las Func. Aadaptación
//              de cada Individuo de la Poblacion.
//@param: a[] : vector de probabilidades de cada individuo, sin acumular
//@param: n :   número de elementos en a[]
double CalculaTotal(int n, double a[]){

	int i;
	double aux=0;
    for (i=0;i<n;i++) {	aux = aux + a[i];}
	return aux;
}

//LanzarDado: devuelve un valor entre 0 y 1 aleatorio
double LanzarDado(){
    srand(time(0));
    double dado = 0+(1-0)*rand()/((double)RAND_MAX);
    return dado;
}



