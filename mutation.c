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
#include <math.h>
#include <inttypes.h>

#define E 2.71828
#define READ_BIAS   50      //BIAS PARA OPERAR EN EL READ
#define ERROR_PER   0.25    //PORCENTAJE DE ERRORES EN EL READ


double exp_1(int,double);
int AdjustKExp(int,double,int);

int main (int argc, char *argv[]) {	
    
    uint32_t  Pos;            //Posición de Matching respecto a la referencia
    char      strand;         //Caractér con el sentido del matching
    uint16_t  *Offsets;       //Arreglo de offsets por cada error 
    uint8_t   *Oper;          //Arreglo con la operación por error   
    uint8_t   *BaseRef;       //Arreglo con la base de la referencia (Read Referencia)
    uint8_t   *BaseRead;      //Arreglo con la base después de la mutación (Read Destino) 
    //Las inserciones no tienen base origen
    //Las deleciones no tienen base destino
    //En ambos casos poner 0
    uint32_t  id;              //Identificador consecutivo para cada uno de los Reads
    uint16_t  lendesc;         //Cantidad de errores total en el Read
    uint16_t  Cnts,CntS,Cnti,CntI,Cntd,CntD,CntT,CntC;    //Cantidad de veces que ocurre una mutación


    //ASUMMING READ LENGTH IS 50
    unsigned int    L   =   50;
	char *Read  =   (char*) malloc(L+READ_BIAS*sizeof(char));
    memcpy(Read,"GGGCGGCGACCTCGCGGGTTTTCGCTATTTATGAAAATTTTCCGGTTTAAGGCGTTTCCGTTCTTCTTCGTCATAACTTAATGTTTTTATTTAAAATACC",L+READ_BIAS);

    //Calcular la cantidad de errores del read
    int MaxK        =   ERROR_PER   *   L;
    double lambda   =   0.0225;
    int start       =   0;
    int t           =   AdjustKExp(start,lambda,MaxK);

    Offsets    =   (uint16_t*)   malloc(t*sizeof(uint16_t));
    Oper       =   (uint8_t*)    malloc(t*sizeof(uint8_t));


    for(int i=0;    i<t;    i++){

    }
}

double exp_1(int t , double lambda){
    return (lambda*pow(E, (-1.0)*lambda*((double)t)));
}
int AdjustKExp(int t, double P0, int MaxK ){
    
    double    aux     =   0.0;
    int       FinalK  =   0.0;
    int       auxt    =   t;

    //Iterar buscando que la suma de Pi, 0<=i<=k  se aproxime más a 0.99999999 y menor o igual que 1.0
    while ((aux<=0.999999999)&&(auxt<=MaxK)){
        aux=aux+exp_1(auxt , P0);
        auxt++;
        printf (" Acumulado en la iter %d es : %f  \n", auxt, aux);
    }      
    return (auxt);
}

//BusqBin_Rul: función que utiliza el procedimiento de la Busq. Binaria
//             para ubicar el elemento seleccionado a través del mecanismo
//	            de la ruleta
//@param: double prob[] : vector de probabilidades acumuladas
//@param: n :   número de elementos en prob[]
//@param: x :   resultado de lanzar los dados (es decir, usar rand para generar un número)
int BusqBin_Rul(double prob[], int n, double x){

	int primero,ultimo,central;
	short encontrado;
    primero     =   0;
    ultimo      =   n-1;
    encontrado  =   0;
	while ((primero <= ultimo) && (encontrado==0)){
   	    central = (primero + ultimo)/2;
        if      (x  ==   prob[central]) encontrado = 1;
        else if (x  >    prob[central]) primero = central + 1;
        else                            ultimo = central - 1;
   }
   return(primero);
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
    double dado = 0+(1-0)*rand()/((double)RAND_MAX);
    return dado;
}
