/*
 ============================================================================
 Name        : 	ArtificialReadsGenerator.c
 Author      : 	Juan Camilo Peña Vahos - Aníbal Guerra - Sebastian Isaza Ramírez
 Version     :
 Copyright   : 	This project is totally opensource
 Description :	Conjunto de funciones que dan soporte a la parte estadística
 ============================================================================
 */

#include <stdlib.h>
#include <math.h>

#include "Stats.h"

//exp1	:   Calcula el valor de la exponencial para un t en particular
//@param:   t       :   Valor de t actual
//@param:   lambda  :   Valor ajustable, es la constante de la exponencial
double exp1(int t, double lambda){
    return (lambda*pow(EXP, (-1.0)*lambda*((double)t)));
}

//AdjustKExp:           Iterar buscando que la suma de Pi, 0<=i<=k e aproxime más a 0.99999999
//@param: start		: 	Valor inicial de t -> Suele ser cero (t es el número de elementos de la ruleta)
//@param: lambda	:   Valor ajustable, es la constante de la exponencial
//@param: MaxK 		:   Número máximo de errores de un read
int AdjustKExp(int start, double lambda, int MaxK ){

    double    aux     	=   0.0;
    int       auxt    	=   start;

    while ((aux	<=	0.9998)	&&	(auxt	<=	MaxK)){
        aux	=	aux	+	exp1(auxt,lambda);
        auxt++;
    }
    return (auxt);
}

//generarRuletaExp      :   Crea una ruleta con cada opción definida por la función de una
//                          distribución exponencial
//@param:   ErrorStat   :   Vector donde se va a guardar la ruleta
//@param:   t   		: 	Casilla de la ruleta a llenar
//@param:   lambda	    :   Valor ajustable, es la constante de la exponencial
void generarRuletaExp(double *ErrorStat,int t,double lambda){

	double aux	=	0.0;

	for(int i=0;	i<t;	i++){

		aux	=	aux	+	exp1(i,lambda);
		ErrorStat[i]	=	aux;
	}
}

//BusqBin_Rul: función que utiliza el procedimiento de la Busq. Binaria
//             para ubicar el elemento seleccionado a través del mecanismo
//	            de la ruleta
//@param: double prob[] :   vector de probabilidades acumuladas
//@param: n             :   número de elementos en prob[]
//@param: x             :   resultado de lanzar los dados (es decir, usar rand para generar
//                          un número)
int BusqBin_Rul(double prob[], int n, double x){

	int primero,ultimo,central;
	short encontrado;
    primero     =   0;
    ultimo      =   n-1;
    encontrado  =   0;

	while ((primero <= ultimo)  &&  (encontrado==0)){

   	    central =   (primero + ultimo)/2;
        if          (x  ==   prob[central]) encontrado = 1;
        else if     (x  >    prob[central]) primero = central + 1;
        else                                ultimo = central - 1;
   }
   return(primero);
}

//CalculaTotal: Calculo de la Sumatoria de las Func. Adaptación
//              de cada Individuo de la Poblacion.
//@param: a[] : vector de probabilidades de cada individuo, sin acumular
//@param: n :   número de elementos en a[]
double CalculaTotal(int n, double a[]){

	double aux	=	0;
    for (int i=0;i<n;i++) {	aux = aux + a[i];}
	return aux;
}

//LanzarDado: devuelve un valor entre 0 y 1 aleatorio
double LanzarDado(){
    double dado = 0+(1-0)*rand()/((double)RAND_MAX);
    return dado;
}

