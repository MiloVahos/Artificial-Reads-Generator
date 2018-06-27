#include <stdio.h>
#include <math.h>
#define E 2.71828

int main(void) {
	//
	// int Max=300;double Lambda=0.0225; //[0.02 (max231Erroers) a 0.4 (max 5 errores)]--> Quiero al menos 20 muestras entonces. 
		//Con un espaciado de : 0.0025 (aprobado 28/05/18)  
		//Esto es para el calculo de la probabilidad de errores
	int  start =0; 
	
	//Ahora para el calculo de los offsets [0,1023]
	double Lambda=0.0057394; 	int Max=1023;
	
	printf ("For lambda : %g , K es : %i ", Lambda, AdjustKExp(start, Lambda, Max )  );
	return 0;
}

double exp_1(int t , double lambda){
    return (lambda*pow(E, (-1.0)*lambda*((double)t)));
}
//reliawiki.org/index.php/The_Exponential_Distribution
int AdjustKExp(int t, double P0, int MaxK ){
      double aux=0.0;
      int  FinalK=0.0, aux_t=t;

       //Iterar buscando que la suma de Pi, 0<=i<=k  se aproxime más a 0.99999999 y menor o igual que 1.0
       while ((aux<=0.999999999)&&(aux_t<=MaxK)){
           aux=aux+exp_1(aux_t , P0);
           aux_t++;
           //printf (" Acumulado en la iter %d es : %f  \n", aux_t, aux);
       }      
     return (aux_t);
}

/*//ShellSort: Ordena los Arreglos de Func. de Apatación
//           y Poblacion en base a la Func. de Adaptación
void ShellSort (Poblacion p, int n, double a[]){
	int intervalo, d, j, k;
	Cromosoma crom;
	double aux_fcrom;

	intervalo = n/2;
	while (intervalo > 0 ){
		for (d = ( intervalo ); d<n; d++){
			j = d - intervalo;
			while (j >= 0){
				k = j + intervalo;
				if  (a[j] >= a[k])  j = -1;
				else {
					aux_fcrom= a[j];
					a[j] = a[k];
					a[k] = aux_fcrom;
					crom= p[j];
					p[j] = p[k];
					p[k] = crom;
				}//fin si->sino
				j = j - intervalo;
			}//fin mientras
		}//fin para d
		intervalo = intervalo / 2;
	}//fin mientras
}//fin Ordenamiento ShellSort
*/
//CalculaTotal: Calculo de la Sumatoria de las Func. Aadaptación
//              de cada Individuo de la Poblacion.
double CalculaTotal(int n, double a[]){
	int i;
	double aux=0;

   for (i=0;i<n;i++) {	aux = aux + a[i];}
	return aux;
}//Fin Calculo Total de la Sumatoria de la FADAPT

//Proba_Ind: Genera un Arreglo con la Probabilidad de cada Individuo
//           de ser escogido por la ruleta.
double Proba_Ind(double prob[], int n, double a[]){
	int i;
	double total_adapt = CalculaTotal(n,a);

	prob[0]= a[0] / total_adapt;
	for (i=1;i<n;i++)
	prob[i]= (a[i] / total_adapt) + prob[i-1];
	return(total_adapt/n) ;
}//Fin Calcular Probabilidad de cada Individuo de ser escogido

//BusqBin_Rul: función que utiliza el procedimiento de la Busq. Binaria
//             para ubicar el elemento seleccionado a través del mecanismo
//	            de la ruleta a fin de seleccionar los Padres.
int BusqBin_Rul(double prob[], int n, double x){
	int primero,ultimo,central;
	short encontrado;

	primero= 1 ; ultimo= n ; encontrado=0;
	while ((primero <= ultimo) && (encontrado==0)){
   	central = (primero + ultimo)/2;
      if (x == prob[central]) encontrado = 1;
      else
      	if (x > prob[central]) primero = central + 1;
         else ultimo = central - 1;
   }//fin mientras
	return(central);
}//fin Busqueda Binaria segun la Ruleta

//Ruleta: escoge un individuo de la poblacion.
//Modificar para que si la ruleta sale fuera del tope se tome com prob 0
/*unsigned int Ruleta(double prob[], int n){
	double aux=(float)rand()/(float)RAND_MAX;

	if (aux<prob[0]) return(0);
     else return(BusqBin_Rul(prob,n,aux));
}//fin Ruleta
*/
//SeleccPadre: Selecciona un cromosoma de la poblacion
//             para utilizarlo en la reproduccion.
/*Cromosoma SeleccPadre(Poblacion P, int pos){
	return (P[pos]);
}*///fin Seleccionar Padre