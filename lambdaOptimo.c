/**
 * @Developer:      Juan Camilo Peña Vahos
 * @Description:    Imprime el rango de lambdas óptimo de acuerdo a la longitud del Read
 *                  para una función de probabilidad del tipo exponencial
 * @LastRevised:    20/06/2018
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define E           2.71828
#define LAMBDA      0.0025
#define TINICIAL    0

double exp_1(int,double);
int AdjustKExp(int,double,int);

int main (int argc, char *argv[]) {	
    
    //Lets assume that our reads have a minimun length of 10 and a Max length of 1024
    //And lets assume that each Read have a maximun number of errors about the 25% of its length  
    int ReadLengthMin       =   10;
    int ReadLengthMax       =   1024;
    double ErrPercentage    =   0.25;
    double lambda           =   0.0;
    FILE *pf;
    pf  =   fopen("Lambda óptimo","w");

    //However, if you want to change those params
	if(argc>1){
		for (int i = 1; i < argc; i++)	{
			if (strcmp(argv[i], "-MIN") == 0)	ReadLengthMin	=	atoi(argv[i+1]);
			if (strcmp(argv[i], "-MAX") == 0)	ReadLengthMax	=	atoi(argv[i+1]);
            if (strcmp(argv[i], "-PER") == 0)	ErrPercentage	=	atof(argv[i+1]);
		}
	}

    //Lets calculate the minimun and the maximun number of errors per read
    int MinK    =   (int) floor(ErrPercentage   *   ReadLengthMin);
    int MaxK    =   (int) floor(ErrPercentage   *   ReadLengthMax);
    fprintf(pf,"MinK    =   %d\n",MinK);
    fprintf(pf,"MaxK    =   %d\n",MaxK); 

    //For every possible number of errors
    for(int i =  MinK;  i <= MaxK;  i++){
        int flag    =   0;
        while(flag  ==  0){
            lambda  =   lambda  +   LAMBDA;
            double    aux     =   0.0;
            int       cont    =   TINICIAL;
            while ( (aux<=0.9998)   &&  (cont<=i) ){
                aux =   aux +   exp_1(cont,lambda);
                cont++;
            }   
            if( (aux>0.9998) &&  (aux<=1) ){
                flag    =   1;
                fprintf(pf,"Para el valor de K  =   %d, se tiene un lambda de %f ",i,lambda);
                fprintf(pf,"y una probabilidad de %f\n",aux);
                aux =   0;
                lambda  =   0;
            }  
            if(aux>1){
                flag    =   1;
                aux =   aux -   exp_1(cont,lambda);
                lambda  =   lambda  -   LAMBDA;   
                fprintf(pf,"Para el valor de K  =   %d, se tiene un lambda de %f ",i,lambda);
                fprintf(pf,"y una probabilidad de %f\n",aux);
                aux = 0;
                lambda  =  0;
            }
        }
    }
    fclose(pf);
    
}

double exp_1(int t , double lambda){
    return (lambda*pow(E, (-1.0)*lambda*((double)t)));
}
