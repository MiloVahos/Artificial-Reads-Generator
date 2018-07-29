/*
 *  @Developers:    Juan Camilo Peña Vahos	- Aníbal Guerra	- Sebastian Isaza Ramirez
 *  @Last Revised:  28/07/2018
 *  @Description:   Módulo para las funciones estadísticas
*/ 

#ifndef STATS_H
#define STATS_H

#define EXP 			2.71828	//VALOR DEL NÚMERO e
#define TINICIAL		0		//VALOR INICIAL DE t


int  		BusqBin_Rul         (double[],int,double);
double 		CalculaTotal        (int,double[]);
double 		LanzarDado          ();
int 		AdjustKExp          (int,double,int);
void		generarRuletaExp    (double*,int,double);
double 		exp1                (int,double);


#endif