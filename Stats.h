/*
 ============================================================================
 Name        : 	Stats.h
 Author      : 	Juan Camilo Peña Vahos - Aníbal Guerra - Sebastian Isaza Ramírez
 Version     :
 Copyright   : 	This project is totally opensource
 Description :	Archivo de cabecera para las funciones que dan soporte a la parte estadística
 ============================================================================
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
