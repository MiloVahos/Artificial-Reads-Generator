#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

void intercambiar (int*, int, int);
void ordenarOffsets (int*, int);

int main (void) {

    int lendesc = 256;
    int L = 1024;
    int *Offsets;
    Offsets =  (int*) malloc (lendesc*sizeof(int));
	srand(time(NULL));								

    for ( int i = 0; i < lendesc; i++ ) {

		int flag = 0;
		int OffsetAux = 0;
		do {
			OffsetAux	=	rand() %(((L-1)-(lendesc - i)) + 1 - 0) + 0;
			if ( i != 0 ) {
				int contador = 0;
				for (int j = 0; j < i; j++ ) {
					if ( OffsetAux	== Offsets[j] ){
						contador++;
					} 
				}
				if( contador == 0){
					flag = 1;
				}
			} else {
				flag = 1;
			}
			
		} while ( flag != 1 );
		Offsets[i]	=	OffsetAux;
	}

    for ( int i = 0; i < lendesc; i++ ) {
		printf("%d ",Offsets[i]);
	}
    printf("\n");
    printf("\n");

    //ORGANIZARLOS
    ordenarOffsets(Offsets,lendesc);

    for ( int i = 0; i < lendesc; i++ ) {
		printf("%d ",Offsets[i]);
	}
    printf("\n");

}

void intercambiar (int *Offsets, int i, int j) {
    int tmp = Offsets[i];
    Offsets[i] = Offsets[j];
    Offsets[j] = tmp;
}
void ordenarOffsets (int *Offsets, int N){

    int i, j, k;
    for (i = 0; i < N - 1; i++) {
        for (k = i, j = i + 1; j < N; j++)
            if (Offsets[j] < Offsets[k]) k = j;
            if (k != i) intercambiar (Offsets, i, k);
    }
}