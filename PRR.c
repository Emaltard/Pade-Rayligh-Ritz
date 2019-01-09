// The Padé-Rayligh-Ritz method
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "lapacke.h"

double frand_a_b(double a, double b){
    return ( rand()/(double)RAND_MAX ) * (b-a) + a;
}

/* Init vector randomly */
void init_vector(double* x, int size){
	
	for(int i = 0; i < size; i++){
		x[i] = frand_a_b(1.0, 1000.0);
		printf("%0.3f  ", x[i]);
	}
}

/* Fonction Algorithme itérative PRR */
void PRR(int m, int* x){

}


int main(int argc, char** argv){
	
	int m = 10;
	double* v = malloc(sizeof(double)*m);
	srand(time(0));
	
	init_vector(v, m);

	return 0;	
}
