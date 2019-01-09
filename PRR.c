// The Padé-Rayligh-Ritz method
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
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
	printf("\n");
}

/* Calcul de la norme d'un vecteur */
double vect_norm(double* x, int size){
	double result = 0;
	for(int i = 0; i < size; i ++){
		result = x[i]*x[i];
	}

	return sqrt(result);
}

/* Normalisation d'un vecteur */
void normalize(double*x, double norm, int size){
	for(int i = 0; i < size; i++){
		x[i] = x[i]/norm;
	}
}

/* Fonction Algorithme itérative PRR */
void PRR(int m, double* x){

	// Normalisation de x + calcul de y0
	double* y = x;
	double norm = vect_norm(y, m);
	normalize(y, norm, m);
	
	// C0 = || y0 ||^2
	double C1, C2;
	norm = vect_norm(y, m);
	C1 = norm*norm;
	if ((int)C1 == 1) {
		printf("OK\n");
	}



}


int main(int argc, char** argv){
	
	int m = 10;
	double* v = malloc(sizeof(double)*m);
	srand(time(0));

	init_vector(v, m);
	PRR(m, v);
	return 0;	
}
