// The Padé-Rayligh-Ritz method
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "lapacke.h"

struct vector
{
	int size;
	double *data;
};
typedef struct vector Vector;

struct matrix
{
	int size[2];
	double **data;
};
typedef struct matrix Matrix;

Vector* init_vector(int size){
	Vector* v = malloc(sizeof(Vector));
	v->size = size;
	v->data = (double*)calloc(size, sizeof(double));
	return v;
}

Matrix* init_matrix(int size_x, int size_y){
	Matrix* m = malloc(sizeof(Matrix));
	m->size[0] = size_x;
	m->size[1] = size_y;
	m->data = (double**)calloc(size_x, sizeof(double*));
	m->data[0] = (double *)calloc(size_x * size_y, sizeof(double));
	for (int i = 0; i < size_x; i++)
	{
		m->data[i] = m->data[0] + i * size_y;
	}
	return m;
}

void free_vector(Vector* v){
	free(v->data);
	free(v);
}

void free_matrix(Matrix* m)
{
	free(m->data);
	free(m);
}

void print_vector(Vector* v){
	for(int i = 0; i < v->size; i++){
		printf("%0.3f", v->data[i]);
	}
	printf("\n");
}

double frand_a_b(double a, double b){
    return ( rand()/(double)RAND_MAX ) * (b-a) + a;
}

/* Init vector randomly */
void fill_vector_with_random_values(Vector* v){
	for(int i = 0; i < v->size; i++){
		v->data[i] = frand_a_b(1.0, 1000.0);
	}
}

/* Calcul de la norme d'un vecteur */
double vect_norm(Vector* x){
	double result = 0;
	for(int i = 0; i < x->size; i ++){
		result = x->data[i]*x->data[i];
	}

	return sqrt(result);
}

/* Normalisation d'un vecteur */
Vector* normalize(Vector* x, double norm){
	Vector* v;
	v = init_vector(x->size);
	for(int i = 0; i < x->size; i++){
		v->data[i] = x->data[i]/norm;
	}
	return v;
}

// void compute_C_and_y(double* C, ){

// }

/* Fonction Algorithme itérative PRR */
void PRR(int m, Vector* x){

	// Normalisation de x + calcul de y0
	double norm = vect_norm(x);
	Vector* y = normalize(x, norm);
	
	// C0 = || y0 ||^2
	double C1, C2;
	norm = vect_norm(y);
	C1 = norm*norm;
	printf("C : %f\n", C1);



}


int main(int argc, char** argv){
	
	int m = 10;
	Vector* v;
	v = init_vector(m);
	srand(time(0));

	fill_vector_with_random_values(v);

	free_vector(v);
	PRR(m, v);
	return 0;	
}
