// The Padé-Rayligh-Ritz method
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
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

double frand_a_b(double a, double b){
    return ( rand()/(double)RAND_MAX ) * (b-a) + a;
}

/* Init vector randomly */
void fill_vector_with_random_values(Vector* v){
	for(int i = 0; i < v->size; i++){
		v->data[i] = frand_a_b(1.0, 1000.0);
		printf("%0.3f", v->data[i]);
	}
	printf("\n");
}

// void compute_C_and_y(double* C, ){

// }

/* Fonction Algorithme itérative PRR */
void PRR(int m, int* x){

}


int main(int argc, char** argv){
	
	int m = 10;
	Vector* v;
	v = init_vector(m);
	srand(time(0));

	fill_vector_with_random_values(v);

	free_vector(v);

	return 0;	
}
