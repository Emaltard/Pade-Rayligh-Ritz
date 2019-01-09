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

void print_matrix(Matrix* m){
	for(int i = 0; i < m->size[0]; i++){
		for(int j = 0; j < m->size[1]; j++){
			printf("%0.3f", m->data[i][j]);
		}
		printf("\n");
	}
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

/* Produit matrice-vector : A*b = c*/
void prod_mat_vect(Matrix* a, Vector* b, Vector* c){
    int i, j;
    
    //omp_set_num_threads(4);
    #pragma omp parallel for private(i, j)
    for ( i = 0; i < b->size; i++){
        c->data[i] = 0;
        for (j = 0; j < b->size; j++){
            c->data[i] = c->data[i] + a->data[i][j]*b->data[j];
        }
    }
}

/* Produit scalaire */
double prodScalaire(Vector* v1, Vector* v2){
    double res = 0;
    int i;

    //omp_set_num_threads(thr);
    #pragma omp parallel for private(i) reduction(+:res)
    for (i = 0; i < v1->size; i++){
        res += v1->data[i]*v2->data[i];
    }
    return res;
}

/* Etape 4 de l'algorithme */
// C doit etre de taille 2m [0,...., 2m-1], C[0] = C0
void step4(Vector* C, Matrix* A, Vector* y, int m){

	Vector* y1;
	y1 = init_vector(y->size);
	prod_mat_vect(A, y, y1);
	for(int i = 1; i <= m-1; i++){
		C->data[2*i-1] = prodScalaire(y1, y);
		C->data[2*i] = prodScalaire(y1, y1);
		y = y1;
		prod_mat_vect(A, y, y1);
	}
	C->data[2*m-1] = prodScalaire(y1, y);
}

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
	PRR(m, v);
	free_vector(v);
	return 0;	
}
