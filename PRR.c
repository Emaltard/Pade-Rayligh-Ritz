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
		printf("%0.3f ", v->data[i]);
	}
	printf("\n");
}

void print_matrix(Matrix* m){
	for(int i = 0; i < m->size[0]; i++){
		for(int j = 0; j < m->size[1]; j++){
			printf("%0.3f ", m->data[i][j]);
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

void fill_matrix_with_random_values(Matrix* m){
	for(int i = 0; i < m->size[0]; i++){
		for(int j = 0; j < m->size[1]; j++){
			m->data[i][j] = frand_a_b(1.0, 1000.0);
		}
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

void inversion_matrix(Matrix* m){
	int info;
	LAPACKE_dgetrf(LAPACK_ROW_MAJOR, m->size[0], m->size[1], m->data[0], m->size[0], &info);
	printf("LAPACKE_dgetrf: %d\n", info);

	LAPACKE_dgetri(LAPACK_ROW_MAJOR, m->size[0], m->data[0], m->size[0], &info);
	printf("LAPACKE_dgetri: %d\n", info);
}

/* Etape 4 de l'algorithme */
// C doit etre de taille 2m [0,...., 2m-1], C[0] = C0
void step4(Vector* C, Matrix* A, Vector* y0, int m, Vector* V[m]){
	V[0] = y0;
	printf("C:\n");
	print_vector(C);
	printf("V0:\n");
	print_vector(V[0]);
	printf("y0:\n");
	print_vector(y0);
	printf("A:\n");
	print_matrix(A);
	Vector* y1;
	y1 = init_vector(y0->size);
	prod_mat_vect(A, y0, y1);
	printf("y1:\n");
	print_vector(y1);
	for(int i = 1; i <= m-1; i++){
		C->data[2*i-1] = prodScalaire(y1, y0);
		C->data[2*i] = prodScalaire(y1, y1);
		y0 = y1;
		V[i] = y0;
		prod_mat_vect(A, y0, y1);
	}
	C->data[2*m-1] = prodScalaire(y1, y0);
}

void fill_B_and_C(Matrix* B, Matrix* C, Vector* V){
	int s = V->size/2;
	for(int i = 0; i < s; i++){
		for(int j = 0; j < s; j++){
			B->data[i][j] = V->data[i+j];
			C->data[i][j] = V->data[i+j+1];
		}
	}
}

/* Fonction Algorithme itérative PRR */
void PRR(int m, Vector* x, Matrix* A){

	// Normalisation de x + calcul de y0
	double norm = vect_norm(x);
	Vector* y = normalize(x, norm);
	
	// C0 = || y0 ||^2
	double C1;
	norm = vect_norm(y);
	C1 = norm*norm;

	// Calcul de C1, C2,....C2m-1
	Vector* C;
	C = init_vector(2*m);
	C->data[0] = C1;

	Vector* V[m];
	for(int i = 0; i<m; i++){
		V[i] = init_vector(y->size);
	}
	step4(C, A, y, m, V);

	// Calcul de B^ et C^.
	// Matrix* B, *Cc;
	// B = init_matrix(m,m);
	// Cc = init_matrix(m,m);

	// fill_B_and_C(B, Cc, C);

	// print_vector(C);
	// print_matrix(Cc);
	// inversion_matrix(B);

	// Calcul de Xm


}


int main(int argc, char** argv){
	
	int m = 10;	// Matrix* B, *Cc;
	// B = init_matrix(m,m);
	// Cc = init_matrix(m,m);

	// fill_B_and_C(B, Cc, C);

	// print_vector(C);
	Vector* v;
	Matrix* A;
	v = init_vector(m);
	A = init_matrix(m, m);

	srand(time(0));

	fill_matrix_with_random_values(A);
	fill_vector_with_random_values(v);

	PRR(m, v, A);

	free_vector(v);
	free_matrix(A);
	return 0;	
}
