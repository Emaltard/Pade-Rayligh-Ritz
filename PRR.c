// The Padé-Rayligh-Ritz method
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "lapacke.h"

#define K 5		// Nombre de valeurs propres
#define P 0.01	// Precision

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
	// On alloue la mémoire pour le tableau data. Ce tableau 2D sera continu en mémoire
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
		v->data[i] = frand_a_b(1.0, 10.0);
	}
}

void fill_matrix_with_random_values_symetric(Matrix* m){
	for(int i = 0; i < m->size[0]; i++){
		for(int j = i; j < m->size[1]; j++){
			double value = frand_a_b(1.0, 10.0);
			m->data[i][j] = value;
			m->data[j][i] = value;
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

void prod_mat_mat(Matrix* A, Matrix* B, Matrix*C){
    int i, j, k;

    //omp_set_num_threads(thr);
    #pragma omp parallel shared(A,B,C) private(i,j,k) 
    {
        #pragma omp for schedule(static)
        for (i = 0; i < A->size[0]; i++){
            for (j = 0; j < B->size[1]; j++){
                C->data[i][j] = 0;
                for (k = 0; k < B->size[1]; k++)
                C->data[i][j] += A->data[i][k]*B->data[k][j];
            }
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

Vector* scalar_mult_vector(double scalar, Vector* v1){
	Vector* v2 = init_vector(v1->size);
	for(int i = 0; i < v1->size; i++){
		v2->data[i] = scalar * v1->data[i];
	}
	return v2;
}

Vector* vector_minus_vector(Vector* v1, Vector* v2){
	Vector* v3 = init_vector(v1->size);
	for(int i = 0; i < v1->size; i++){
		v3->data[i] = v1->data[i] - v2->data[i];
	}
	return v3;
}

double max_in_vector(Vector* v1){
	double max = 0.0;
	for(int i = 0; i < v1->size; i++){
		if(v1->data[i] > max){
			max = v1->data[i];
		}
	}
	return max;
}

void inversion_matrix(Matrix* m){
	int info;
	Matrix* res = init_matrix(m->size[0], m->size[1]);

	LAPACKE_dgetrf(LAPACK_ROW_MAJOR, m->size[0], m->size[1], m->data[0], m->size[0], &info);
	printf("LAPACKE_dgetrf: %d\n", info);

	LAPACKE_dgetri(LAPACK_ROW_MAJOR, m->size[0], m->data[0], m->size[0], &info);
	printf("LAPACKE_dgetri: %d\n", info);
}

/* Etape 4 de l'algorithme */
// C doit etre de taille 2m [0,...., 2m-1], C[0] = C0
void step4(Vector* C, Matrix* A, Vector* y0, int m, Vector* V[m]){
	V[0] = y0;
	Vector* y1;
	y1 = init_vector(y0->size);
	prod_mat_vect(A, y0, y1);

	for(int i = 1; i < m-1; i++){
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

void step5(int m, Matrix* Xm, Matrix* Vm, Vector* y[m], Vector* val_ritz, Vector* vect_ritz[m]){

	int i_max = 0;
	for(int i = 0; i < m; i++){
		// Calcul de vecteurs de Ritz
		prod_mat_vect(Vm, y[i], vect_ritz[i]); 
		// Calcul de valeurs de Ritz
		val_ritz->data[i] = 1; 	
	}
}

Vector*  step6(Matrix* A, int i, Vector* ritz_vectors[i], Vector* val_ritz){
	Vector* residus = init_vector(i);
	for(int j = 0; j < i ; j++){
		Vector* v1 = init_vector(ritz_vectors[j]->size);
		prod_mat_vect(A, ritz_vectors[j], v1);
		Vector* v2 = scalar_mult_vector(val_ritz->data[j], ritz_vectors[j]);

		
		residus->data[j] = vect_norm(vector_minus_vector(v1, v2)); 
	}
	return residus;
}

Matrix* convert_vector_array_to_matrix(int m, Vector* y[m]){
	Matrix* res = init_matrix(y[0]->size, m);
	for(int i = 0; i < res->size[0]; i++){
		for(int j = 0; j < m; j++){
			res->data[i][j] = y[j]->data[i];
		}
	}
	return res;
}

/* Fonction Algorithme itérative PRR */
void PRR(int m, Vector* x, Matrix* A){
	int N = A->size[0];
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
		V[i] = init_vector(N);
	}
	step4(C, A, y, m, V);

	// Calcul de B^ et C^.
	Matrix* B, *Cc;
	B = init_matrix(m,m);
	Cc = init_matrix(m,m);
	fill_B_and_C(B, Cc, C);

	// Calcul de Xm
	inversion_matrix(B);
	Matrix* Xm = init_matrix(m, m);
	prod_mat_mat(B, Cc, Xm);
	print_matrix(B);

	// Calcul des valeurs propres et vecteurs propres de Xm
	Matrix* Vm = convert_vector_array_to_matrix(m, V);
	Vector* val_ritz = init_vector(m);
	Vector* vect_ritz[m];
	for(int i = 0; i<m; i++){
		vect_ritz[i] = init_vector(N);
	}
	step5(m, Xm, Vm, V, val_ritz, vect_ritz);

	// Test pour la projection ls
	Vector* residus = step6(A, m, vect_ritz, val_ritz);
	double max = max_in_vector(residus); 
	if(max_in_vector(residus) <= P){
		// C EST BON
	}else{
		// ON RESTART
		int i;
		for(i=0; i< residus->size; i++){
			if(residus->data[i] == max){
				break;
			}
		}
		PRR(m, vect_ritz[i], A);
	}
	
	print_vector(val_ritz);
	free_vector(y);
	free_matrix(B);
	free_matrix(Cc);
	free_vector(C);
}	


int main(int argc, char** argv){
	
	int m = 2;
	int N = 100;

	Vector* v;
	Matrix* A;
	v = init_vector(N);
	A = init_matrix(N, N);

	srand(time(0));

	fill_matrix_with_random_values_symetric(A);
	fill_vector_with_random_values(v);
	printf("\n\n");

	PRR(m, v, A);

	free_vector(v);
	free_matrix(A);
	return 0;	
}
