// The Padé-Rayligh-Ritz method
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <omp.h>
#include <mpi.h>
#include "lapacke.h"
#include "mmio.h"

#define K 5	// Nombre de valeurs propres
#define P 0.1 // Precision
#define thr 4  // Nombre de Threads
#define MAX_ITER 2 

/* DGEEV prototype */
/*
   Description.
   ============

   The routine computes for an n-by-n real nonsymmetric matrix A, the
   eigenvalues and, optionally, the left and/or right eigenvectors. The right
   eigenvector v(j) of A satisfies

   A*v(j)= lambda(j)*v(j)
*/
extern void dgeev(char *jobvl, char *jobvr, int *n, double *a,
				  int *lda, double *wr, double *wi, double *vl, int *ldvl,
				  double *vr, int *ldvr, double *work, int *lwork, int *info);

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

Vector *init_vector(int size)
{
	Vector *v = malloc(sizeof(Vector));
	v->size = size;
	v->data = (double *)calloc(size, sizeof(double));
	return v;
}

Matrix *init_matrix(int size_x, int size_y)
{
	Matrix *m = malloc(sizeof(Matrix));
	m->size[0] = size_x;
	m->size[1] = size_y;
	// On alloue la mémoire pour le tableau data. Ce tableau 2D sera continu en mémoire
	m->data = (double **)calloc(size_x, sizeof(double *));
	m->data[0] = (double *)calloc(size_x * size_y, sizeof(double));
	for (int i = 0; i < size_x; i++)
	{
		m->data[i] = m->data[0] + i * size_y;
	}
	return m;
}

void free_vector(Vector *v)
{
	free(v->data);
	free(v);
}

void free_matrix(Matrix *m)
{
	free(m->data);
	free(m);
}

void print_vector(Vector *v)
{
	for (int i = 0; i < v->size; i++)
	{
		printf("%0.3f ", v->data[i]);
	}
	printf("\n");
}

void print_matrix(Matrix *m)
{
	for (int i = 0; i < m->size[0]; i++)
	{
		for (int j = 0; j < m->size[1]; j++)
		{
			printf("%0.3f ", m->data[i][j]);
		}
		printf("\n");
	}
}

//////////////////////////////////////////////////////////////

double frand_a_b(double a, double b)
{
	return (rand() / (double)RAND_MAX) * (b - a) + a;
}

/* Init vector randomly */
void fill_vector_with_random_values(Vector *v)
{
	for (int i = 0; i < v->size; i++)
	{
		v->data[i] = 1; //frand_a_b(1.0, 10.0);
	}
}

void fill_matrix_with_random_values_symetric(Matrix *m)
{
	for (int i = 0; i < m->size[0]; i++)
	{
		for (int j = i; j < m->size[1]; j++)
		{
			double value = frand_a_b(1.0, 10.0);
			m->data[i][j] = value;
			m->data[j][i] = value;
		}
	}
}

/////////////////////////////////////////////////////////////////

/* Calcul de la norme d'un vecteur */
double vect_norm(Vector *x)
{
	double result = 0;
	int i;

	for(i = 0; i < x->size; i ++){
		result += x->data[i]*x->data[i];
	}

	return sqrt(result);
}

/* Normalisation d'un vecteur */
Vector *normalize(Vector *x, double norm)
{
	int i;
	Vector *v;
	v = init_vector(x->size);

	for(i = 0; i < x->size; i++){
		v->data[i] = x->data[i]/norm;
	}
	return v;
}

////// MPI ///////
Matrix *split_matrix(Matrix *A)
{
	int sub_row, sub_size;
	int comm_size, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	sub_row = A->size[0] / comm_size;
	sub_size = sub_row * A->size[1];

	// The matrix is ​​subdivided on each processor
	Matrix *sub_mat = init_matrix(sub_row, A->size[1]);
	MPI_Scatter(&(A->data[0][0]), sub_size, MPI_DOUBLE, &(sub_mat->data[0][0]), sub_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	return sub_mat;
}

void gather_vector(Vector *A, Vector *sub_a)
{
	int sub_size;
	int comm_size, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	sub_size = A->size / comm_size;

	MPI_Gather(&(sub_a->data[0]), sub_size, MPI_DOUBLE, &(A->data[0]), sub_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

/////// END MPI ///////:

/* Produit matrice-vector : A*b = c*/
void prod_mat_vect(Matrix* a, Vector* b, Vector* c){
    int i, j;
    
    // omp_set_num_threads(thr);
    // #pragma omp parallel for private(i, j)
    for ( i = 0; i < c->size; i++){
        c->data[i] = 0;
        for (j = 0; j < b->size; j++){
            c->data[i] = c->data[i] + a->data[i][j]*b->data[j];
        }
    }
}

void prod_mat_vect_mpi(Matrix* a, Vector* b, Vector* c){
    int i, j, rank, size;

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// Le rang Master envoi aux autre le vecteur B
	if (rank == 0)
	{
		for (i = 1; i < size; i++)
		{
			MPI_Send(&(b->data[0]), b->size, MPI_DOUBLE, i, 100, MPI_COMM_WORLD);
		}
	}
	else
	{
		MPI_Recv(&(b->data[0]), b->size, MPI_DOUBLE, 0, 100, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	Matrix *sub_a = split_matrix(a);
	Vector *sub_c = init_vector(sub_a->size[0]);
	prod_mat_vect(sub_a, b, sub_c);
	gather_vector(c, sub_c);
	free_matrix(sub_a);
}

void prod_mat_mat(Matrix *A, Matrix *B, Matrix *C)
{
	int i, j, k;

    // omp_set_num_threads(thr);
    // #pragma omp parallel shared(A,B,C) private(i,j,k) 
    // {
    //     #pragma omp for schedule(static)
        for (i = 0; i < A->size[0]; i++){
            for (j = 0; j < B->size[1]; j++){
                C->data[i][j] = 0;
                for (k = 0; k < B->size[1]; k++)
                C->data[i][j] += A->data[i][k]*B->data[k][j];
            }
        }
    // }
}

/* Produit scalaire */
double prodScalaire(Vector *v1, Vector *v2)
{
	double res = 0;
	int i;

	// omp_set_num_threads(thr);
	// #pragma omp parallel for private(i) reduction(+:res)
	for (i = 0; i < v1->size; i++)
	{
		res += v1->data[i] * v2->data[i];
	}
	return res;
}

Vector *scalar_mult_vector(double scalar, Vector *v1)
{
	Vector *v2 = init_vector(v1->size);
	int i;
	// omp_set_num_threads(thr);
    // #pragma omp parallel for private(i)
	for(i = 0; i < v1->size; i++){
		v2->data[i] = scalar * v1->data[i];
	}
	return v2;
}

Vector *vector_minus_vector(Vector *v1, Vector *v2)
{
	Vector *v3 = init_vector(v1->size);
	for (int i = 0; i < v1->size; i++)
	{
		v3->data[i] = v1->data[i] - v2->data[i];
	}
	return v3;
}

double max_in_vector(Vector *v1)
{
	double max = 0.0;
	for (int i = 0; i < v1->size; i++)
	{
		if (v1->data[i] > max)
		{
			max = v1->data[i];
		}
	}
	return max;
}

////////////////////////////////////////////////////////////
int inversion_matrix(Matrix *m)
{
<<<<<<< HEAD
	int ipiv[m->size[0] + 1];
	LAPACKE_dgetrf(LAPACK_ROW_MAJOR, m->size[0], m->size[1], m->data[0], m->size[0], ipiv);
	LAPACKE_dgetri(LAPACK_ROW_MAJOR, m->size[0], m->data[0], m->size[0], ipiv);
=======
	int info1, info2;
    int ipiv[m->size[0]+1];
    info1 = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, m->size[0], m->size[1], m->data[0], m->size[0], ipiv);
    info2 = LAPACKE_dgetri(LAPACK_ROW_MAJOR, m->size[0], m->data[0], m->size[0], ipiv);

	return info1*info2;
>>>>>>> f7f8895b460f3b66c54d5d7d58ec2c94f416663b
}

/* Etape 4 de l'algorithme */
// C doit etre de taille 2m [0,...., 2m-1], C[0] = C0
void step4(int N, Vector *C, Matrix *A, Vector *y0, int m, Vector *V[m], int rank)
{
	V[0] = y0;
	Vector *y1;
	y1 = init_vector(N);

	prod_mat_vect_mpi(A, y0, y1);
	printf("RANk : %d, size %d\n", rank, y0->size);

	for (int i = 1; i <= m - 1; i++)
	{
		if (rank == 0)
		{
			C->data[2 * i - 1] = prodScalaire(y1, y0);
			C->data[2 * i] = prodScalaire(y1, y1);
			y0 = y1;
			V[i] = y0;
		}
		prod_mat_vect_mpi(A, y0, y1);
	}
	C->data[2 * m - 1] = prodScalaire(y1, y0);
}

void fill_B_and_C(Matrix *B, Matrix *C, Vector *V)
{
	int i, j, s = V->size / 2;

// 	omp_set_num_threads(thr);
// #pragma omp parallel shared(B, C) private(i, j)
	for (i = 0; i < s; i++)
	{
		for (j = 0; j < s; j++)
		{
			B->data[i][j] = V->data[i + j];
			C->data[i][j] = V->data[i + j + 1];
		}
	}
}

/**
* Function to compute eigenvalues and eigenvectors with LAPACK functions
* @param 	m		Integer which represent the size of the matrix mat
* @param	mat		Matrix mat where we want to eigenvalues and eigenvectors
* @param	wr		Array with eigenvalues
* @param 	vr		Array 1D with all eigenvectors of size m
*/
void compute_eigenvalues_and_vectors(int m, Matrix *mat, double *wr, double *vr)
{
	int info;
	double wi[m], vl[m];

	/* Solve eigenproblem */
	info = LAPACKE_dgeev(LAPACK_ROW_MAJOR, 'N', 'V', m, mat->data[0], m, wr, wi, vl, m, vr, m);
	/* Check for convergence */
	if (info > 0)
	{
		printf("The algorithm failed to compute eigenvalues.\n");
		exit(1);
	}
}

void step5(int m, Matrix *Xm, Matrix *Vm, Vector *val_ritz, Vector *vect_ritz[m])
{

	int i_max = 0;
	Vector *y = init_vector(m);
	Vector *eigen_vect = init_vector(m * m);
	// Compute eigenvalues and eigenvectors of Xm (eigenvalues == ritz value)
	compute_eigenvalues_and_vectors(m, Xm, val_ritz->data, eigen_vect->data);

	// Compute vectors of ritz
	int i; // PAS DE BESOINS DE PRAGMA CAR m TRES PETIT
    for(i = 0; i < m; i++){
		y->data = &(eigen_vect->data[i*m]);
		prod_mat_vect(Vm, y, vect_ritz[i]); 
	}

	// print_vector(val_ritz);
}

Vector *step6(Matrix *A, int i, Vector *ritz_vectors[i], Vector *val_ritz)
{
	Vector *residus = init_vector(i);
	for (int j = 0; j < i; j++)
	{
		Vector *v1 = init_vector(ritz_vectors[j]->size);
		prod_mat_vect(A, ritz_vectors[j], v1);
		Vector *v2 = scalar_mult_vector(val_ritz->data[j], ritz_vectors[j]);

		residus->data[j] = vect_norm(vector_minus_vector(v1, v2));
		free_vector(v1);
		free_vector(v2);
	}
	return residus;
}

Matrix *convert_vector_array_to_matrix(int m, Vector *y[m])
{
	Matrix *res = init_matrix(y[0]->size, m);
	for (int i = 0; i < res->size[0]; i++)
	{
		for (int j = 0; j < m; j++)
		{
			res->data[i][j] = y[j]->data[i];
		}
	}
	return res;
}

Matrix *extract_matrix_from_mm_file(FILE *f)
{
	int ret_code;
	MM_typecode matcode;
	int M, N, nz;
	int i, I_index, J_index;
	double val;

	if (mm_read_banner(f, &matcode) != 0)
	{
		printf("Could not process Matrix Market banner.\n");
		exit(1);
	}

	if (mm_is_complex(matcode) && mm_is_matrix(matcode) &&
		mm_is_sparse(matcode))
	{
		printf("Sorry, this application does not support ");
		printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
		exit(1);
	}

	if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) != 0)
		exit(1);

	if (M != N)
	{
		printf("Matrix must be symetric.\n");
		exit(1);
	}

	Matrix *A = init_matrix(M, N);

	for (i = 0; i < nz; i++)
	{
		fscanf(f, "%d %d %lf\n", &I_index, &J_index, &val);
		I_index--; /* adjust from 1-based to 0-based */
		J_index--;
		A->data[I_index][J_index] = val;
	}

	return A;
}

Matrix *extract_matrix_from_txt_file(FILE *f)
{
	int M, N;
	double val;

	fscanf(f, "%d %d\n", &M, &N);

	if (M != N)
	{
		printf("Matrix must be symetric.\n");
		exit(1);
	}

	Matrix *A = init_matrix(M, N);

	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			fscanf(f, "%lf ", &val);
			A->data[i][j] = val;
		}
		fscanf(f, "\n");
	}

	return A;
}

Matrix *extract_matrix_from_file(int argc, char **argv)
{
	if (argc < 3)
	{
		fprintf(stderr, "Usage: %s [file-type] [filename]\n", argv[0]);
	}
	FILE *f;
	if ((f = fopen(argv[2], "r")) == NULL)
		exit(1);

	Matrix *A;

	if (strcmp(argv[1], "mm") == 0)
	{
		A = extract_matrix_from_mm_file(f);
	}
	else if (strcmp(argv[1], "txt") == 0)
	{
		A = extract_matrix_from_txt_file(f);
	}
	else
	{
		fprintf(stderr, "Usage: %s [file-type: mm or txt] [filename]\n", argv[0]);
		exit(1);
	}

	printf("Matrix size: %d, %d\n", A->size[0], A->size[1]);

	fclose(f);

	return A;
}

/* Fonction Algorithme itérative PRR */
void PRR(int m, Vector *x, Matrix *A)
{
	// Size and rank of MPI
	int size, rank, iter = 0;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// STEP 1
	int N;
	if (rank == 0){
		N = A->size[0];
		for(int i = 1; i < size; i++){
			MPI_Send(&N, 1, MPI_INT, i, 111, MPI_COMM_WORLD);
		}
	}
	else{
		MPI_Recv(&N, 1, MPI_INT, 0, 111, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	// Normalisation de x + calcul de y0
	double norm, max, C1;
	Vector *y, *residus, *C, *val_ritz;
	Vector *V[m], *vect_ritz[m];
	Matrix *B, *Cc, *Xm, *Vm;
	val_ritz = init_vector(m);

	for (int i = 0; i < m; i++)
	{
		vect_ritz[i] = init_vector(N);
	}

	// START ITERATION
	do {
		if (rank == 0){
			norm = vect_norm(x);
			// printf("%f\n", norm);
			y = normalize(x, norm);
		}
<<<<<<< HEAD
		// print_vector(y);

		C1 = 0;
=======
		else
		{
			y = init_vector(N);
		}
		
>>>>>>> f7f8895b460f3b66c54d5d7d58ec2c94f416663b
		// C0 = || y0 ||^2
		if (rank == 0){
			norm = vect_norm(y);
			C1 = norm * norm;
		}

		// Calcul de C1, C2,....C2m-1
		C = init_vector(2 * m);
		if (rank == 0)
		{
			C->data[0] = C1;
		}
		for (int i = 0; i < m; i++)
		{
			V[i] = init_vector(N);
		}
<<<<<<< HEAD
		step4(C, A, y, m, V, rank);

=======
		// printf("C : \n");
		step4(N, C, A, y, m, V, rank);
		// print_vector(C);

		
>>>>>>> f7f8895b460f3b66c54d5d7d58ec2c94f416663b
		// Calcul de B^ et C^.
		B = init_matrix(m, m);
		Cc = init_matrix(m, m);
		Xm = init_matrix(m, m);

		if (rank == 0)
		{
			fill_B_and_C(B, Cc, C);
			inversion_matrix(B);
			// print_matrix(B)
			prod_mat_mat(B, Cc, Xm);
<<<<<<< HEAD

			// print_vector(Xm);

=======
			// printf("MATRIX \n");
			// print_matrix(Xm);
>>>>>>> f7f8895b460f3b66c54d5d7d58ec2c94f416663b
			// Calcul des valeurs propres et vecteurs propres de Xm
			Vm = convert_vector_array_to_matrix(m, V);
		}

		val_ritz = init_vector(m);
		if (rank == 0)
		{
			step5(m, Xm, Vm, val_ritz, vect_ritz);
		}

		// Test pour la projection ls

		residus = step6(A, m, vect_ritz, val_ritz);
		max = max_in_vector(residus);

<<<<<<< HEAD
		printf("max = %f\n", max);

=======
>>>>>>> f7f8895b460f3b66c54d5d7d58ec2c94f416663b
		// ON RESTART
		int i;
		for (i = 0; i < residus->size; i++)
		{
			if (residus->data[i] == max)
			{
				break;
			}
		}
<<<<<<< HEAD
		// FREE MEMORY
		// print_vector(val_ritz);
		free_vector(y);
		free_matrix(B);
		free_matrix(Cc);
		free_vector(C);
=======
			// FREE MEMORY 
			free_vector(y);
			free_matrix(B);
			free_matrix(Cc);
			free_vector(C);
>>>>>>> f7f8895b460f3b66c54d5d7d58ec2c94f416663b
		iter++;
		x = vect_ritz[i];
		printf("VAL RITZ : \n");
		print_vector(val_ritz);

	}while((max_in_vector(residus) > P) && (iter < MAX_ITER));
	
}
int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);
	int size; MPI_Comm_size(MPI_COMM_WORLD, &size);
	int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	int k, m;
	if (rank == 0){
		printf("Choose number of eigenvalues compute : ");
		scanf("%d", &k);
		m = 2*k; // Size subspace of kylov
		for(int i = 1; i < size; i++){
			MPI_Send(&m, 1, MPI_INT, i, 88, MPI_COMM_WORLD);
		}
	}	
	else {
		MPI_Recv(&m, 1, MPI_INT, 0, 88, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	Vector *v;
	Matrix *A;

	srand(time(0));

	if (rank == 0)
	{
		A = extract_matrix_from_file(argc, argv);
		v = init_vector(A->size[0]);

		fill_vector_with_random_values(v);
		printf("\n");
	}

	// double wr[m];
	// double vr[m*m];
	// compute_eigenvalues_and_vectors(A->size[0], A, wr, vr);

	// for(int i=0; i < m; i++){
	// 	printf("%f ", wr[i]);
	// }
	// printf("\n");

	// for(int i=0; i < m*m; i++){
	// 	printf("%f ", vr[i]);
	// }
	// printf("\n");

	double start, end;

	MPI_Barrier(MPI_COMM_WORLD);
	start = MPI_Wtime();

	PRR(m, v, A);

	MPI_Barrier(MPI_COMM_WORLD);
	end = MPI_Wtime();

	MPI_Finalize();

	if (rank == 0)
	{ /* use time on master node */
		printf("Runtime = %fs\n", end - start);
		free_vector(v);
		free_matrix(A);
	}

	return 0;
}
