#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>
#include <math.h>
#include <lapacke.h>

//resolution du système de yule-walker système
void solveLS(double *a, int m, double *b, int nrhs)
{
    int lda = m;
    int *y = calloc(m, sizeof(int));
    int info;
    LAPACK_dgetrf(&m, &m, a, &lda, y, &info); //la factorisation de la matrice
    char trans = 'N';
    int ldb = m;
    LAPACK_dgetrs(&trans, &m, &nrhs, a, &lda, y, b, &ldb, &info); //la résolution du système linéaire
    free(y);
}

void computeEV(double *B, int n, double *wr, double *wi, double *vr)
{
    char balanc = 'B';
    char jobvl = 'N';
    char jobvr = 'V';
    char sense = 'V';
    //int n = m;
    int ldb = n;
    //double* wr = calloc(n,sizeof(double));
    //double* wi = calloc(n,sizeof(double));
    int ldvl = n;
    double *vl = calloc(ldvl * n, sizeof(double));
    int ldvr = n;
    //double* vr = calloc(ldvr*n,sizeof(double));
    int ilo = 0, ihi = 0;
    double *scale = calloc(n, sizeof(double));
    double abnrm = 0.;
    double *rconde = calloc(n, sizeof(double));
    double *rcondv = calloc(n, sizeof(double));
    int *iwork = calloc(2 * n - 2, sizeof(int));
    double wkopt = 0.;
    int lwork = -1;
    int info;

    LAPACK_dgeevx(&balanc, &jobvl, &jobvr, &sense, &n, B, &ldb, wr, wi, vl, &ldvl, vr, &ldvr, &ilo, &ihi, scale, &abnrm, rconde, rcondv, &wkopt, &lwork, iwork, &info);
    //LAPACK_dgeev(&jobvl,&jobvr,&n,B,&ldb,wr,wi,vl,&ldvl,vr,&ldvr,&wkopt,&lwork,&info);
    lwork = (int)wkopt;
    double *work = calloc(lwork, sizeof(double));
    LAPACK_dgeevx(&balanc, &jobvl, &jobvr, &sense, &n, B, &ldb, wr, wi, vl, &ldvl, vr, &ldvr, &ilo, &ihi, scale, &abnrm, rconde, rcondv, work, &lwork, iwork, &info);
    //LAPACK_dgeev(&jobvl,&jobvr,&n,B,&ldb,wr,wi,vl,&ldvl,vr,&ldvr,work,&lwork,&info);
    if (info > 0)
    {
        printf("The algorithm failed to compute eigenvalues.\n");
        exit(1);
    }
    //valeurs propres
    //print_valeurs("valeurs propres de H",n,wr,wi);
    //vecteurs propres
    //	print_vectors("vecteurs propres approchés de H(droite)",n,wi,vr,ldvr);
    free(vl);
    free(scale);
    free(rconde);
    free(rcondv);
    free(iwork);
    free(work);
}

void computeSortEV(double *B, int n, double *wr, double *wi, double *vs)
{
    char jobvs = 'V';
    char sort = 'S';
    LAPACK_D_SELECT2 select;
    char sense = 'V';
    int lda = n;
    int sdim;
    int ldvs = n;
    double rconde, rcondv;
    double wkopt = 0.;
    int iwkopt = 0;
    int lwork = -1;
    int liwork = -1;
    int *bwork = calloc(n, sizeof(int));
    int info;

    LAPACK_dgeesx(&jobvs, &sort, select, &sense, &n, B, &lda, &sdim, wr, wi, vs, &ldvs, &rconde, &rcondv, &wkopt, &lwork, &iwkopt, &liwork, bwork, &info);
    lwork = (int)wkopt; // printf("lwork = %d",lwork);
    double *work = calloc(lwork, sizeof(double));
    liwork = iwkopt; //printf("liwork = %d",liwork);
    int *iwork = calloc(liwork, sizeof(int));
    LAPACK_dgeesx(&jobvs, &sort, select, &sense, &n, B, &lda, &sdim, wr, wi, vs, &ldvs, &rconde, &rcondv, work, &lwork, iwork, &liwork, bwork, &info);
    if (info > 0)
    {
        printf("The algorithm failed to compute eigenvalues.\n");
        exit(1);
    }
    //print_valeurs("valeurs propres de H",n,wr,wi);
    //print_vectors("vecteurs propres approchés de H(droite)",n,wi,vs,ldvs);
    free(bwork);
    free(work);
    free(iwork);
}