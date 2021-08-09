#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#define BLOCKSIZE 64

void do_block(int n, int si, int sj, int sk, double *A, double *B, double *C) {
    for (int i = si; i < si + BLOCKSIZE; i++) {
        for (int j = sj; j < sj + BLOCKSIZE; j++) {
            double cij = C[i + j * n]; /*Cij=C[i][j]*/
            for (int k = sk; k < sk + BLOCKSIZE; k++) {
                cij += A[i + k * n] * B[k + j * n];
            }
            C[i + j * n] = cij;
        }
    }
}

void blocked_dgemm(int n, double *A, double *B, double *C) {
    for (int sj = 0; sj < n; sj += BLOCKSIZE)
        for (int si = 0; si < n; si += BLOCKSIZE)
            for (int sk = 0; sk < n; sk += BLOCKSIZE)
                do_block(n, si, sj, sk, A, B, C);
}

void dgemm(int n, double *A, double *B, double *C) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double cij = C[i + j * n];
            for (int k = 0; k < n; k++) {
                cij += A[i + k * n] * B[k + j * n];
            }
            C[i + j * n] = cij; /*C[i][j] = cij */
        }
    }
}

int main() {
    int n = 960;
    double *A;
    double *B;
    double *C;
    struct timeval start, end;
    long time = 0;
    A = (double *) malloc(n * sizeof(double));
    B = (double *) malloc(n * sizeof(double));
    C = (double *) malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        A[i] = 1;
        B[i] = 2;
        C[i] = 0;
    }

    printf("For block size = %d and matrices %d x %d: \n",BLOCKSIZE,n,n);

    gettimeofday(&start, NULL);
    dgemm(n, A, B, C);
    gettimeofday(&end, NULL);
    time = end.tv_sec - start.tv_sec;
    time = ((time * 1000000) + end.tv_usec) - (start.tv_usec);
    printf("Naive matrix multiply time = %ld microseconds\n", time);

    gettimeofday(&start, NULL);
    blocked_dgemm(n, A, B, C);
    gettimeofday(&end, NULL);
    time = end.tv_sec - start.tv_sec;
    time = ((time * 1000000) + end.tv_usec) - (start.tv_usec);
    printf("Cache blocked matrix multiply time = %ld microseconds\n", time);
    return 0;
}
