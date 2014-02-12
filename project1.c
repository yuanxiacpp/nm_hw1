#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

void printMatrix(double *a, int n) {
  printf("***************** Matrix %d x %d *********************\n", n, n);
  int i = 0;
  for (i = 0; i < n * n; ++i) {
    printf("%8.5f ", a[i]);
    if ((i+1) % n == 0)
      printf("\n");
  }
  return;
}

void transpose(double *a, int n) {
  int i, j;
  for (i = 0; i < n; ++i) {
    for (j = i + 1; j < n; ++j) {
      double tmp = a[i*n + j];
      a[i*n + j] = a[j*n + i];
      a[j*n + i] = tmp;
    }
  }
  return;
}

double* multiply(double *a, double *b, int p, int q, int r) {
  double *result = (double*)malloc(p*r*sizeof(double));
  int i, j, k;
  for (i = 0; i < p; i++) {
    for (j = 0; j < r; j++) {
      double sum = 0;
      for (k = 0; k < q; k++)
        sum += a[i*q+k] * b[k*r+j];
      result[i*r+j] = sum;
    }
  }
  return result;
}

void normalizeCol(double *a, int n, int col) {
  double sum = 0;
  int i;
  for (i = 0; i < n; i++) {
    sum += a[n*i+col] * a[n*i+col];
  }
  sum = sqrt(sum);

  if (sum == 0)
    return;
  for (i = 0; i < n; i++)
    a[n*i+col] /= sum;
  return;
}
double dotProductForCols(double *a, int n, int col1, int col2) {
  int i;
  double result = 0.0;
  for (i = 0; i < n; i++)
    result += a[n*i+col1]*a[n*i+col2];
  return result;
}
void updateFromPrevCols(double *u, int n, int col) {
  int i, j;
  double* tmp = (double *)malloc(n*sizeof(double));
  for (i = 0; i < n; i++)
    tmp[i] = 0.0;
  
  //sum of all prev dot products
  for (i = 0; i < col; i++) {
    double dot = dotProductForCols(u, n, col, i);
    for (j = 0; j < n; j++) 
      tmp[j] += dot * u[j*n+i];
  }

  //update current col
  for (i = 0; i < n; i++)
    u[i*n+col] -= tmp[i];
  
  return;
}


void updateFollowingCols(double *u, int n, int col) {
  int i, j;
  for (j = col; j < n; j++) {
    double dot = dotProductForCols(u, n, j, col-1);
    //update j-th col with (col-1)-th unit vector
    for (i = 0; i < n; i++)
      u[i*n+j] -= dot * u[i*n+col-1];
  }

  return;
}

void inv_double_gs(double *a, int n, double *u, double *b) {
  
  int i, j;
  for (i = 0; i < n; i++) {
    updateFromPrevCols(u, n, i);
    normalizeCol(u, n, i);
    updateFollowingCols(u, n, i+1);
  }


  //printMatrix(a, n);
  //printMatrix(u, n);

  double* ut = (double*)malloc(n*n*sizeof(double));
  memcpy(ut, u, n*n*sizeof(double));

  //transpose(ut, n);
  
  printMatrix(multiply(u, ut, n, n, n), n);
    

  
  return;
}
void problem1(int n) {
  int i, j;

  double *a = (double *)malloc(n*n*sizeof(double));
  double *u = (double *)malloc(n*n*sizeof(double));
  double *b = (double *)malloc(n*n*sizeof(double));

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      a[i*n+j] = (double)rand()/(double)RAND_MAX;
      //a[i*n+j] = rand() % 10;

  //printMatrix(a, n);
  
  memcpy(u, a, n*n*sizeof(double));
  
  inv_double_gs(a, n, u, b);

}

int main() {
  problem1(10);
  return 0;
}
