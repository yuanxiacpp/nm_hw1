#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void printMatrix(double *a, int n) {
  printf("***************** Matrix %d x %d *********************\n", n, n);
  int i = 0;
  for (i = 0; i < n * n; ++i) {
    printf("%20.15f ", a[i]);
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

void normalizeVector(double *a, int n) {
	double sum = 0;
	int i;
	for (i = 0; i < n; i++)
		sum += a[i] * a[i];
	sum = sqrt(sum);

	if (sum == 0)
		return;
	for (i = 0; i < n; i++)
		a[i] /= sum;
	return;
}

void inv_double_gs(double *a, int n, double *u, double *b) {
	return;
}
void problem1(int n) {
	int i, j;

	double *a = (double *)malloc(n*n*sizeof(double));

	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			a[i*n+j] = (double)rand()/(double)RAND_MAX;


  
}

int main() {
	
	return 0;
}