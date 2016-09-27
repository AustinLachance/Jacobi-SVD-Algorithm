/*****************************************************************
Name: Austin Lachance                                                            
Email: austin.lachance@yale.edu                                                  
Date: 02/26/16                                                                   
                                                                                 
CPSC440                                                                          
Jacobi Algorithm SVD                                                
Description: This program performs SVD on a matrix a
using the Jacobi Method                                                                    
*****************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void jacobi(double * a, int n, double * s, double * u, double * v);
void *identity(double* matrix, int n);
void *transpose(double *matrix, int n);
void matrix_print(double *matrix, int n);
int error_check(double a, double b, double error);
int is_diagonal(double *matrix, int n);
void neg_check(double * a, double* u, double *s, int n);
void column_change(double * matrix, int col1, int col2, int n);
void sort(double* s, double *a, double *u, double *v, int n);


int main(int argc, char **argv) {
  int i, k;
  int n = 3;
  double *a = malloc(sizeof(double) * n * n); //Matrix to be diagonalized
  double *s = malloc(sizeof(double) * n); //Spectrum of a
  double *u = calloc(n*n, sizeof(double)); //Left singular values
  double *v = calloc(n*n, sizeof(double)); //Right singular values

  a[0] = -3.4;
  a[1] = 2;
  a[2] = 4.56;
  a[3] = -2;
  a[4] = 2;
  a[5] = 1.89;
  a[6] = 0;
  a[7] = 4;
  a[8] = 10000.34;
  
  jacobi(a, n, s, u, v);

  matrix_print(u, n); 
  printf("\n");
  matrix_print(v, n); 
  printf("\n");
  for(i = 0; i < n; i++) {
  	printf("%f\t", s[i]);
  }

  free(a); 
  free(s); 
  free(u); 
  free(v); 
  return 0;
}


//Perform Jacobi algorithm to find SVD of matrix a of size n
void jacobi(double * a, int n, double * s, double * u, double * v) {
	double temp, temp1, temp2, rot1, rot2;
	int row, col, i, j;
	int num_iterations = 0;
	identity(v, n);
	identity(u, n);

	while((num_iterations < 1) || !is_diagonal(a, n)) {

		for(row = 0; row < n-1; row++) {
			for(col = row + 1; col < n; col++) {
				
				temp1 = atan((a[col*n + row]
				 - a[row*n + col])/(a[row*n + row] + a[col*n + col]));
				temp2 = atan((a[row*n + col]
				 + a[col*n + row])/(a[row*n + row] - a[col*n + col]));

				//Rotation angles
				rot1 = 0.5*(temp1 + temp2);
				rot2 = 0.5*(temp2 - temp1);

				//Transform rows & cols of a,u and a,v
				for(i = 0; i < n; i++) {
					temp = a[row*n + i];
					a[row*n + i] = a[row*n + i]*cos(rot1)
					 + a[col*n + i]*sin(rot1);
					a[col*n + i] = a[col*n + i]*cos(rot1)
					 - temp*sin(rot1);

					 temp = u[row*n + i];
					u[row*n + i] = u[row*n + i]*cos(rot1)
					 + u[col*n + i]*sin(rot1);
					u[col*n + i] = u[col*n + i]*cos(rot1)
					 - temp*sin(rot1);
				}
                
				for(i = 0; i < n; i++) {
					temp = a[i*n + row];
					a[i*n + row] = a[i*n + row]*cos(rot2)
					 + a[i*n + col]*sin(rot2);
					a[i*n + col] = a[i*n + col]*cos(rot2)
					 - temp*sin(rot2);

					 temp = v[i*n + row];
					v[i*n + row] = v[i*n + row]*cos(rot2)
					 + v[i*n + col]*sin(rot2);
					v[i*n + col] = v[i*n + col]*cos(rot2)
					 - temp*sin(rot2);
				}
			}
		}
		num_iterations++;
	}

	//transpose U
	transpose(u,n);
	
	//make sure s isn't negative. Fix columns in other matrices accordingly
	neg_check(a, u, s, n);

	//sort s and adjust columns in other matrices
	sort(s, a, u, v, n);	
}


//Create an indentity matrix "matrix" of size n
void *identity(double* matrix, int n) {
  int i = 0;
  for(i = 0; i < n; i++) 
    {
      matrix[i*n + i] = (double)1;
    }
}


//Transpose matrix of size n
void *transpose(double *matrix, int n) {
  double *transpose = (double*)malloc(n*n*sizeof(double));
  int i, j = 0;
  for(i = 0; i < n; i++) {
  	for(j = 0; j < n; j++) {
	  transpose[i*n + j] = matrix[j*n + i];
	}
  }
  for(i = 0; i < n*n; i++) {
  	matrix[i] = transpose[i];
  }
  free(transpose);
}


//Prints matrix of size n
void matrix_print(double *matrix, int n) {
  int i, j = 0;
  for(i = 0; i < n; i++) {
  	for(j = 0; j < n; j++) {
	  printf("%e\t", matrix[i*n + j]);
	}
	printf("\n");
  }
}


//Checks if a and b are within error of one another
int error_check(double a, double b, double error) {
  return (fabs(a-b) < error);
}


//Determines if matrix of size n is diagonal by ensuring only diags are > 0.0
int is_diagonal(double *matrix, int n) {
  int row, col;
  for(row=0; row<n; row++) {
    for(col=0; col<n; col++) {
      if(row != col) {
        if(!error_check(matrix[row*n + col], 0.0, 1e-16)) {
        	return 0;
        }
      }
    }
  }
  return 1;
}


//Makes all values in s positive and alters the a and u matrices accordingly
void neg_check(double * a, double* u, double *s, int n) {
	int i, j;
	for(i = 0; i < n; i++) {
		if(a[i*n + i] < 0.0) {
			a[i] = -1.0*a[i*n + i];
			for(j = 0; j < n; j++) {
				u[j*n + i] = -1.0*u[j*n + i];
			}
		}
		s[i] = a[i*n + i];
	}
}


//swaps col1 and col2 in matrix of size n
void column_change(double * matrix, int col1, int col2, int n) {
	int i, j;
	double temp;
	for(i = 0; i < n; i++) {
		temp = matrix[i*n + col1];
		matrix[i*n + col1] = matrix[i*n + col2];
		matrix[i*n + col2] = temp;
	}
}


//sorts array s or size n and swaps appropriate cols in a, u, and v
void sort(double* s, double *a, double *u, double *v, int n) {
	double temp=0;
	int i, j;
	j=0;
	for(i = 1; i < n; i++) 
	{
		temp = s[i];
		j = i-1;
		while(j>=0 && (temp > s[j])) 
		{
			s[j+1] = s[j];
			column_change(a, j+1, j, n);
			column_change(u, j+1, j, n);
			column_change(v, j+1, j, n);
			j--;
		}
		s[j+1] = temp;
	}
}
