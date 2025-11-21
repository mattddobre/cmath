#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include "cmath.h"
#include <math.h>

int main(void) {


	
	int k = 4;
	
	cnum ***A_array = cnum_array_init(k, 1);

	A_array[0][0] = cnum_init(-7, 21);
	A_array[1][0] = cnum_init(3, 33);
	A_array[2][0] = cnum_init(40, -5);
	A_array[3][0] = cnum_init(0, -39);

	/*
	A_array[1][0] = cnum_init(1, -2);  // conjugate of A[0][1]
	A_array[1][1] = cnum_init(4, 0);
	A_array[1][2] = cnum_init(0.5, 3);
	A_array[1][3] = cnum_init(-1, 2);

	A_array[2][0] = cnum_init(-2, 1);  // conjugate of A[0][2]
	A_array[2][1] = cnum_init(0.5, -3); // conjugate of A[1][2]
	A_array[2][2] = cnum_init(5, 0);
	A_array[2][3] = cnum_init(1, -1);

	A_array[3][0] = cnum_init(0, -1);  // conjugate of A[0][3]
	A_array[3][1] = cnum_init(-1, -2); // conjugate of A[1][3]
	A_array[3][2] = cnum_init(1, 1);   // conjugate of A[2][3]
	A_array[3][3] = cnum_init(2, 0);
	*/
	


	cmat *A = cmat_init(k, 1, A_array);
	cmat *A = cmat_init(k, 1, A_array);
	printf("Matrix A:\n");	
	print_cmat(A);
	float start_angle_rad = -M_PI / 2;
	float end_angle_rad = M_PI / 2;
	float step_rad = M_PI / 180; // 1 degree step
	float max_theta = find_max_music(k, A, start_angle_rad, end_angle_rad, step_rad); 
	printf("Estimated angle (degrees): %f\n", max_theta);

	/*
	
	cmat *R = NULL;
	cmat *Q = QR_decomposition(A, &R);

	printf("Q matrix from QR decomposition:\n");
	print_cmat(Q);
	printf("R matrix from QR decomposition:\n");
	print_cmat(R);

	cmat *prod = cmat_mul(Q, R);
	printf("Product Q*R:\n");
	print_cmat(prod);
	
	
	// eigenvalues of A
	
	cmat *eigvals = eigenvalues(A, 100);
	printf("Eigenvalues of A:\n");
	print_cmat(eigvals);

	cmat *eigenvects = eigenvectors(A, 100);
	printf("Eigenvectors of A:\n");
	print_cmat(eigenvects);
	*/
	delete_cmat(A);


	return 0;
}