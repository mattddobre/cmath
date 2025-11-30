#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include "cmath.h"
#include <math.h>

#define M_PI 3.14159265358979323846

int main(void) {


	
	int k = 4;
	
	cnum ***A_array = cnum_array_init(k, 1);

	A_array[0][0] = cnum_init(-7, 21);
	A_array[1][0] = cnum_init(3, 33);
	A_array[2][0] = cnum_init(40, -5);
	A_array[3][0] = cnum_init(0, -39);

	cnum ***B_array = cnum_array_init(k, 1);
	B_array[0][0] = cnum_init(-7, 21);
	B_array[1][0] = cnum_init(3, 33);
	B_array[2][0] = cnum_init(40, -5);
	B_array[3][0] = cnum_init(0, -39);

	
	cnum ***D_array = cnum_array_init(k, 1);
	D_array[0][0] = cnum_init(.490, 0);
	D_array[1][0] = cnum_init(.672, -.249);
	D_array[2][0] = cnum_init(-.385, -.805);
	D_array[3][0] = cnum_init(-.819, .273);int n = 4;

	cnum ***C_array = cnum_array_init(n, n);

	// --- Proper diagonal ---
	C_array[0][0] = cnum_init(10, 0);
	C_array[1][1] = cnum_init(20, 0);
	C_array[2][2] = cnum_init(5, 0);
	C_array[3][3] = cnum_init(13, 0);

	// Upper triangle chosen using your vector
	C_array[0][1] = cnum_init(-7, 21);    // v0
	C_array[0][2] = cnum_init(3, 33);     // v1
	C_array[0][3] = cnum_init(40, -5);    // v2
	C_array[1][2] = cnum_init(0, -39);    // v3
	C_array[1][3] = cnum_init(2, 11);     // arbitrary
	C_array[2][3] = cnum_init(-8, 4);     // arbitrary

	// Lower triangle = conjugates of upper
	C_array[1][0] = cnum_conj(C_array[0][1]);
	C_array[2][0] = cnum_conj(C_array[0][2]);
	C_array[3][0] = cnum_conj(C_array[0][3]);
	C_array[2][1] = cnum_conj(C_array[1][2]);
	C_array[3][1] = cnum_conj(C_array[1][3]);
	C_array[3][2] = cnum_conj(C_array[2][3]);

	cmat *C = cmat_init(n, n, C_array);

	printf("Hermetian matrix C:\n");
	print_cmat(C);
	

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
	cmat *B = cmat_init(k, 1, B_array);
	cmat *D = cmat_init(k, 1, D_array);
	printf("Matrix A:\n");	
	print_cmat(A);
	printf("Matrix B:\n");	
	print_cmat(B);
	printf("Matrix D:\n");	
	print_cmat(D);

	cmat *B_t = cmat_hermetian(B);
	cmat *Rxx = cmat_mul(A, B_t);
	printf("Correlation matrix Rxx = A * B^H:\n");
	print_cmat(Rxx);

	cnum *scale_factor = cnum_init(1 / 100.0f, 0);
	cmat *Rxx_scaled = cmat_const_mul(scale_factor, Rxx);
	printf("Scaled correlation matrix Rxx (scaled by 0.01):\n");
	print_cmat(Rxx_scaled);
		

	cmat *h = householder_reflection_x(D);
	printf("Householder reflection matrix H:\n");
	print_cmat(h);
	
	cmat *R_ = cmat_mul(h, D);
	printf("Reflected matrix R = H * D:\n");
	print_cmat(R_);

	/*
	cmat *embed = houseolder_embedding(h, 5);
	printf("Embedded Householder matrix:\n");
	print_cmat(embed);	


	cmat *A_copy = cmat_copy(A);
	printf("Copy of matrix A:\n");
	print_cmat(A_copy);


	delete_cmat(A_copy);
	delete_cmat(embed);
	delete_cmat(h);	
	*/

	

	cmat *R  = NULL; 
	cmat *Q = QR_decomposition(C, &R);
	cmat *eigenvals = eigenvalues(C);
	cmat *eigenvects = eigenvectors(C);

	
	printf("Q matrix from QR decomposition:\n");
	print_cmat(Q);
	printf("R matrix from QR decomposition:\n");
	print_cmat(R);

	printf("Eigenvalues of C:\n");
	print_cmat(eigenvals);
	printf("Eigenvectors of C:\n");
	print_cmat(eigenvects);

	float max_theta = find_max_music(k, Rxx_scaled, -M_PI / 2, M_PI / 2, M_PI / 180);
	printf("Estimated angle (degrees): %f\n", max_theta);
	

	delete_cmat(Q);
	delete_cmat(R);	
	delete_cmat(Rxx);
	delete_cmat(A);
	delete_cmat(B);	
	delete_cmat(B_t);
	delete_cmat(D);
	delete_cmat(h);
	delete_cmat(R_);
	delete_cnum(scale_factor);	
	delete_cmat(Rxx_scaled);
	delete_cmat(C);
	delete_cmat(eigenvals);
	delete_cmat(eigenvects);
	

	/*
	float start_angle_rad = -M_PI / 2;
	float end_angle_rad = M_PI / 2;
	float step_rad = M_PI / 180; // 1 degree step
	float max_theta = find_max_music(k, A, start_angle_rad, end_angle_rad, step_rad); 
	printf("Estimated angle (degrees): %f\n", max_theta);
	*/


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
	return 0;
}