#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include "cmath.h"
#include <math.h>

#define M_PI 3.14159265358979323846

int main(void) {


	
	int k = 4;
	
	cnum *c = cnum_init(4, 4);
	cnum *dest = cnum_init(0,0);
	cnum_conj(c, dest);
	print_cnum(dest, 1);
	cmul(c, dest, dest);
	print_cnum(c, 1);
	print_cnum(dest, 1);

	cnum ***A_array = cmat_array_init(k, 1);
	A_array[0][0] = cnum_init(-64, -27);
	A_array[1][0] = cnum_init(-18, -25);
	A_array[2][0] = cnum_init(26, 8);
	A_array[3][0] = cnum_init(-38, 23);
	cmat *A = cmat_init(k, 1, A_array);
	printf("Matrix A:\n");
	print_cmat(A);

	cnum ***B_array = cmat_array_init(k, 1);
	B_array[0][0] = cnum_init(-1, 0);
	B_array[1][0] = cnum_init(2, 3);
	B_array[2][0] = cnum_init(4, 5);
	B_array[3][0] = cnum_init(6, 7);
	cmat *B = cmat_init(k, 1, B_array);

	cnum ***C_array = cmat_array_init(k, 1);
	cmat *C = cmat_init(k, 1, C_array);

	cnum ***D_array = cmat_array_init(k, k);
	cmat *D = cmat_init(k, k, D_array);

	
	
	cnum ***Rxx_array = cmat_array_init(k, k);
	cmat *Rxx = cmat_init(k, k, Rxx_array);
	cnum ***Rxx_array_scaled = cmat_array_init(k, k);
	cmat *Rxx_scaled = cmat_init(k, k, Rxx_array_scaled);
	cmat_hermetian(A, C);
	cmat_mul(A, C, Rxx);
	cmat_const_mul(cnum_init(1.0/100, 0), Rxx, Rxx_scaled);
	printf("Matrix Rxx_scaled (A * A^H):\n");
	print_cmat(Rxx_scaled);


	/*
	float max_theta = find_max_music(k, Rxx_scaled, -M_PI / 2, M_PI / 2, M_PI / 180);
	printf("Estimated angle (degrees): %f\n", max_theta);
	*/
	
	beamforming_array(k, M_PI/6, B);
	printf("Beamforming vector at 30 degrees:\n");
	print_cmat(B);

	float music_val = music_algorithm_theta(k, Rxx_scaled, B);
	printf("MUSIC algorithm output P: %f\n", music_val);
	
	



	/*
	cnum ***R_array = cmat_array_init(k, k);
	cmat *R = cmat_init(k, k, R_array);

	cnum ***Q_array = cmat_array_init(k, k);
	cmat *Q = cmat_init(k, k, Q_array);

	QR_decomposition(Rxx_scaled, R, Q);

	printf("Matrix Rxx_scaled (A * A^H):\n");
	print_cmat(Rxx_scaled);
	printf("Matrix R (from QR Decomposition):\n");
	print_cmat(R);		
	printf("Matrix Q (from QR Decomposition):\n");
	print_cmat(Q);

	cmat_mul(Q, R, D);
	printf("Product Q * R:\n");
	print_cmat(D);
	*/
	

	/*
	eigenvalues(Rxx_scaled, D);
	printf("Eigenvalues of Rxx_scaled:\n");
	print_cmat(D);

	eigenvectors(Rxx_scaled, D);
	printf("Eigenvectors of Rxx_scaled:\n");
	print_cmat(D);
	
	
	
	
	householder_reflection_x(B, D);
	printf("Householder Reflection H:\n");
	print_cmat(D);

	cmat_mul(D, B, C);
	printf("H * B:\n");
	print_cmat(C);
	*/
	


	/*
	cnum ***T_array = cmat_array_init(5, 5);
	cmat *T = cmat_init(5, 5, T_array);
	
	householder_embedding(D, 5, T);
	printf("Householder Embedding P:\n");
	print_cmat(T);
	

	
	cmat_hermetian(A, C);
	cmat_hermetian(C, B);
	printf("Hermetian of A:\n");
	print_cmat(C);

	cmat_mul(A, C, D);
	printf("A * A^H:\n");
	print_cmat(D);

	cmat_add(A, A, B);
	printf("A + A:\n");
	print_cmat(B);
	

	cmat_sub(A, A, B);
	printf("A - A:\n");
	print_cmat(B);
	*/
	
	




		
	


	return 0;
}