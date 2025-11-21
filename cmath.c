#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include "cmath.h"

#define QR_ITER 10 

struct cnum_implementation {
	float a; 
	float b;
};

struct cmat_implementation {
	int n; // number of rows
	int m; // number of columns
	cnum ***cmat;
};

cnum *cnum_init(float a, float b){ 
    cnum *cinst = malloc(sizeof(cnum));
    cinst->a = a; 
    cinst->b = b; 
    return cinst;
}

cnum ***cnum_array_init(float m, float n){ 
	cnum ***array = malloc(m * sizeof(cnum **));
	for (int i = 0; i < m; i++) {
		array[i] = malloc(n * sizeof(cnum *));
	}
	return array;
}

// generate identity matrix of size n x n
cmat *cmat_identity(float n){ 
	cnum ***identity_array = cnum_array_init(n, n);
	for(int i = 0; i < n; i++){ 
		for(int j = 0; j < n; j++){ 
			if(i == j){ 
				identity_array[i][j] = cnum_init(1, 0);
			} else { 
				identity_array[i][j] = cnum_init(0, 0);
			}
		}
	}
	cmat *identity = cmat_init(n, n, identity_array);
	return identity;
}

cmat *cmat_init(float n, float m, cnum ***array){ 
	cmat *cmatinst = malloc(sizeof(cmat));

    cmatinst->n = n; 
    cmatinst->m = m; 	
	cmatinst->cmat = array;
    return cmatinst;
}

cnum *complex_exp(float r, float theta){ 
	return cnum_init(r * cos(theta), r * sin(theta));
}

float cnum_mag(cnum *c1) { 
	return sqrt ( pow(c1->a, 2) + pow(c1->b, 2) );
}

cnum *cadd(cnum *c1, cnum *c2){ 
    return cnum_init(c1->a + c2->a, c1->b + c2->b);
}

cnum *csub(cnum *c1, cnum *c2){ 
    return cnum_init(c1->a - c2->a, c1->b - c2->b);
}

cnum *cmul(cnum *c1, cnum *c2){ 
    return cnum_init(
        c1->a * c2->a - c1->b * c2->b,
        c1->a * c2->b + c1->b * c2->a
    );
}

cnum *cnum_conj(cnum *c) { 
    return cnum_init(c->a, -c->b);
}

cmat *cmat_mul(cmat *a, cmat *b){ 
	// dimensions: a is m x n, b is p x q; require n == p
	if (a->m != b->n) { 
		printf("Dimension mismatch\n");
		return NULL;
	}

	/* result will be a->m x b->n */
	cnum ***result_array = cnum_array_init(a->n, b->m);

	cmat *result = cmat_init(a->n, b->m, result_array);
    
	int i, j, k;

	for(i = 0; i < a->n; i++){ 
		for(j = 0; j < b->m; j++){ 
			result_array[i][j] = cnum_init(0, 0);
			for(k = 0; k < a->m; k++){ 
				cnum *prod = cmul(a->cmat[i][k], b->cmat[k][j]);
				cnum *sum = cadd(result_array[i][j], prod);
				free(result_array[i][j]);
				free(prod);
				result_array[i][j] = sum;
			}
		}
	}

	return result; 
}

cmat *cmat_const_mul(cnum *a, cmat *c){ 
	cnum ***result_array = cnum_array_init(c->n, c->m);

	for(int i = 0; i < c->n; i++){ 
		for(int j = 0; j < c->m; j++){ 
			cnum *prod = cmul(a, c->cmat[i][j]);
			result_array[i][j] = prod;
		}
	}

	cmat *result = cmat_init(c->n, c->m, result_array);
	return result;
	
}

// take norm of nx1 vector
float cmat_norm(cmat *m){ 
	if (m->m != 1) { 
		printf("Input must be a column vector\n");
		return -1;
	}

	float s = 0.0;
	for(int i = 0; i < m->n; i++){ 
		s += pow(cnum_mag(m->cmat[i][0]), 2);
	}

	return sqrt(s);

}

cmat *cmat_normalize(cmat *m){ 
	float norm = cmat_norm(m);
	cnum *norm_cnum = cnum_init(1 / norm, 0);
	cmat *normalized = cmat_const_mul(norm_cnum, m);
	free(norm_cnum);
	return normalized;
}


cmat *cmat_add(cmat *a, cmat *b){ 
	// check dimension
	if (a->m != b->m || a->n != b->n) { 
		printf("Dimension mismatch\n");
		return NULL;
	}

	cnum ***result_array = cnum_array_init(a->n, a->m);

	for(int i = 0; i < a->n; i++){ 
		for(int j = 0; j < a->m; j++){ 
			cnum *cij = cadd(a->cmat[i][j], b->cmat[i][j]);
			result_array[i][j] = cij;
		}
	}

	cmat *result = cmat_init(a->n, a->m, result_array);
	return result; 
}

cmat *cmat_sub(cmat *a, cmat *b){ 
	// check dimension
	if (a->m != b->m || a->n != b->n) { 
		printf("Dimension mismatch\n");
		return NULL;
	}

	cnum ***result_array = cnum_array_init(a->n, a->m);

	for(int i = 0; i < a->n; i++){ 
		for(int j = 0; j < a->m; j++){ 
			cnum *cij = csub(a->cmat[i][j], b->cmat[i][j]);
			result_array[i][j] = cij;
		}
	}

	cmat *result = cmat_init(a->n, a->m, result_array);
	return result; 
}

cmat *cmat_transpose(cmat *m){ 
	cnum ***transposed_array = cnum_array_init(m->m, m->n);	

	int i, j;
	for(i = 0; i < m->n; i++){ 
		for(j = 0; j < m->m; j++){ 
			transposed_array[j][i] = cnum_init(m->cmat[i][j]->a, m->cmat[i][j]->b);
		}
	}

	cmat *transposed = cmat_init(m->m, m->n, transposed_array);
	return transposed;
}

cmat *cmat_hermetian(cmat *m){ 
	cnum ***transposed_array = cnum_array_init(m->m, m->n);	

	int i, j;
	for(i = 0; i < m->n; i++){ 
		for(j = 0; j < m->m; j++){ 
			cnum *cji = cnum_init(m->cmat[i][j]->a, m->cmat[i][j]->b);
			transposed_array[j][i] = cnum_conj(cji);
			//free(cji);
		}
	}

	cmat *transposed = cmat_init(m->m, m->n, transposed_array);
	return transposed;
}



void delete_cmat(cmat *m){ 
	for(int i = 0; i < m->n; i++){ 
		for(int j = 0; j < m->m; j++){ 
			free(m->cmat[i][j]);
		}
		free(m->cmat[i]);
	}
	free(m->cmat);
	free(m);
}

void print_cnum(cnum *c, int newline){ 
    if (newline)
        printf("%f + %fj\n", c->a, c->b);
    else
        printf("%f + %fj ", c->a, c->b);
}

void print_cmat(cmat *m){ 
    for(int i = 0; i < m->n; i++){ 
        for(int j = 0; j < m->m; j++){ 
            print_cnum(m->cmat[i][j], 0);
        }
        printf("\n");
    }
}

// computes proejction p that reflects vector about a plane
// projects a nx1 vector onto another nx1 vector b that is parallel 
// to the x axis or the (1, 0, 0 ... ) axis, to be used in QR decomposition

cmat* householder_reflection_x(cmat *m){ 

	cmat *I = cmat_identity(m->n);
	// construct vector b which is just 1, 0, 0, ...
	cmat *b = cmat_init(m->n, 1, cnum_array_init(m->n, 1));
	b->cmat[0][0] = cnum_init(1, 0);
	for(int i = 1; i < m->n; i++){ 
		b->cmat[i][0] = cnum_init(0, 0);
	}

	cnum *norm_m = cnum_init(cmat_norm(m), 0);
	float mag = cnum_mag(m->cmat[0][0]);
	cnum *beta = cnum_init(m->cmat[0][0]->a/mag, m->cmat[0][0]->b / mag);

	cmat *scaled_b = NULL;

	if (m->cmat[0][0]->a == 0 && m->cmat[0][0]->b == 0) {
		scaled_b = cmat_const_mul(norm_m, b);
	}else { 
		cnum *coeff = cmul(beta, norm_m);
		scaled_b = cmat_const_mul(coeff, b); // add in appropriate scale to the b vector
		free(coeff);
	}


	cmat *u = NULL;
	if (m->cmat[0][0]->a >= 0 && m->cmat[0][0]->b == 0) { 
		u = cmat_add(m, scaled_b);
	} else { 
		u = cmat_sub(m, scaled_b);
	}

	cmat *n = cmat_normalize(u);
	cmat *n_t = cmat_hermetian(n);

	cmat *outer_product = cmat_mul(n, n_t);
	cnum *two = cnum_init(2, 0);
	cmat *scaled_outer = cmat_const_mul(two, outer_product);
	cmat *result = cmat_sub(I, scaled_outer);

	// free what is used for computation
	delete_cmat(scaled_b);
	free(norm_m);
	delete_cmat(b);
	delete_cmat(u);
	delete_cmat(n);
	delete_cmat(n_t);
	delete_cmat(outer_product);
	free(two);
	delete_cmat(I);
	delete_cmat(scaled_outer);

	
	return result;

}

// embed Householder reflection matrix p of size nxn into the lower part of an identity matrix of size mxm
cmat *houseolder_embedding(cmat *p, int m){ 
	if(m < p->n){ 
		printf("Identity dimension must be larger than Householder matrix dimension\n");
		return NULL;
	}

	cmat *I = cmat_identity(m);

	int i, j; 
	for(i = 0; i < p->n; i++){ 
		for(j = 0; j < p->n; j++){ 
			//free(I->cmat[m - p->n + i][m - p->n + j]);
			I->cmat[m - p->n + i][m - p->n + j] = cnum_init(p->cmat[i][j]->a, p->cmat[i][j]->b);
		}
	}
	return I;
}

cmat *cmat_copy(cmat *m) {
	cnum ***copy_array = cnum_array_init(m->n, m->m);
	for (int i = 0; i < m->n; i++) {
		for (int j = 0; j < m->m; j++) {
			copy_array[i][j] = cnum_init(m->cmat[i][j]->a, m->cmat[i][j]->b);
		}
	}
	cmat *copy = cmat_init(m->n, m->m, copy_array);
	return copy;
}

cmat *QR_decomposition(cmat *A, cmat **R_out){ 
	int n = A->n;
	cnum ***A_copy = cmat_copy(A)->cmat;
	
	cmat *R = cmat_init(A->n, A->m, A_copy);
	cmat *Q_total = cmat_identity(A->n);

	for(int k = 0; k < n-1; k++){ 
		// extract subvector from R
		cnum ***subvector_array = cnum_array_init(R->n - k, 1);
		for (int i = k; i < R->n; i++) {
			/* deep-copy each element into the subvector to avoid aliasing */
			subvector_array[i - k][0] = cnum_init(R->cmat[i][k]->a, R->cmat[i][k]->b);
		}
		cmat *subvector = cmat_init(R->n - k, 1, subvector_array);

		cmat *P_k = householder_reflection_x(subvector);
		cmat *E_k = houseolder_embedding(P_k, A->n);

		cmat *R_new = cmat_mul(E_k, R);
		cmat *Q_total_new = cmat_mul(Q_total, E_k);

		delete_cmat(R);
		delete_cmat(Q_total);
		delete_cmat(subvector);
		delete_cmat(P_k);
		delete_cmat(E_k);

		R = R_new;
		Q_total = Q_total_new;
	}

	*R_out = R;
	return Q_total;
}

// return eigenvalues as computed using QR decomposition
// oh we're modifying A
cmat *eigenvalues(cmat *A, int max_iters){ 
	int i; 
	cmat *A_0 = cmat_copy(A);
	for (i = 0; i < max_iters; i++) {
		cmat *R = NULL;
		cmat *Q = QR_decomposition(A_0, &R);

		delete_cmat(A_0);
		cmat *A_next = cmat_mul(R, Q);
		A_0 = A_next;
		delete_cmat(R);
		delete_cmat(Q);

	}
	// placeholder function
	return A_0;
}


// return eigenvectors as computed using QR decomposition
cmat *eigenvectors(cmat *A, int max_iters) {
	int i; 
	cmat *A_0= cmat_copy(A);
	cmat *E_0 = cmat_identity(A->n);
	for (i = 0; i < max_iters; i++) {
		cmat *R = NULL;
		cmat *Q = QR_decomposition(A_0, &R);

		cmat *A_next = cmat_mul(R, Q);
		cmat *E_next = cmat_mul(E_0, Q);
		
		delete_cmat(A_0);
		delete_cmat(E_0);
		A_0 = A_next; 
		E_0 = E_next;

		delete_cmat(R);
		delete_cmat(Q);

	}
	// placeholder function
	return E_0;
}

int find_largest_eigenvalue_index(cmat *eigvals) { 
	int idx = 0;
	float min_mag = cnum_mag(eigvals->cmat[0][0]);
	for (int j = 1; j < eigvals->m; j++) { 
		float mag = cnum_mag(eigvals->cmat[0][j]);
		if (mag > min_mag) { 
			min_mag = mag;
			idx = j;
		}
	}
	return idx;
}

cmat *cmat_remove_lowest_eigenvector(int idx, cmat *eigenvects){ 
	// placeholder function
	cnum ***eigenvects_reduced = cnum_array_init(eigenvects->n, eigenvects->m - 1);
	cmat *eigenvects_reduced_mat = cmat_init(eigenvects->n, eigenvects->m - 1, eigenvects_reduced);

	int i, j; 
	for (i = 0; i < eigenvects->n; i++){ 
		int col_idx = 0;
		for (j = 0; j < eigenvects->m; j++){ 
			if (j == idx) {
				continue;
			}
			eigenvects_reduced[i][col_idx] = cnum_init(eigenvects->cmat[i][j]->a, eigenvects->cmat[i][j]->b);
			col_idx++;
		}
	}
	return eigenvects_reduced_mat;
}

// assumes a lambda/2 spacing
cmat *beamforming_array(int num_elements, float angle_rad) {
	cnum ***array = cnum_array_init(num_elements, 1);
	cmat *bf_array = cmat_init(num_elements, 1, array);
	float d = 0.5; // element spacing in wavelengths
	for (int n = 0; n < num_elements; n++) { 
		float phase_shift = 2 * (3.141592653589793) * d * n * sin(angle_rad);
		array[n][0] = complex_exp(1.0, phase_shift);
	}


	return bf_array;

}

float music_algorithm_theta(int num_elements, cmat *X, float theta_rad) { 
	printf("1\n");
	cmat *bf_theta = beamforming_array(num_elements, theta_rad);
	printf("2\n");
	cmat *bf_theta_herm = cmat_hermetian(bf_theta);	
	printf("3\n");

	cmat *X_herm = cmat_hermetian(X);
	printf("4\n");

	cmat *Rxx = cmat_mul(X, X_herm);
	printf("5\n");
	cmat *eigvect = eigenvectors(Rxx, QR_ITER);
	printf("6\n");
	cmat *eigvals = eigenvalues(Rxx, QR_ITER);
	printf("7\n");

	int largest_eig_idx = find_largest_eigenvalue_index(eigvals);
	cmat *eigvect_reduced = cmat_remove_lowest_eigenvector(largest_eig_idx, eigvect);
	printf("8\n");
	cmat *eigvect_reduced_herm = cmat_hermetian(eigvect_reduced);
	printf("9\n");
	cmat *denom_mat = cmat_mul(bf_theta_herm, eigvect_reduced);
	printf("10\n");
	denom_mat = cmat_mul(denom_mat, eigvect_reduced_herm);
	printf("11\n");
	denom_mat = cmat_mul(denom_mat, bf_theta);		
	printf("12\n");
	float denom_mag = cnum_mag(denom_mat->cmat[0][0]);
	float P_music = 1 / denom_mag;

	// free everything not used for computation
	free(denom_mat);
	delete_cmat(bf_theta);
	delete_cmat(bf_theta_herm);
	delete_cmat(X_herm);
	delete_cmat(Rxx);
	delete_cmat(eigvect);
	delete_cmat(eigvals);
	delete_cmat(eigvect_reduced);
	delete_cmat(eigvect_reduced_herm);	
	printf("13\n");

	return P_music;
}

float find_max_music(int num_elements, cmat *X, float start_angle_rad, float end_angle_rad, float step_rad) { 
	float max_P = -1.0;
	float max_angle = start_angle_rad;
	for (float theta = start_angle_rad; theta <= end_angle_rad; theta += step_rad) { 
		float P_music = music_algorithm_theta(num_elements, X, theta);
		if (P_music > max_P) { 
			max_P = P_music;
			max_angle = theta;
		}
	}
	float max_angle_deg = 180 * max_angle / 3.141592653589793;
	printf("Max MUSIC P: %f at angle (deg): %f\n", max_P, max_angle_deg);
	return max_angle_deg;
}




