#ifndef CMATH_POOL_H
#define CMATH_POOL_H

#include <stdint.h>

// Forward declarations
typedef struct cnum_implementation cnum;
typedef struct cmat_implementation cmat;

// -------------------- cnum --------------------
cnum *cnum_init(float a, float b);
float cnum_mag(cnum *c);
float cnum_mag_squared(cnum *c);

int cadd(cnum *c1, cnum *c2, cnum *dest);
int csub(cnum *c1, cnum *c2, cnum *dest);
int cmul(cnum *c1, cnum *c2, cnum *dest);
int cnum_conj(cnum *c, cnum *dest);

// -------------------- cmat --------------------
cnum ***cmat_array_init(int n, int m);
cmat *cmat_init(int n, int m, cnum ***array);
int cmat_copy(cmat *m, cmat *dest);
int cmat_identity(int n, cmat *dest);
int cmat_zero(int n, int m, cmat *dest);

int cmat_add(cmat *a, cmat *b, cmat *dest);
int cmat_sub(cmat *a, cmat *b, cmat *dest);
int cmat_mul(cmat *a, cmat *b, cmat *dest);
int cmat_const_mul(cnum *a, cmat *c, cmat *dest);
int cmat_hermetian(cmat *src, cmat *dest);
float cmat_norm(cmat *m);
int cmat_normalize(cmat *m, cmat *dest);

// -------------------- Pool markers --------------------
uint16_t cnum_mark(void);
void cnum_reset(uint16_t mark);

uint16_t cnum_ptr_mark(void);
void cnum_ptr_reset(uint16_t mark);

uint16_t cnum_ptrptr_mark(void);
void cnum_ptrptr_reset(uint16_t mark);

uint16_t cmat_mark(void);
void cmat_reset(uint16_t mark);

// -------------------- Householder & QR --------------------
int householder_reflection_x(cmat *m, cmat *result);
int householder_embedding(cmat *p, int m, cmat *dest);
int QR_decomposition(cmat *A, cmat *R, cmat *Q);

// -------------------- Eigenvalues / Eigenvectors --------------------
int eigenvalues(cmat *A, cmat *dest);
int eigenvectors(cmat *A, cmat *dest);

// -------------------- Beamforming / MUSIC --------------------
int beamforming_array(int num_elements, float angle_rad, cmat *bf);
float music_algorithm_theta(int num_elements, cmat *Rxx, cmat *bf_theta);
float find_max_music(int num_elements, cmat *X, float start_angle_rad, float end_angle_rad, float step_rad);

// -------------------- Helpers --------------------
cnum *complex_exp(float r, float theta);
int find_largest_eigenvalue_index(cmat *eigvals);
int cmat_remove_lowest_eigenvector(int idx, cmat *eigenvects, cmat *dest);

// -------------------- Debug --------------------
void print_cnum(cnum *c, int newline);
void print_cmat(cmat *m);

#endif // CMATH_POOL_H
