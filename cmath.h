#include <stdint.h>

typedef struct cnum_implementation cnum; 
typedef struct cmat_implementation cmat;

cnum *cnum_init(float a, float b);
cmat *cmat_init(float n, float m, cnum ***cnum_arr);
cnum *cadd(cnum *c1, cnum *c2);
cnum *mul(cnum *c1, cnum *c2);
cnum *cnum_conj(cnum *c); 
void print_cnum(cnum *c, int newline);
void print_cmat(cmat *m);
float cnum_mag(cnum *c1);
void delete_cmat(cmat *m);
cmat *cmat_mul(cmat *a, cmat *b);
cnum ***cnum_array_init(float n, float m);
cmat *cmat_add(cmat *a, cmat *b);
cmat *cmat_sub(cmat *a, cmat *b);
cmat *cmat_identity(float n);
cmat *cmat_transpose(cmat *m);
cmat *cmat_hermetian(cmat *m);
cnum *complex_exp(float r, float theta);
cnum *csub(cnum *c1, cnum *c2);
cmat* householder_reflection_x(cmat *m);
cmat *cmat_const_mul(cnum *a, cmat *c);
float cmat_norm(cmat *m);
cmat *cmat_normalize(cmat *m);
cmat *houseolder_embedding(cmat *p, int m);
cmat *QR_decomposition(cmat *A, cmat **R_out);
cmat *eigenvalues(cmat *A);
cmat *eigenvectors(cmat *A);
cmat *cmat_copy(cmat *m);
float cnum_mag_squared(cnum *c1);
int find_lowest_eigenvalue_index(cmat *eigvals);
cmat *cmat_remove_lowest_eigenvector(int idx, cmat *eigenvects);
float music_algorithm_theta(int num_elements, cmat *X, cmat *bf_theta);
float find_max_music(int num_elements, cmat *X, float start_angle_rad, float end_angle_rad, float step_rad);
void delete_cnum(cnum *c);