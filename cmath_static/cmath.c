// cmath_pool.c
// Pool-based implementation with TEMP (mark/reset) and PERM (persistent) pools.
// No malloc/free. Designed for embedded (nRF52 etc).

#include "cmath.h"
#include <stdio.h>
#include <stdint.h>
#include <math.h>



#define M_PI   3.14159265358979323846
#define MAX_ITERS 5

// -----------------------------
// Configuration: pool sizes
// -----------------------------


#define CNUM_PERM_POOL_SIZE   512
#define CNUM_TEMP_POOL_SIZE   5000
#define PTR_PERM_POOL_SIZE    512// pointer slots (for cnum*/cnum** arrays)
#define PTR_TEMP_POOL_SIZE    5000
#define CMAT_PERM_POOL_SIZE   512
#define CMAT_TEMP_POOL_SIZE   5000



struct cnum_implementation { float a; float b; };
struct cmat_implementation {
    int n;
    int m;
    cnum ***cmat;   // pointer to pointer-to-pointer layout
};

// Use typedefs from header
// typedef struct cnum_implementation cnum;
// typedef struct cmat_implementation cmat;

// -----------------------------
// Permanent pools (never reset)
// -----------------------------
static cnum   cnum_perm_pool[CNUM_PERM_POOL_SIZE];
static uint16_t cnum_perm_index = 0;

static uintptr_t ptr_perm_pool[PTR_PERM_POOL_SIZE]; // store pointers generically
static uint16_t   ptr_perm_index = 0;

static cmat   cmat_perm_pool[CMAT_PERM_POOL_SIZE];
static uint16_t cmat_perm_index = 0;

// -----------------------------
// Temporary pools (stack-style, mark/reset)
// -----------------------------
static cnum   cnum_temp_pool[CNUM_TEMP_POOL_SIZE];
static uint16_t cnum_temp_index = 0;

static uintptr_t ptr_temp_pool[PTR_TEMP_POOL_SIZE];
static uint16_t   ptr_temp_index = 0;

static cmat   cmat_temp_pool[CMAT_TEMP_POOL_SIZE];
static uint16_t cmat_temp_index = 0;

// -----------------------------
// Mark/reset helpers for temp pools
// -----------------------------
static inline uint16_t cnum_temp_mark(void) { return cnum_temp_index; }
static inline void cnum_temp_reset(uint16_t mark) { cnum_temp_index = mark; }

static inline uint16_t ptr_temp_mark(void) { return ptr_temp_index; }
static inline void ptr_temp_reset(uint16_t mark) { ptr_temp_index = mark; }

static inline uint16_t cmat_temp_mark(void) { return cmat_temp_index; }
static inline void cmat_temp_reset(uint16_t mark) { cmat_temp_index = mark; }

// -----------------------------
// Basic allocators
// -----------------------------
// Permanent allocators (for objects that must survive across calls)
static inline cnum* cnum_perm_alloc(void) {
    if (cnum_perm_index >= CNUM_PERM_POOL_SIZE) {
		printf("cnum_perm_alloc: out of memory\n");
		return NULL;
    }
    return &cnum_perm_pool[cnum_perm_index++];
}
static inline void* ptr_perm_alloc(int count) {
    if (ptr_perm_index + count >= PTR_PERM_POOL_SIZE) { 
		printf("ptr_perm_alloc: out of memory\n");
		return NULL;
    }	
    void* block = (void*)&ptr_perm_pool[ptr_perm_index];
    ptr_perm_index += count;
    return block;
}
static inline cmat* cmat_perm_alloc(void) {
    if (cmat_perm_index >= CMAT_PERM_POOL_SIZE) { 
		printf("cmat_perm_alloc: out of memory\n");
		return NULL;
	}
    // zero-initialize pointers for safety
    cmat *m = &cmat_perm_pool[cmat_perm_index++];
    for (int i=0;i< (int)sizeof(m->cmat)/sizeof(void*); ++i) { /*noop*/ }
    return m;
}

// Temporary allocators (for short lived temporaries)
static inline cnum* cnum_temp_alloc(void) {
    if (cnum_temp_index >= CNUM_TEMP_POOL_SIZE) { 
		printf("cnum_temp_alloc: out of memory\n");
		return NULL;
	}
    return &cnum_temp_pool[cnum_temp_index++];
}
static inline void* ptr_temp_alloc(int count) {
    if (ptr_temp_index + count >= PTR_TEMP_POOL_SIZE) { 
		printf("ptr_temp_alloc: out of memory\n");
		return NULL;
	}
    void* block = (void*)&ptr_temp_pool[ptr_temp_index];
    ptr_temp_index += count;
    return block;
}
static inline cmat* cmat_temp_alloc(void) {
    if (cmat_temp_index >= CMAT_TEMP_POOL_SIZE) { 
		printf("cmat_temp_alloc: out of memory\n");
		return NULL;
	}
    return &cmat_temp_pool[cmat_temp_index++];
}

// -----------------------------
// Public "API-compatible" constructors
// By default these create TEMP objects (most internal temporaries).
// Use `promote_*` to copy them into PERM pools if needed.
// -----------------------------
cnum *cnum_init(float a, float b) {
    cnum *c = cnum_temp_alloc();
    if (!c) return NULL;
    c->a = a; c->b = b;
    return c;
}

// allocate pointer-arrays and cnum entries as TEMP
cnum ***cmat_array_init(int n, int m) {
    if (n <= 0 || m <= 0) return NULL;
    // allocate top-level pointer array (n entries) from ptr_temp_pool
    cnum ***mat = (cnum ***) ptr_temp_alloc(n);
    if (!mat) return NULL;
    for (int i = 0; i < n; ++i) {
        mat[i] = (cnum **) ptr_temp_alloc(m);
        if (!mat[i]) return NULL;
        for (int j = 0; j < m; ++j) {
            mat[i][j] = cnum_temp_alloc();
            if (!mat[i][j]) return NULL;
            mat[i][j]->a = 0.0f;
            mat[i][j]->b = 0.0f;
        }
    }
    return mat;
}

// cmat_init creates a TEMP cmat (caller may promote to PERM)
cmat *cmat_init(int n, int m, cnum ***array) {
    cmat *mat = cmat_temp_alloc();
    if (!mat) return NULL;
    mat->n = n;
    mat->m = m;
    mat->cmat = array;
    return mat;
}

// -----------------------------
// Promote helpers - deep copy TEMP object into PERM pools
// Must be used when you want to return an object that survives reset.
// -----------------------------
static cnum *cnum_promote_copy(cnum *src) {
    cnum *dst = cnum_perm_alloc();
    if (!dst) return NULL;
    dst->a = src->a;
    dst->b = src->b;
    return dst;
}

cnum ***cmat_array_promote_copy(cnum ***src, int n, int m) {
    cnum ***dst = (cnum ***) ptr_perm_alloc(n);
    if (!dst) return NULL;
    for (int i = 0; i < n; ++i) {
        dst[i] = (cnum **) ptr_perm_alloc(m);
        if (!dst[i]) return NULL;
        for (int j = 0; j < m; ++j) {
            dst[i][j] = cnum_promote_copy(src[i][j]);
            if (!dst[i][j]) return NULL;
        }
    }
    return dst;
}

cmat *cmat_promote_copy(cmat *src) {
    if (!src) return NULL;
    cnum ***array_copy = cmat_array_promote_copy(src->cmat, src->n, src->m);
    if (!array_copy) return NULL;
    cmat *dst = cmat_perm_alloc();
    if (!dst) return NULL;
    dst->n = src->n;
    dst->m = src->m;
    dst->cmat = array_copy;
    return dst;
}



// -----------------------------
// Utility mark/reset functions for callers (expose if needed)
// -----------------------------
uint16_t cnum_mark(void) { return cnum_temp_mark(); }
void cnum_reset(uint16_t mark) { cnum_temp_reset(mark); }

uint16_t ptr_mark(void) { return ptr_temp_mark(); }
void ptr_reset(uint16_t mark) { ptr_temp_reset(mark); }

uint16_t cmat_mark(void) { return cmat_temp_mark(); }
void cmat_reset(uint16_t mark) { cmat_temp_reset(mark); }

// -----------------------------
// Math helpers (use TEMP allocations internally)
// -----------------------------

// define complex exponential and magnitude functions
cnum *complex_exp(float r, float theta) {
    return cnum_init(r * cosf(theta), r * sinf(theta));
}

float cnum_mag(cnum *c1) {
    return sqrtf(c1->a * c1->a + c1->b * c1->b);
}
float cnum_mag_squared(cnum *c1) {
    return c1->a*c1->a + c1->b*c1->b;
}


// change these to be inplace operations to reduce allocations
int cadd(cnum *c1, cnum *c2, cnum *dest) { 
	if (!dest) { 
		return -1; 
	}
	float tempa = c1->a + c2->a;
	float tempb = c1->b + c2->b;
	dest->a = tempa;
	dest->b = tempb;
	return 0;
}

int csub(cnum *c1, cnum *c2, cnum *dest) { 
	if (!dest) { 
		return -1; 
	}
	float tempa = c1->a - c2->a;
	float tempb = c1->b - c2->b;
	dest->a = tempa;
	dest->b = tempb;
	return 0; 
}

int cmul(cnum *c1, cnum *c2, cnum *dest) { 
	if (!dest) { 
		return -1; 
	}
	float tempa = c1->a*c2->a - c1->b*c2->b;
	float tempb = c1->a*c2->b + c1->b*c2->a;
	dest->a = tempa;
	dest->b = tempb;
	return 0;
}

int cnum_conj(cnum *c, cnum *dest) { 
	if (!dest) { 
		return -1; 
	}
	float tempa = c->a;
	float tempb = -c->b;
	dest->a = tempa;
	dest->b = tempb;
	return 0;
}

// -----------------------------
// Matrix operations all done in place
// -----------------------------
int cmat_copy(cmat *src, cmat *dest) {
	if (!dest) { 
		return -1;
	}
    dest->n = src->n;
    dest->m = src->m;

    

    for (int i = 0; i < src->n; i++) {
        for (int j = 0; j < src->m; j++) {
            dest->cmat[i][j]->a = src->cmat[i][j]->a;
            dest->cmat[i][j]->b = src->cmat[i][j]->b;
        }
    }
	return 0;
}

// turns destination into an identity matrix
int cmat_identity(int n, cmat *dest) {

    if (!dest || dest->n != n || dest->m != n) {
		printf("cmat_identity: dimension mismatch or null destination\n");
		return -1;
    }

    for (int i=0;i<n;i++) {
		for (int j=0;j<n;j++) {
			if (i == j) {
        		dest->cmat[i][j]->a = 1.0;
        		dest->cmat[i][j]->b = 0.0;
			} else {
				dest->cmat[i][j]->a = 0.0;
				dest->cmat[i][j]->b = 0.0;
			}
		}
	}
    return 0;
}

// note A and B cannot be the same
int cmat_mul(cmat *a, cmat *b, cmat *dest) {
    if (a->m != b->n || !dest) {
        printf("Dimension mismatch or null destination\n");
        return -1;
    }

    cnum temp;
    for (int i = 0; i < a->n; i++) {
        for (int j = 0; j < b->m; j++) {
            dest->cmat[i][j]->a = 0;
            dest->cmat[i][j]->b = 0;
            for (int k = 0; k < a->m; k++) {
                cmul(a->cmat[i][k], b->cmat[k][j], &temp);
                cadd(dest->cmat[i][j], &temp, dest->cmat[i][j]);
            }
        }
    }
    return 0;
}


// multiply by constant
int cmat_const_mul(cnum *a, cmat *c, cmat *dest) {

    for (int i = 0; i < c->n; i++) {
        for (int j = 0; j < c->m; j++) {

            // read source
            float ar = a->a;
            float ai = a->b;
            float br = c->cmat[i][j]->a;
            float bi = c->cmat[i][j]->b;

            // write result into *existing* cnum in dest
            dest->cmat[i][j]->a = ar*br - ai*bi;
            dest->cmat[i][j]->b = ar*bi + ai*br;
        }
    }

    return 0;
}

float cmat_norm(cmat *v) {
    if (v->m != 1) { printf("Input must be column vector\n"); return -1; }
    float s=0;
    for (int i=0;i<v->n;i++) s += cnum_mag_squared(v->cmat[i][0]);
    return sqrtf(s);
}

int cmat_normalize(cmat *v, cmat *dest) {
    float n = cmat_norm(v);
    if (n == 0.0) {
		cmat_copy(v, dest);
		return 0;
	} 

	// generate an extra cnum here
    cnum *tmp = cnum_init(1.0/n, 0.0);
	cmat_const_mul(tmp, v, dest);
    return 0;
}

int cmat_add(cmat *a, cmat *b, cmat *dest) {
    if (a->n!=b->n||a->m!=b->m) { 
		printf("Dimension mismatch\n"); 
		return -1; 
	}
	for (int i=0;i<a->n;i++) { 
		for (int j=0;j<a->m;j++) {
			cadd(a->cmat[i][j], b->cmat[i][j], dest->cmat[i][j]);
		}
	}
	return 0;
}
int cmat_sub(cmat *a, cmat *b, cmat *dest) {
    if (a->n!=b->n||a->m!=b->m) { 
		printf("Dimension mismatch\n"); 
		return -1; 
	}
	for (int i=0;i<a->n;i++) { 
		for (int j=0;j<a->m;j++) {
			csub(a->cmat[i][j], b->cmat[i][j], dest->cmat[i][j]);
		}
	}
	return 0;
}



int cmat_hermetian(cmat *src, cmat *dest) {
	if (!dest) { 
		return -1;
	}
    dest->n = src->m;
    dest->m = src->n;

    //dest->cmat = cmat_array_init(src->m, src->n);

    for (int i = 0; i < src->n; i++) {
        for (int j = 0; j < src->m; j++) {
			float tempa = src->cmat[i][j]->a;
			float tempb = -src->cmat[i][j]->b;
            dest->cmat[j][i]->a = tempa;
            dest->cmat[j][i]->b = tempb;
        }
    }
	return 0;
}


// -----------------------------
// Printing
// -----------------------------
void print_cnum(cnum *c, int newline) {
    if (newline) printf("%f + %fj\n", c->a, c->b);
    else {
		printf("%f + %fj ", c->a, c->b);
	}
}
void print_cmat(cmat *m) {
    for (int i=0;i<m->n;i++){
        for (int j=0;j<m->m;j++) print_cnum(m->cmat[i][j], 0);
        printf("\n");
    }
}



// -----------------------------
// Householder, QR, Eigen, MUSIC
// All functions allocate temporaries in TEMP pools.
// Before returning matrices that must survive, we PROMOTE copies into PERM pools.
// -----------------------------

// requires about 5 matrices of allocation
int householder_reflection_x(cmat *m, cmat *result) {
	if (!result) { 
		printf("householder_reflection_x: null destination\n");
		return -1;
	}

	// I: identity matrix
	cnum *** I_array = cmat_array_init(m->n, m->n);
	cmat* I = cmat_init(m->n, m->n, I_array);
	cmat_identity(m->n, I);
	//print_cmat(I);

	// b: unit vector
    cnum ***b_arr = cmat_array_init(m->n, 1);
    cmat *b = cmat_init(m->n, 1, b_arr);

    // set first element

    //cnum *norm_m = cnum_init(cmat_norm(m), 0);
    float mag_0 = cnum_mag(m->cmat[0][0]);
	float mag_m = cmat_norm(m);

	//printf("m->cmat[0][0]: %f + %fj\n", m->cmat[0][0]->a, m->cmat[0][0]->b);
    if (m->cmat[0][0]->b == 0) {

		
		if (m->cmat[0][0]->a >= 0) {
			b->cmat[0][0]->a = mag_m;
			b->cmat[0][0]->b = 0;
		} else {
			b->cmat[0][0]->a = -mag_m;
			b->cmat[0][0]->b = 0;
		}
		
		
    } else {
		//printf("test comp\n");
		b->cmat[0][0]->a = mag_m * (m->cmat[0][0]->a / mag_0);
		b->cmat[0][0]->b = mag_m * (m->cmat[0][0]->b / mag_0);

    }
	//printf("Vector b:\n");
	//print_cmat(b);
	cnum ***u_array = cmat_array_init(m->n, 1);
	cnum ***u_temp_array = cmat_array_init(m->n, 1);
	cnum ***u_t_array = cmat_array_init(1, m->n);
	cnum ***proj_array = cmat_array_init(m->n, m->n);
	cnum ***temp_array = cmat_array_init(m->n, m->n);

	cmat *u = cmat_init(m->n, 1, u_array);	
	cmat *u_temp = cmat_init(m->n, 1, u_temp_array);
	cmat *u_t = cmat_init(1, m->n, u_t_array);	
	cmat *proj  = cmat_init(m->n, m->n, proj_array);
	cmat *temp  = cmat_init(m->n, m->n, temp_array);;

	
	

    cmat_sub(m, b, u_temp);
    cmat_normalize(u_temp, u);
	//printf("Vector u:\n");
	//print_cmat(u);
	
    cmat_hermetian(u, u_t);
	cmat_mul(u, u_t, proj);

	//printf("Vector b:\n");
	//print_cmat(b);	
	//printf("Proj matrix:\n");
	//print_cmat(proj);
	


	
    cnum *two = cnum_init(2,0);
    cmat_const_mul(two, proj, temp);
    cmat_sub(I, temp, result);

    // promote temp result into permanent pools (so caller can use it)
    //cmat *result = cmat_promote_copy(result_temp);

    // restore temps (frees temporaries)
	
	return 0;
   
}


int householder_embedding(cmat *p, int m, cmat *dest) {
    // produce PERM result because likely used outside - create temp then promote
   // uint16_t m1 = cnum_mark(); uint16_t p1 = ptr_mark(); uint16_t q1 = cmat_mark();
	if (!dest) { 
		return -1;
	}	

	cmat_identity(m, dest);

    for (int i=0;i<p->n;i++) { 
		for (int j=0;j<p->n;j++) {
			dest->cmat[m - p->n + i][m - p->n + j] = p->cmat[i][j];
		}
	}
    //cmat *perm = cmat_promote_copy(I);
    //cnum_reset(m1); ptr_reset(p1); cmat_reset(q1);
    return 0;
}


int QR_decomposition(cmat *A, cmat *R, cmat *Q) {


    int n = A->n;
    cmat_copy(A, R);    
	//print_cmat(R);        // temp
    cmat_identity(n, Q);  // temp
	//print_cmat(Q);        // temp

	cnum ***emb_array = cmat_array_init(n, n);
	cmat *E_k = cmat_init(n, n, emb_array);

	cnum ***tmp_R_array = cmat_array_init(n, n);
	cmat *tmp_R = cmat_init(n, n, tmp_R_array);

	cnum ***tmp_Q_array = cmat_array_init(n, n);
	cmat *tmp_Q = cmat_init(n, n, tmp_Q_array);


    for (int k=0;k<n-1;k++) {
		//printf("Iteration k=%d\n", k);
		cnum ***householder_result_array = cmat_array_init(n-k, n-k);
		cmat *P_k = cmat_init(n-k, n-k, householder_result_array);
        cnum ***subvec_arr = cmat_array_init(R->n - k, 1);
		//printf("R matrix at start of iteration %d:\n", k);
		//print_cmat(R);
        for (int i=k; i<R->n; i++) {
			subvec_arr[i-k][0]->a = R->cmat[i][k]->a;
			subvec_arr[i-k][0]->b = R->cmat[i][k]->b;
		}

			cmat *subvec = cmat_init(R->n - k, 1, subvec_arr);
			//printf("Subvector:\n");
			//print_cmat(subvec);

			householder_reflection_x(subvec, P_k); 
			//printf("Householder P_k:\n");
			//print_cmat(P_k);
			householder_embedding(P_k, n, E_k);  
			//printf("Embedded E_k:\n");
			//print_cmat(E_k);  


			cmat_copy(R, tmp_R);
			cmat_copy(Q, tmp_Q);

			/*
			printf("Tmp R before multiply:\n");
			print_cmat(tmp_R);	
			printf("Tmp Q before multiply:\n");
			print_cmat(tmp_Q);
			*/
			
			cmat_mul(E_k, tmp_R, R);
			cmat_mul(tmp_Q, E_k, Q);
			/*
			printf("R after E_k * R:\n");
			print_cmat(R);	
			printf("Q after E_k * Q:\n");	
			print_cmat(Q);
			*/

    }

    // zero lower triangle 
	for (int i=0;i<R->n;i++) {
		for (int j=0;j<i;j++) {
			R->cmat[i][j]->a = 0.0;
			R->cmat[i][j]->b = 0.0;
		}
	}
	



   
    return 0;
}



int eigenvalues(cmat *A, cmat *dest) {
	if (!dest) { 
		return -1;
	}

	// pass in reduntant arrays between functions
    cmat_copy(A, dest);  
	
	cnum ***Q_array = cmat_array_init(A->n, A->n);
	cmat *Q = cmat_init(A->n, A->n, Q_array);
	cnum***R_array = cmat_array_init(A->n, A->n);
	cmat *R = cmat_init(A->n, A->n, R_array);

    for (int i=0;i<MAX_ITERS;i++) {
        QR_decomposition(dest, R, Q);
        cmat_mul(R, Q, dest);
		
	
        // promote A_next to permanent to keep for next iter (or copy into A0)
        //cmat *A_next_perm = cmat_promote_copy(A_next);
        // note: temporaries R, Q, A_next are temp and will be discarded by QR internals
    }

    return 0;
}


int eigenvectors(cmat *A, cmat *dest) {
	if (!dest) { 
		return -1;
	}

	cnum ***A_temp_array = cmat_array_init(A->n, A->n);
	cmat *A_temp = cmat_init(A->n, A->n, A_temp_array);
    cmat_copy(A, A_temp);

	cnum ***Q_array = cmat_array_init(A->n, A->n);
	cmat *Q = cmat_init(A->n, A->n, Q_array);
	cnum***R_array = cmat_array_init(A->n, A->n);
	cmat *R = cmat_init(A->n, A->n, R_array);

	cnum ***E0_array = cmat_array_init(A->n, A->n);
	cmat *E0 = cmat_init(A->n, A->n, E0_array);
	cmat_identity(A->n, E0);

    for (int i=0;i<MAX_ITERS;i++) {
        QR_decomposition(A_temp, R, Q);
        cmat_mul(R, Q, A_temp);

		/*
		printf("Iteration %d:\n", i);
		printf("R:\n");
		print_cmat(R);		
		printf("Q:\n");
		print_cmat(Q);
		*/

		cmat_mul(E0, Q, dest);
		cmat_copy(dest, E0);

    }
	
    return 0;
}


int find_largest_eigenvalue_index(cmat *eigvals) {
    int idx = 0;
    float max_mag = cnum_mag(eigvals->cmat[0][0]);
    for (int j=1;j<eigvals->m;j++) {
        float mag = cnum_mag(eigvals->cmat[j][j]);
        if (mag > max_mag) { max_mag = mag; idx = j; }
    }
    return idx;
}


int cmat_remove_lowest_eigenvector(int idx, cmat *eigenvects, cmat *dest) {
	if (!dest) { 
		return -1;
	}
	int n = eigenvects->n;
	int m = eigenvects->m;
	int new_m = m - 1;

	for (int i=0;i<n;i++) {
		int dest_col = 0;
		for (int j=0;j<m;j++) {
			if (j == idx) continue;
			dest->cmat[i][dest_col]->a = eigenvects->cmat[i][j]->a;
			dest->cmat[i][dest_col]->b = eigenvects->cmat[i][j]->b;
			dest_col++;
		}
	}
	return 0;
}

int beamforming_array(int num_elements, float angle_rad, cmat *bf) {
	if (!bf) { 
		return -1;
	}

    float d = 0.5;
    for (int n=0;n < num_elements;n++) { 
		bf->cmat[n][0]->a = cosf(2.0*M_PI * d * n * sinf(angle_rad));
		bf->cmat[n][0]->b = sinf(2.0*M_PI * d * n * sinf(angle_rad));
	}
    return 0;
}


float music_algorithm_theta(int num_elements, cmat *Rxx, cmat *bf_theta) {
	cnum ***bfH_array = cmat_array_init(1, num_elements);
	cmat *bfH = cmat_init(1, num_elements, bfH_array);
	cmat_hermetian(bf_theta, bfH);

	cnum ***eigvals_array = cmat_array_init(Rxx->n, Rxx->n);
	cmat *eigvals = cmat_init(Rxx->n, Rxx->n, eigvals_array);
	cnum ***eigv_array = cmat_array_init(Rxx->n, Rxx->n);
	cmat *eigv = cmat_init(Rxx->n, Rxx->n, eigv_array);

	cnum ***reduced_eigv_array = cmat_array_init(Rxx->n, Rxx->n - 1);
	cmat *reduced_eigv = cmat_init(Rxx->n, Rxx->n - 1, reduced_eigv_array);

	cnum ***reduced_eigv_H_array = cmat_array_init(Rxx->n-1, Rxx->n);
	cmat *reduced_eigv_H = cmat_init(Rxx->n, Rxx->n - 1, reduced_eigv_H_array);	

	cnum ***prod1_array = cmat_array_init(Rxx->n-1 , 1);
	cmat *prod1 = cmat_init(Rxx->n-1, 1, prod1_array);

	cnum ***prod2_array = cmat_array_init(Rxx->n , 1);
	cmat *prod2 = cmat_init(Rxx->n, 1, prod2_array);

	cnum ***denom_array = cmat_array_init(1, 1);
	cmat *denom = cmat_init(1, 1, denom_array);

    eigenvalues(Rxx, eigvals);
    eigenvectors(Rxx, eigv);

	/*
	printf("Eigenvalue matrix:\n");
	print_cmat(eigvals);
	printf("Eigenvector matrix:\n");
	print_cmat(eigv);
	*/

    int largest_idx = find_largest_eigenvalue_index(eigvals);
    cmat_remove_lowest_eigenvector(largest_idx, eigv, reduced_eigv);
    cmat_hermetian(reduced_eigv, reduced_eigv_H);

	/*
	printf("Reduced eigenvector matrix:\n");
	print_cmat(reduced_eigv);	
	printf("Reduced eigenvector Hermetian matrix:\n");
	print_cmat(reduced_eigv_H);
	*/


    cmat_mul(reduced_eigv_H, bf_theta, prod1);
	//printf("Product 1:\n");
	//print_cmat(prod1);
    cmat_mul(reduced_eigv, prod1, prod2);
	//printf("Product 2:\n");
	//print_cmat(prod2);
    cmat_mul(bfH, prod2, denom);
	//printf("bfH :\n");
	//print_cmat(bfH);
	//printf("Denominator:\n");
	//print_cmat(denom);

    float denom_mag = cnum_mag(denom->cmat[0][0]);
    float P = 1.0f / denom_mag;
    return P;
}


float find_max_music(int num_elements, cmat *X, float start_angle_rad, float end_angle_rad, float step_rad) {
    float maxP = -1.0f;
    float best = start_angle_rad;
	cnum ***bf_array = cmat_array_init(num_elements, 1);
	cmat *bf = cmat_init(num_elements, 1, bf_array);

	printf("Scanning angles \n");
    for (float th = start_angle_rad; th <= end_angle_rad; th += step_rad) {
		printf("iteration at angle %f rad\n", th);
		beamforming_array(num_elements, th, bf);
		printf("Beamforming array:\n");
		print_cmat(bf);
        float p = music_algorithm_theta(num_elements, X, bf);
		printf("P = %f\n", p);
        if (p > maxP) { maxP = p; best = th; }
        // free temporaries by resetting temp pools between iterations:
		
        cnum_temp_index = 0;
        ptr_temp_index = 0;
        cmat_temp_index = 0;
		cnum_perm_index = 0;
		ptr_perm_index = 0;
		cmat_perm_index = 0;
		

    }
    return best * 180.0f / M_PI;
}




