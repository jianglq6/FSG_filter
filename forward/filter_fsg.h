#ifndef FILTER_FSG_H
#define FILTER_FSG_H


// -----------------------------------------------------------------------------//
//            vel is on integer grid     
// -----------------------------------------------------------------------------//
void
filter_stress_intv(
float *restrict Txx1,  float *restrict Tzz1,  float *restrict Txz1,
float *restrict Txx2,  float *restrict Tzz2,  float *restrict Txz2,
int lni1, int lni2, int lnj1, int lnj2, size_t siz_line_half,
int half_len, float *coef);

void 
filter_velocity_intv(
float *restrict Vx1,  float *restrict Vz1,
float *restrict Vx2,  float *restrict Vz2,
int lni1, int lni2, int lnj1, int lnj2, size_t siz_line_half,
int half_len, float *coef);


// lower left
int
filter_lower_left135( float *u0, float *u1, size_t iptr, int half_len, float *coef,
		size_t dim2, size_t dim1);

int
filter_lower_left45( float *u0, float *u1, size_t iptr, int half_len, float *coef,
		size_t dim2, size_t dim1);



// upper left
int
filter_upper_left135( float *u0, float *u1, size_t iptr, int half_len, float *coef,
		size_t dim2, size_t dim1);

int
filter_upper_left45( float *u0, float *u1, size_t iptr, int half_len, float *coef,
		size_t dim2, size_t dim1);


// upper right
int
filter_upper_right135( float *u0, float *u1, size_t iptr, int half_len, float *coef,
		size_t dim2, size_t dim1);

int
filter_upper_right45( float *u0, float *u1, size_t iptr, int half_len, float *coef,
		size_t dim2, size_t dim1);


// lower right
int
filter_lower_right135( float *u0, float *u1, size_t iptr, int half_len, float *coef,
		size_t dim2, size_t dim1);

int
filter_lower_right45( float *u0, float *u1, size_t iptr, int half_len, float *coef,
		size_t dim2, size_t dim1);

#endif
