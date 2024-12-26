#ifndef ABS_FUNCS_FSG_H
#define ABS_FUNCS_FSG_H

// set exp for lebedev

void
abs_exp_velocity
(int ni1, int ni2, int nk1, int nk2,
 size_t siz_line, int *boundary_layer_number,
 float *restrict Vx1, float *restrict Vz1,
 float *restrict Vx2, float *restrict Vz2);

void
abs_exp_stress
(int ni1, int ni2, int nk1, int nk2,
 size_t siz_line,  int *boundary_layer_number,
 float *restrict Txx1,  float *restrict Tzz1, float *restrict Txz1,
 float *restrict Txx2,  float *restrict Tzz2, float *restrict Txz2);


float
cal_exp(int N);

void
abs_exp_x1
(int ni1, int ni2, int nk1, int nk2,
 size_t siz_line, int boundary_layer_number, float *u);


void
abs_exp_z1
(int ni1, int ni2, int nk1, int nk2,
 size_t siz_line, int boundary_layer_number, float *u);

void
abs_exp_x2
(int ni1, int ni2, int nk1, int nk2,
 size_t siz_line,  int boundary_layer_number, float *u);

void
abs_exp_z2
(int ni1, int ni2, int nk1, int nk2,
 size_t siz_line, int boundary_layer_number, float *u);

#endif
