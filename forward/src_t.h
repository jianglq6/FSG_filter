#ifndef SRC_FUNCS_FSG_H
#define SRC_FUNCS_FSG_H

#include "constants.h"
#include "gd_info.h"
#include "gd_t.h"

#define SRC_IMPULSE_DIRECT 0
#define SRC_IMPULSE_AVERAGE 1

#define SRC_SIG_RICKER 0
#define SRC_SIG_GAUSS 1
#define SRC_ZERO 1e-20

typedef struct {
  int max_nt; // max nt of stf and mrf per src
  // for getting value in calculation
  int it;

  float stf_freqfactor;
  float stf_timefactor;

  // time independent
  int si; // local i index 
  int sk; // local k index 

  // time dependent
  // force stf
  float Fx; 
  float Fz;

  // moment rate
  float Mxx;
  float Mzz;
  float Mxz;

} src_t;


// -----------------------------------------------------------------------------//
//            vel is on integer grid for LG    
// -----------------------------------------------------------------------------//
  
void src_init_bypar(src_t *src, int *src_indx, float *force_vector,
    float source_peak_time, float source_freq);

void src_print(src_t *src, gdinfo_t *gdinfo, gd_t *gdcurv);

void src_force_sinc_intv(
		float *hVx1, float *hVx2, float *hVz1, float *hVz2,
    float *restrict slw2d1, float *restrict slw2d4,
		int lni1, int lni2, int lnj1, int lnj2, 
		size_t siz_line_half, size_t siz_slice_qr,
		int fd_lebedev_half_len, float *fd_lebedev_coef,
		float *x2d_lebedev, float *z2d_lebedev,
		float dx, float dz, 
		float current_time, src_t *src);

float 
cal_source_time_function(int flag_stf_type, float t, float t0, float f0);

float 
fun_ricker(float t, float fc, float t0);

float 
fun_gauss(float t, float a, float t0);


float cal_delta(float x, float dx,  int fd_lebedev_half_len, float *fd_lebedev_coef);

float sinc_func(float x);

#endif
