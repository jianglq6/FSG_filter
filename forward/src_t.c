// test for lebedev
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "src_t.h"
#include "gd_t.h" 
#include "gd_info.h"

#ifndef PI
#define PI 3.141592653589793
#endif

  
// -----------------------------------------------------------------------------//
//            vel is on integer grid for LG    
// -----------------------------------------------------------------------------//

void src_init_bypar(src_t *src, int *src_indx, float *force_vector,
    float source_peak_time, float source_freq)
{
  src->si = src_indx[0];
  src->sk = src_indx[1];

  src->Fx = force_vector[0];
  src->Fz = force_vector[1];

  src->stf_timefactor = source_peak_time;
  src->stf_freqfactor = source_freq;

}

void src_print(src_t *src, gdinfo_t *gdinfo, gd_t *gdcurv)
{
  int si = src->si * 2 + gdinfo->ni1;
  int sk = src->sk * 2 + gdinfo->nk1;

  size_t indx = si + gdinfo->siz_iz * sk;

  float coordx = gdcurv->x2d[indx];
  float coordz = gdcurv->z2d[indx];
  fprintf(stdout, "  source coord: (%g, %g) m.\n", coordx, coordz);

}


// force src for vel is integer
// ref: FD-consistent point source (Koene et al., 2020, 2021)
void src_force_sinc_intv(
		float *hVx1, float *hVx2, float *hVz1, float *hVz2,
    float *restrict slw2d1, float *restrict slw2d4,
		int lni1, int lni2, int lnj1, int lnj2, 
		size_t siz_line_half, size_t siz_slice_qr,
		int fd_lebedev_half_len, float *fd_lebedev_coef,
		float *x2d_lebedev, float *z2d_lebedev,
		float dx, float dz, float current_time,
    src_t *src) 
{
  float *restrict x2d1  = x2d_lebedev + 0 * siz_slice_qr;
  float *restrict x2d2  = x2d_lebedev + 1 * siz_slice_qr;
  float *restrict x2d3  = x2d_lebedev + 2 * siz_slice_qr;
  float *restrict x2d4  = x2d_lebedev + 3 * siz_slice_qr;

  float *restrict y2d1  = z2d_lebedev + 0 * siz_slice_qr;
  float *restrict y2d2  = z2d_lebedev + 1 * siz_slice_qr;
  float *restrict y2d3  = z2d_lebedev + 2 * siz_slice_qr;
  float *restrict y2d4  = z2d_lebedev + 3 * siz_slice_qr;

	int si = src->si; 
  int sk = src->sk;
	float Fx = src->Fx;
  float Fz = src->Fz;
  int indx_xs, indx_zs, indx_source;
  float stf;
	float slw1, slw4;

  stf = cal_source_time_function(SRC_SIG_RICKER, current_time,
              src->stf_timefactor, src->stf_freqfactor);
  
  indx_xs = si + lni1;
  indx_zs = sk + lnj1;

  indx_source = indx_zs*siz_line_half+indx_xs;

  for (int j=indx_zs-30; j<=indx_zs+30; j++) { 
		for (int i=indx_xs-30; i<=indx_xs+30; i++) {

		  int iptr = i + j * siz_line_half;

		  float x1 = x2d1[iptr] - x2d1[indx_source];
		  float z1 = y2d1[iptr] - y2d1[indx_source];

		  float x4 = x2d4[iptr] - x2d1[indx_source];
		  float z4 = y2d4[iptr] - y2d1[indx_source];

		  float delta_x1 = cal_delta(x1, dx, fd_lebedev_half_len, fd_lebedev_coef);
		  float delta_z1 = cal_delta(z1, dz, fd_lebedev_half_len, fd_lebedev_coef);

		  float delta_x4 = cal_delta(x4, dx, fd_lebedev_half_len, fd_lebedev_coef);
		  float delta_z4 = cal_delta(z4, dz, fd_lebedev_half_len, fd_lebedev_coef);

      slw1 = slw2d1[iptr];
      slw4 = slw2d4[iptr];

		  hVx1[iptr] += slw1 * delta_x1 * delta_z1 * Fx * stf;
		  hVz2[iptr] += slw1 * delta_x1 * delta_z1 * Fz * stf;

		  hVx2[iptr] += slw4 * delta_x4 * delta_z4 * Fx * stf;
		  hVz1[iptr] += slw4 * delta_x4 * delta_z4 * Fz * stf;

	  }
	}
}



float cal_delta(float x, float dx,  int fd_lebedev_half_len, float *fd_lebedev_coef)
{
	float delta = 0.0;
	float fdx_coef[fd_lebedev_half_len];
	float sinc1, sinc2;
	float del;

  for (int i=0; i < fd_lebedev_half_len; i++) {

    fdx_coef[i] = fd_lebedev_coef[i+fd_lebedev_half_len];

	  sinc1 = sinc_func(x/dx-(i+1-0.5));
	  sinc2 = sinc_func(x/dx+(i+1-0.5));
	
	  del = fdx_coef[i] * (i+1-0.5)/dx * ( sinc1 + sinc2 );
	  delta += del;

	}
	
	return delta;
}


float sinc_func(float x){

	float sinc = sin(PI*x)/(PI*x);

	if (x == 0.0) sinc = 1.0;

	if (fabs(sinc) > 1e-15) {
		return sinc;
	} else {
		return 0.0;
	}
}


float cal_source_time_function(int flag_stf_type, float t, float t0, float f0)
{
    float stf;

    if (flag_stf_type == SRC_SIG_RICKER) {
        stf = fun_ricker(t,f0,t0);
    } else if (flag_stf_type == SRC_SIG_GAUSS) {
        stf = fun_gauss(t,f0,t0);
    }

    return stf;

}


float fun_ricker(float t, float fc, float t0)
{
    float u, f0,v;
    if(t<=0.0) {
        v = 0.0;
    } else {
        f0 = sqrt(PI)/2.0;
        u = (t-t0)*2.0*PI*fc;
        v = (u*u/4-0.5)*exp(-u*u/4)*f0;
    }

    return v;
}

float fun_gauss(float t, float a, float t0)
{
    float f;
    if(fabs(t0)>SRC_ZERO && (t<=0.0 || t>=2*t0)) {
        f = 0;
    } else {
        f = exp(-(t-t0)*(t-t0)/(a*a))/(sqrt(PI)*a);
    }

    return f;
}
