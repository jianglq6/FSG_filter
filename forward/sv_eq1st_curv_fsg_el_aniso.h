#ifndef SV_EQ1ST_CURV_LEBEDEV_EL_ANISO_H
#define SV_EQ1ST_CURV_LEBEDEV_EL_ANISO_H

#include "fd_t.h"
#include "gd_info.h"
#include "gd_t.h"
#include "md_t.h"
#include "wav_t.h"

// -----------------------------------------------------------------------------//
//            vel is on integer grid     
// -----------------------------------------------------------------------------//
void
sv_eliso1st_curv_lebedev_allstep_intv(
	fd_t     *fd,
	gdinfo_t *gdinfo,
	gd_t     *gdcurv,
	md_t     *m2d,
  src_t    *src,
	wav_t    *wav,
	char *output_dir,
	float dx, float dz,
	float dt, int nt_total, float t0);

void
sv_eliso1st_curv_lebedev_aniso_cal_hook_intv(
float *restrict Vx1, float *restrict Vz1,
float *restrict Vx2, float *restrict Vz2, 
float *restrict hTxx1, float *restrict hTzz1, float *restrict hTxz1,
float *restrict hTxx2, float *restrict hTzz2, float *restrict hTxz2,
int lni1, int lni2, int lnk1,int lnk2,
float dx, float dz,
size_t siz_line_half, size_t siz_slice_qr,
int fd_lebedev_half_len, float *fd_lebedev_coef,
float *restrict c2_11, float *restrict c2_13, float *restrict c2_15, 
                       float *restrict c2_33, float *restrict c2_35,
                                              float *restrict c2_55,
float *restrict c3_11, float *restrict c3_13, float *restrict c3_15, 
                       float *restrict c3_33, float *restrict c3_35,
                                              float *restrict c3_55);

void
sv_eliso1st_curv_lebedev_cal_momentum_intv(
float *restrict hVx1, float *restrict hVz1,
float *restrict Txx1, float *restrict Tzz1, float *restrict Txz1,
float *restrict hVx2, float *restrict hVz2,
float *restrict Txx2, float *restrict Tzz2, float *restrict Txz2,
float *restrict slw2d1, float *restrict slw2d4,
int lni1, int lni2, int lnk1,int lnk2,
float dx, float dz,
size_t siz_line_half,
int fd_lebedev_half_len, float *fd_lebedev_coef);



void
sv_eliso1st_curv_lebedev_update_velocity(
float *restrict hVx1, float *restrict hVz1,
float *restrict Vx1,  float *restrict Vz1,
float *restrict hVx2, float *restrict hVz2,
float *restrict Vx2,  float *restrict Vz2,
float dt, int lnx, int lnz,
size_t siz_line_half);



void
sv_eliso1st_curv_lebedev_update_stress(
float *restrict hTxx1, float *restrict hTzz1, float *restrict hTxz1,
float *restrict Txx1,  float *restrict Tzz1,  float *restrict Txz1,
float *restrict hTxx2, float *restrict hTzz2, float *restrict hTxz2,
float *restrict Txx2,  float *restrict Tzz2,  float *restrict Txz2,
float dt, int lnx, int lny,
size_t siz_line_half);



#endif
