#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "fdlib_math.h"
#include "fdlib_mem.h"
#include "fd_t.h"
#include "bdry_abs.h"


// exp coef //
#define ALPHA 0.015

#ifndef PI
#define PI 3.1415926535898
#endif

//-----------------------------
// lebedev schem exp abs//
//-----------------------------

void
abs_exp_velocity
(int ni1, int ni2, int nk1, int nk2,
 size_t siz_line, int *boundary_layer_number,
 float *restrict Vx1, float *restrict Vz1,
 float *restrict Vx2, float *restrict Vz2)
{
    if( boundary_layer_number[0] > 0 ) {
        abs_exp_x1(ni1, ni2, nk1, nk2, siz_line,
                   boundary_layer_number[0], Vx1);

        abs_exp_x1(ni1, ni2, nk1, nk2, siz_line,
                   boundary_layer_number[0], Vx2);


        abs_exp_x1(ni1, ni2, nk1, nk2, siz_line, 
                   boundary_layer_number[0], Vz1);

        abs_exp_x1(ni1, ni2, nk1, nk2, siz_line, 
                   boundary_layer_number[0], Vz2);

    }

    if( boundary_layer_number[1] > 0 ) {
        abs_exp_x2(ni1, ni2, nk1, nk2, siz_line,
                   boundary_layer_number[1], Vx1);

        abs_exp_x2(ni1, ni2, nk1, nk2, siz_line,
                   boundary_layer_number[1], Vx2);

        abs_exp_x2(ni1, ni2, nk1, nk2, siz_line,
                   boundary_layer_number[1], Vz1);

        abs_exp_x2(ni1, ni2, nk1, nk2, siz_line,
                   boundary_layer_number[1], Vz2);

    }

    if( boundary_layer_number[2] > 0 ) {
        abs_exp_z1(ni1, ni2, nk1, nk2, siz_line,
                   boundary_layer_number[2], Vx1);

        abs_exp_z1(ni1, ni2, nk1, nk2, siz_line,
                   boundary_layer_number[2], Vx2);


        abs_exp_z1(ni1, ni2, nk1, nk2, siz_line,
                   boundary_layer_number[2], Vz1);

        abs_exp_z1(ni1, ni2, nk1, nk2, siz_line,
                   boundary_layer_number[2], Vz2);

    }


    if( boundary_layer_number[3] > 0 ) {
        abs_exp_z2(ni1, ni2, nk1, nk2, siz_line,
                   boundary_layer_number[3], Vx1);

        abs_exp_z2(ni1, ni2, nk1, nk2, siz_line,
                   boundary_layer_number[3], Vx2);


        abs_exp_z2(ni1, ni2, nk1, nk2, siz_line,
                   boundary_layer_number[3], Vz1);

        abs_exp_z2(ni1, ni2, nk1, nk2, siz_line,
                   boundary_layer_number[3], Vz2);
    }

}


void
abs_exp_stress
(int ni1, int ni2, int nk1, int nk2,
 size_t siz_line, int *boundary_layer_number,
 float *restrict Txx1,  float *restrict Tzz1, float *restrict Txz1,  
 float *restrict Txx2,  float *restrict Tzz2, float *restrict Txz2)
{
    if( boundary_layer_number[0] > 0 ) {
        abs_exp_x1(ni1, ni2, nk1, nk2, siz_line,
                   boundary_layer_number[0], Txx1);

        abs_exp_x1(ni1, ni2, nk1, nk2, siz_line,
                   boundary_layer_number[0], Tzz1);

        abs_exp_x1(ni1, ni2, nk1, nk2, siz_line,
                   boundary_layer_number[0], Txz1);

        abs_exp_x1(ni1, ni2, nk1, nk2, siz_line,
                   boundary_layer_number[0], Txx2);

        abs_exp_x1(ni1, ni2, nk1, nk2, siz_line,
                   boundary_layer_number[0], Tzz2);

        abs_exp_x1(ni1, ni2, nk1, nk2, siz_line,
                   boundary_layer_number[0], Txz2);

    }


    if( boundary_layer_number[1] > 0 ) {
        abs_exp_x2(ni1, ni2, nk1, nk2, siz_line,
                   boundary_layer_number[1], Txx1);
                   
        abs_exp_x2(ni1, ni2, nk1, nk2, siz_line,
                   boundary_layer_number[1], Tzz1);

        abs_exp_x2(ni1, ni2, nk1, nk2, siz_line,
                   boundary_layer_number[1], Txz1);
                   
        abs_exp_x2(ni1, ni2, nk1, nk2, siz_line,
                   boundary_layer_number[1], Txx2);

        abs_exp_x2(ni1, ni2, nk1, nk2, siz_line,
                   boundary_layer_number[1], Tzz2);

        abs_exp_x2(ni1, ni2, nk1, nk2, siz_line,
                   boundary_layer_number[1], Txz2);
    }


    if( boundary_layer_number[2] > 0 ) {
        abs_exp_z1(ni1, ni2, nk1, nk2, siz_line,
                   boundary_layer_number[2], Txx1);

        abs_exp_z1(ni1, ni2, nk1, nk2, siz_line,
                   boundary_layer_number[2], Tzz1);

        abs_exp_z1(ni1, ni2, nk1, nk2, siz_line,
                   boundary_layer_number[2], Txz1);

        abs_exp_z1(ni1, ni2, nk1, nk2, siz_line,
                   boundary_layer_number[2], Txx2);

        abs_exp_z1(ni1, ni2, nk1, nk2, siz_line,
                   boundary_layer_number[2], Tzz2);

        abs_exp_z1(ni1, ni2, nk1, nk2, siz_line,
                   boundary_layer_number[2], Txz2);
    }

    if( boundary_layer_number[3] > 0 ) {
        abs_exp_z2(ni1, ni2, nk1, nk2, siz_line,
                   boundary_layer_number[3], Txx1);

        abs_exp_z2(ni1, ni2, nk1, nk2, siz_line,
                   boundary_layer_number[3], Tzz1);

        abs_exp_z2(ni1, ni2, nk1, nk2, siz_line,
                   boundary_layer_number[3], Txz1);

        abs_exp_z2(ni1, ni2, nk1, nk2, siz_line,
                   boundary_layer_number[3], Txx2);

        abs_exp_z2(ni1, ni2, nk1, nk2, siz_line,
                   boundary_layer_number[3], Tzz2);

        abs_exp_z2(ni1, ni2, nk1, nk2, siz_line,
                   boundary_layer_number[3], Txz2);
    }

}


float
cal_exp(int N)
{
    float factor;
    factor = exp( - ALPHA*N *ALPHA*N );
    return factor;
}



void
abs_exp_x1
(int ni1, int ni2, int nk1, int nk2,
 size_t siz_line, int boundary_layer_number, float *u)
{
    int bni1, bni2, bnk1, bnk2;
    float damp = 1.0;
    bni1 = ni1;
    bni2 = ni1 - 1 + boundary_layer_number;
    bnk1 = nk1;
    bnk2 = nk2;

    for (int ix = bni1; ix <bni2; ix++) {
      for (int iz = bnk1; iz < bnk2; iz++) {
        size_t iptr = iz*siz_line+ix;
        damp = cal_exp( boundary_layer_number-(ix-bni1));
        u[iptr] *= damp;
      }
    }
}


void
abs_exp_z1
(int ni1, int ni2, int nk1, int nk2,
 size_t siz_line, int boundary_layer_number, float *u)
{
    int bni1, bni2, bnk1, bnk2;
    float damp = 1.0;
    bni1 = ni1;
    bni2 = ni2;
    bnk1 = nk1;
    bnk2 = nk1 - 1 + boundary_layer_number;

    for (int ix = bni1; ix < bni2; ix++) {
      for (int iz = bnk1; iz < bnk2; iz++) {
        size_t iptr = iz*siz_line+ix;
        damp = cal_exp( boundary_layer_number-(iz-bnk1));
        u[iptr] *= damp;
      }
        
    }
}



void
abs_exp_x2
(int ni1, int ni2, int nk1, int nk2,
 size_t siz_line, int boundary_layer_number, float *u)
{
    int bni1, bni2, bnk1, bnk2;
    float damp = 1.0;
    bni1 = ni2 - boundary_layer_number;
    bni2 = ni2 - 1;
    bnk1 = nk1;
    bnk2 = nk2;

    for (int ix = bni1; ix <bni2; ix++) {
      for (int iz = bnk1; iz < bnk2; iz++) {
        size_t iptr = iz*siz_line+ix;
        damp = cal_exp( boundary_layer_number-(bni2-ix));
        u[iptr] *= damp;
      }
    }
}


void
abs_exp_z2
(int ni1, int ni2, int nk1, int nk2,
 size_t siz_line, int boundary_layer_number, float *u)
{
    int bni1, bni2, bnk1, bnk2;
    float damp = 1.0;
    bni1 = ni1;
    bni2 = ni2;
    bnk1 = nk2 - boundary_layer_number;
    bnk2 = nk2-1;

    for (int ix = bni1; ix <bni2; ix++) {
      for (int iz = bnk1; iz < bnk2; iz++) {
        size_t iptr = iz*siz_line+ix;
        damp = cal_exp( boundary_layer_number-(bnk2-iz));
        u[iptr] *= damp;
      }
    }
}






//
// abl exp type
//
/*
int abs_set_ablexp(size_t nx, size_t ny, size_t nz, 
    size_t ni1, size_t ni2, size_t nj1, size_t nj2, size_t nk1, size_t nk2, 
    int *boundary_itype, // input
    int *in_abs_numbers, //
    float *abs_alpha, //
    float *abs_beta, //
    float *abs_velocity, //
    int *abs_numbers, // output
    size_t *abs_coefs_dimpos, 
    float **p_abs_coefs)
{
  int ivar;
  
  float *Ax, *Bx, *Dx;
  float *Ay, *By, *Dy;
  float *Az, *Bz, *Dz;

  int num_of_coefs = 1; // damping
  
  size_t abs_ceofs_size = 0;
  
  // copy input to struct
  memcpy(abs_numbers, in_abs_numbers, FD_NDIM_2 * sizeof(int));

  // size
  for (i=0; i<2; i++) { // x1,x2
    abs_coefs_dimpos[i] = abs_ceofs_size;
    abs_ceofs_size += abs_numbers[i] * nj * nk; 
  }
  for (i=2; i<4; i++) { // y1,y2
    abs_coefs_dimpos[i] = abs_ceofs_size;
    abs_ceofs_size += abs_numbers[i] * ni * nk; 
  }
  for (i=4; i<6; i++) { // z1,z2
    abs_coefs_dimpos[i] = abs_ceofs_size;
    abs_ceofs_size += abs_numbers[i] * ni * nj; 
  }

  *p_abs_coef_size = abs_coefs_size;

  // vars
  *p_abs_coefs = (float *) fdlib_mem_calloc_1d_float( 
               abs_coefs_size, 0.0, "abs_set_ablexp");
  if (*p_abs_coefs == NULL) {
      fprintf(stderr,"Error: failed to alloc ablexp coefs\n");
      fflush(stderr);
      ierr = -1;
  }

  // set damping values

  return ierr;
}

int abs_ablexp_cal_damp(i,Vs,ah,nb)
{
  int ierr = 0;

  integer,intent(in) :: i,nb
  real(SP),intent(in) :: Vs,ah
  real(SP) :: d
  real(SP) :: ie
  integer m,n
  ie=i
  !Vs=5000.0_SP
  m=(nb*ah)/(Vs*stept)
  d=0.0_SP
  do n=1,m
     d=d+(n*stept*Vs)**2/(nb*ah)**2
  end do
  d=0.8_SP/d*1.1_SP
  d=exp(-d*(ie/nb)**2)

  return ierr;
}
*/
