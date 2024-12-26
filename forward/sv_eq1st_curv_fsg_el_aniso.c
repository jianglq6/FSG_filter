/*******************************************************************************
 * solver of isotropic elastic 1st-order eqn using 
 * fully staggered cuvgrids 
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "netcdf.h"

#include "blk_t.h"
#include "filter_fsg.h"
#include "sv_eq1st_curv_fsg_el_aniso.h"
#include "bdry_abs.h"

//#define M_NCERR(ierr) {fprintf(stderr,"sv_ nc error: %s\n", nc_strerror(ierr)); exit(1);}
#define Nt 100
#ifndef PI
#define PI 3.141592653589793
#endif


#ifndef M_NCERR
#define M_NCERR {fprintf(stderr,"sv_ nc error\n"); exit(1);}
#endif

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
  float dt, int nt_total, float t0)
{
  int ni = gdinfo->ni;
  int nk = gdinfo->nk;
  int lni1 = gdinfo->lni1;
  int lnk1 = gdinfo->lnk1;
  int lni2 = gdinfo->lni2;
  int lnk2 = gdinfo->lnk2;
  int lni = gdinfo->lni;
  int lnj = gdinfo->lnk;
  int lnx = gdinfo->lnx;
  int lny = gdinfo->lnz;

  int is_filter = fd->is_filter;

  size_t siz_line = gdinfo->siz_iz;
  size_t siz_line_half = gdinfo->siz_iz * 0.5;
  size_t siz_slice = gdinfo->siz_icmp;
  size_t siz_slice_qr = gdinfo->siz_icmp * 0.25;

  // fd and filter 
  int    fd_lebedev_half_len = fd->lebedev_half_len;
  float *fd_lebedev_coef = fd->lebedev_coef;
  int    fd_filter_half_len = fd->filter_coef_half_len;
  float *fd_filter_coef = fd->filter_coef;

  float *x2d_lebedev = gdcurv->x2d_lebedev;
  float *z2d_lebedev = gdcurv->z2d_lebedev;

  
  float *w2d = wav->v5d;

  // local pointer
  float *restrict w_cur;
  float *restrict w_rhs;
  float *restrict w_end;

  float t_cur;
  float t_end;  // time after this loop, for nc output
  size_t w3d_size_per_level =  siz_slice_qr * 5 * 2;
  // get wavefield
  w_cur = w2d + 0 * w3d_size_per_level;
  w_rhs = w2d + 1 * w3d_size_per_level;
  w_end = w2d + 2 * w3d_size_per_level;

  // boundary num
  int abs_num_of_layers[4] = {30, 30, 30, 30};
  for (int i=0; i<4; i++) {
      fprintf(stdout,"-> abs_num_of_layers[%d] = %d\n", i, abs_num_of_layers[i]);
  }

  for (int i=0; i<5; i++) {
      fprintf(stdout,"-> fd_filter_coef[%d] = %f\n", i, fd_filter_coef[i]); }
  printf("dx = %f, dz = %f\n", dx, dz);


  // snapshot create
  fprintf(stdout, "== Start create nc file\n");
  // for velocity
  int ncid;
  int timeid;
  int varid[2];
  int dimid[3];
  char fname[50];
  sprintf(fname, "%s/volume_vel.nc",output_dir);
  if (nc_create(fname, NC_CLOBBER | NC_64BIT_OFFSET, &ncid)) M_NCERR;
  if (nc_def_dim(ncid, "time", NC_UNLIMITED, &dimid[0])) M_NCERR;
  if (nc_def_dim(ncid, "k", nk, &dimid[1])) M_NCERR;
  if (nc_def_dim(ncid, "i", ni, &dimid[2])) M_NCERR;
  if (nc_def_var(ncid, "time", NC_FLOAT, 1, dimid+0, &timeid)) M_NCERR;
  if (nc_def_var(ncid, "Vx", NC_FLOAT, 3, dimid, &varid[0])) M_NCERR;
  if (nc_def_var(ncid, "Vz", NC_FLOAT, 3, dimid, &varid[1])) M_NCERR;
  if (nc_enddef(ncid)) M_NCERR;


  float *restrict  Vx2 = w_cur + 0 * siz_slice_qr;
  float *restrict  Vz2 = w_cur + 1 * siz_slice_qr;
  float *restrict Txx2 = w_cur + 2 * siz_slice_qr;
  float *restrict Tzz2 = w_cur + 3 * siz_slice_qr;
  float *restrict Txz2 = w_cur + 4 * siz_slice_qr;
  float *restrict  Vx1 = w_cur + 5 * siz_slice_qr;
  float *restrict  Vz1 = w_cur + 6 * siz_slice_qr;
  float *restrict Txx1 = w_cur + 7 * siz_slice_qr;
  float *restrict Tzz1 = w_cur + 8 * siz_slice_qr;
  float *restrict Txz1 = w_cur + 9 * siz_slice_qr;

  float *restrict  hVx2 = w_rhs + 0 * siz_slice_qr;
  float *restrict  hVz2 = w_rhs + 1 * siz_slice_qr;
  float *restrict hTxx2 = w_rhs + 2 * siz_slice_qr;
  float *restrict hTzz2 = w_rhs + 3 * siz_slice_qr;
  float *restrict hTxz2 = w_rhs + 4 * siz_slice_qr;
  float *restrict  hVx1 = w_rhs + 5 * siz_slice_qr;
  float *restrict  hVz1 = w_rhs + 6 * siz_slice_qr;
  float *restrict hTxx1 = w_rhs + 7 * siz_slice_qr;
  float *restrict hTzz1 = w_rhs + 8 * siz_slice_qr;
  float *restrict hTxz1 = w_rhs + 9 * siz_slice_qr;

  float *restrict  Vx_end = w_end + 0 * siz_slice;
  float *restrict  Vz_end = w_end + 1 * siz_slice;
  float *restrict Txx_end = w_end + 2 * siz_slice;
  float *restrict Tzz_end = w_end + 3 * siz_slice;
  float *restrict Txz_end = w_end + 4 * siz_slice;

  float *c1_11 = m2d->c11_lebedev + 0 * siz_slice_qr;
  float *c2_11 = m2d->c11_lebedev + 1 * siz_slice_qr;
  float *c3_11 = m2d->c11_lebedev + 2 * siz_slice_qr;
  float *c4_11 = m2d->c11_lebedev + 3 * siz_slice_qr;

  float *c1_13 = m2d->c13_lebedev + 0 * siz_slice_qr;
  float *c2_13 = m2d->c13_lebedev + 1 * siz_slice_qr;
  float *c3_13 = m2d->c13_lebedev + 2 * siz_slice_qr;
  float *c4_13 = m2d->c13_lebedev + 3 * siz_slice_qr;

  float *c1_15 = m2d->c15_lebedev + 0 * siz_slice_qr;
  float *c2_15 = m2d->c15_lebedev + 1 * siz_slice_qr;
  float *c3_15 = m2d->c15_lebedev + 2 * siz_slice_qr;
  float *c4_15 = m2d->c15_lebedev + 3 * siz_slice_qr;

  float *c1_33 = m2d->c33_lebedev + 0 * siz_slice_qr;
  float *c2_33 = m2d->c33_lebedev + 1 * siz_slice_qr;
  float *c3_33 = m2d->c33_lebedev + 2 * siz_slice_qr;
  float *c4_33 = m2d->c33_lebedev + 3 * siz_slice_qr;

  float *c1_35 = m2d->c35_lebedev + 0 * siz_slice_qr;
  float *c2_35 = m2d->c35_lebedev + 1 * siz_slice_qr;
  float *c3_35 = m2d->c35_lebedev + 2 * siz_slice_qr;
  float *c4_35 = m2d->c35_lebedev + 3 * siz_slice_qr;

  float *c1_55 = m2d->c55_lebedev + 0 * siz_slice_qr;
  float *c2_55 = m2d->c55_lebedev + 1 * siz_slice_qr;
  float *c3_55 = m2d->c55_lebedev + 2 * siz_slice_qr;
  float *c4_55 = m2d->c55_lebedev + 3 * siz_slice_qr;

  float *slw2d1 = m2d->rho_lebedev + 0  * siz_slice_qr;
  float *slw2d2 = m2d->rho_lebedev + 1  * siz_slice_qr;
  float *slw2d3 = m2d->rho_lebedev + 2  * siz_slice_qr;
  float *slw2d4 = m2d->rho_lebedev + 3  * siz_slice_qr;

  //
  // time loop
  //
  fprintf(stdout, "Start time loop\n");

  for (int it=0; it<nt_total; it++)
  {

    // 
    // update Vi
    //

    t_cur = it * dt + t0;
    t_end = t_cur + 0.5*dt;

    fprintf(stdout,"------> it=%d, t=\%f\n", it, t_cur);
    sv_eliso1st_curv_lebedev_cal_momentum_intv( 
                  hVx1, hVz1, Txx1, Tzz1, Txz1,
                  hVx2, hVz2, Txx2, Tzz2, Txz2,
                  slw2d1, slw2d4,
                  lni1, lni2, lnk1, lnk2,
                  dx, dz, siz_line_half,
                  fd_lebedev_half_len, fd_lebedev_coef);


    src_force_sinc_intv(
        hVx1, hVx2, hVz1, hVz2,
        slw2d1, slw2d4,
        lni1, lni2, lnk1, lnk2, 
        siz_line_half, siz_slice_qr,
        fd_lebedev_half_len, fd_lebedev_coef,
        x2d_lebedev, z2d_lebedev,
        dx, dz, t_cur, src);
               

    sv_eliso1st_curv_lebedev_update_velocity(
                   hVx1, hVz1, Vx1, Vz1,
                   hVx2, hVz2, Vx2, Vz2,
                   dt,  lnx, lny,
                   siz_line_half);


    abs_exp_velocity( lni1, lni2, lnk1, lnk2,
                      siz_line_half, abs_num_of_layers,
                      Vx1, Vz1, Vx2, Vz2);


    if (is_filter == 1) {
      filter_velocity_intv(
            Vx1, Vz1, Vx2, Vz2,
            lni1, lni2, lnk1, lnk2, siz_line_half,
            fd_filter_half_len, fd_filter_coef);
    }


    //-- print snapshot: velocity 
    size_t startp[] = { it, 0, 0 };
    size_t countp[] = { 1, nk, ni};
    size_t start_tdim = it;
    int snap_line = ni;
    // pack the wavefield into one array accroding to the variable layout
    for (int k = lnk1; k <= lnk2; k++) {
      for (int i = lni1; i <= lni2; i++) {
        int iptr = i + k * siz_line_half;
        // integer-integer grid point
        int iptr_00 = 2*(i-lni1) + 2 * (k-lnk1) * snap_line;
        Vx_end[iptr_00] = Vx1[iptr];
        Vz_end[iptr_00] = Vz2[iptr];

        // half-half grid point
        int iptr_11 = iptr_00 + snap_line + 1;
        Vx_end[iptr_11] = Vx2[iptr];
        Vz_end[iptr_11] = Vz1[iptr];

        // half-integer grid points: no value, 
        //  compute it by averaging surrounding values
//        int iptr_10 = iptr_00 + 1;
//        Vx_end[iptr_10] = (Vx1[iptr] + Vx1[iptr+1] 
//                         + Vx2[iptr] + Vx2[iptr-siz_line_half]) * 0.25;
//        Vz_end[iptr_10] = (Vz2[iptr] + Vz2[iptr+1] 
//                         + Vz1[iptr] + Vz1[iptr-siz_line_half]) * 0.25;
//
//        int iptr_01 = iptr_00 + snap_line;
//        Vx_end[iptr_01] = (Vx1[iptr] + Vx1[iptr+siz_line_half]
//                          +Vx2[iptr] + Vx2[iptr-1]) * 0.25;
//        Vz_end[iptr_01] = (Vz2[iptr] + Vz2[iptr+siz_line_half]
//                          +Vz1[iptr] + Vz1[iptr-1]) * 0.25;
//
      }
    }

    nc_put_var1_float(ncid,timeid,&start_tdim,&t_end);
    nc_put_vara_float(ncid,varid[0],startp,countp,Vx_end);
    nc_put_vara_float(ncid,varid[1],startp,countp,Vz_end);


    //
    //-- updata Tij
    //

    t_cur = t_end;
    t_end = t_cur + 0.5*dt;

    sv_eliso1st_curv_lebedev_aniso_cal_hook_intv(
      Vx1, Vz1, Vx2, Vz2,
      hTxx1, hTzz1, hTxz1,
      hTxx2, hTzz2, hTxz2,
      lni1, lni2, lnk1, lnk2,
      dx, dz,
      siz_line_half, siz_slice_qr,
      fd_lebedev_half_len, fd_lebedev_coef, 
      c2_11, c2_13, c2_15,
      c2_33, c2_35, c2_55,
      c3_11, c3_13, c3_15,
      c3_33, c3_35, c3_55);


    sv_eliso1st_curv_lebedev_update_stress( 
                   hTxx1, hTzz1, hTxz1,
                   Txx1,  Tzz1,  Txz1,
                   hTxx2, hTzz2, hTxz2,
                   Txx2,  Tzz2,  Txz2,
                   dt, lnx, lny,
                   siz_line_half);


    abs_exp_stress( lni1, lni2, lnk1, lnk2,
                    siz_line_half, abs_num_of_layers,
                    Txx1, Tzz1, Txz1, Txx2, Tzz1, Txz2);


    if (is_filter == 1) {
      filter_stress_intv(
            Txx1, Tzz1, Txz1, 
            Txx2, Tzz2, Txz2,
            lni1, lni2, lnk1, lnk2, siz_line_half,
            fd_filter_half_len, fd_filter_coef);
    }

    
  }


  if (nc_close(ncid)) M_NCERR;
//  if (nc_close(ncid_s)) M_NCERR;

  fprintf(stdout,"Successfully done!\n");
}


// -----------------------------------------------------------------------------//
//            vel is on integer grid     
// -----------------------------------------------------------------------------//
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
int fd_lebedev_half_len, float *fd_lebedev_coef)
{
  float fdx_coef[fd_lebedev_half_len];
  float fdz_coef[fd_lebedev_half_len];

  float slw1, slw4;
  float xix1, xiz1, ztx1, ztz1;
  float xix4, xiz4, ztx4, ztz4;

  float DxTxx1, DxTzz1, DxTxz1, DzTxx1, DzTxz1, DzTzz1;
  float DxTxx4, DxTzz4, DxTxz4, DzTxx4, DzTxz4, DzTzz4;

  for (int i=0; i < fd_lebedev_half_len; i++) {
    fdx_coef[i] = fd_lebedev_coef[i+fd_lebedev_half_len];
    fdz_coef[i] = fd_lebedev_coef[i+fd_lebedev_half_len];
  }


  for (size_t k=lnk1; k<=lnk2; k++)
  {
    size_t iptr_k = k * siz_line_half;

    size_t iptr = iptr_k + lni1;

    for (size_t i=lni1; i<=lni2; i++)      
    {

        DxTxx1 = 0.0; DxTzz1 = 0.0;  DxTxz1 = 0.0; DzTxx1 = 0.0; DzTzz1 = 0.0; DzTxz1 = 0.0;  
        DxTxx4 = 0.0; DxTzz4 = 0.0;  DxTxz4 = 0.0; DzTxx4 = 0.0; DzTzz4 = 0.0; DzTxz4 = 0.0;  

        for (int n=1; n <= fd_lebedev_half_len; n++) {
          // NO.1
          DxTxx1 += Dx0(Txx1, fdx_coef, n, i, k, siz_line_half);
          DxTzz1 += Dx0(Tzz1, fdx_coef, n, i, k, siz_line_half);
          DxTxz1 += Dx0(Txz2, fdx_coef, n, i, k, siz_line_half);
  
          DzTxx1 += Dz0(Txx2, fdz_coef, n, i, k, siz_line_half);
          DzTzz1 += Dz0(Tzz2, fdz_coef, n, i, k, siz_line_half);
          DzTxz1 += Dz0(Txz1, fdz_coef, n, i, k, siz_line_half);
  
          // NO.4
          DxTxx4 += Dx1(Txx2, fdx_coef, n, i, k, siz_line_half);
          DxTzz4 += Dx1(Tzz2, fdx_coef, n, i, k, siz_line_half);
          DxTxz4 += Dx1(Txz1, fdx_coef, n, i, k, siz_line_half);
  
          DzTxx4 += Dz1(Txx1, fdz_coef, n, i, k, siz_line_half);
          DzTzz4 += Dz1(Tzz1, fdz_coef, n, i, k, siz_line_half);
          DzTxz4 += Dz1(Txz2, fdz_coef, n, i, k, siz_line_half);
        }


        // medium
        slw1 = slw2d1[iptr];
        slw4 = slw2d4[iptr];


        // cart grid
        DxTxx1 = DxTxx1/dx;
        DxTxz1 = DxTxz1/dx;
        DzTxz1 = DzTxz1/dz;
        DzTzz1 = DzTzz1/dz;
    
        DxTxx4 = DxTxx4/dx;
        DxTxz4 = DxTxz4/dx;
        DzTxz4 = DzTxz4/dz;
        DzTzz4 = DzTzz4/dz;
        // the first staggered grid
        hVx1[iptr] = slw1*( DxTxx1 + DzTxz1 );
        hVz1[iptr] = slw4*( DxTxz4 + DzTzz4 );
    
        // the second staggered grid
        hVx2[iptr] = slw4*( DxTxx4 + DzTxz4 );
        hVz2[iptr] = slw1*( DxTxz1 + DzTzz1 );

        iptr += 1;


    } 
  } 
}


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
                                              float *restrict c3_55)
{

  float fdx_coef[fd_lebedev_half_len];
  float fdz_coef[fd_lebedev_half_len];


  float DxVx2, DzVx2, DxVz2, DzVz2;
  float DxVx3, DzVx3, DxVz3, DzVz3;

  float xix2, xiz2, ztx2, ztz2;
  float xix3, xiz3, ztx3, ztz3;

  float c11_2, c13_2, c15_2, c33_2, c35_2, c55_2;
  float c11_3, c13_3, c15_3, c33_3, c35_3, c55_3;


  for (int i=0; i < fd_lebedev_half_len; i++) {
    fdx_coef[i] = fd_lebedev_coef[i+fd_lebedev_half_len];
    fdz_coef[i] = fd_lebedev_coef[i+fd_lebedev_half_len];
  }


  for (size_t k=lnk1; k<=lnk2; k++)
  {
    size_t iptr_k = k * siz_line_half;

    size_t iptr = iptr_k + lni1;

    for (size_t i=lni1; i<=lni2; i++)      
    {

      DxVx2=0.0; DzVx2=0.0; DxVz2=0.0; DzVz2=0.0;
      DxVx3=0.0; DzVx3=0.0; DxVz3=0.0; DzVz3=0.0;

      for (int n=1; n <= fd_lebedev_half_len; n++) {
        // NO.2
        DxVx2 += Dx1(Vx1, fdx_coef, n, i, k, siz_line_half);
        DxVz2 += Dx1(Vz2, fdx_coef, n, i, k, siz_line_half);
        DzVx2 += Dz0(Vx2, fdz_coef, n, i, k, siz_line_half);
        DzVz2 += Dz0(Vz1, fdz_coef, n, i, k, siz_line_half);

        //NO.3
        DxVx3 += Dx0(Vx2, fdx_coef, n, i, k, siz_line_half);
        DxVz3 += Dx0(Vz1, fdx_coef, n, i, k, siz_line_half);
        DzVx3 += Dz1(Vx1, fdz_coef, n, i, k, siz_line_half);
        DzVz3 += Dz1(Vz2, fdz_coef, n, i, k, siz_line_half);
      }        

  
      c11_2 = c2_11[iptr]; 
      c13_2 = c2_13[iptr]; 
      c15_2 = c2_15[iptr]; 
      c33_2 = c2_33[iptr]; 
      c35_2 = c2_35[iptr]; 
      c55_2 = c2_55[iptr]; 
      
      c11_3 = c3_11[iptr]; 
      c13_3 = c3_13[iptr]; 
      c15_3 = c3_15[iptr]; 
      c33_3 = c3_33[iptr]; 
      c35_3 = c3_35[iptr]; 
      c55_3 = c3_55[iptr]; 

      // cart grid
      DxVx2 = DxVx2/dx;
      DzVx2 = DzVx2/dz;
      DzVz2 = DzVz2/dz;
      DxVz2 = DxVz2/dx;
  
      DxVx3 = DxVx3/dx;
      DzVx3 = DzVx3/dz;
      DzVz3 = DzVz3/dz;
      DxVz3 = DxVz3/dx;
  
  
      // first
      // cI_ij I denotes which point
      hTxx1[iptr] = c11_2 * DxVx2 + c13_2 * DzVz2 + c15_2 * (DxVz2 + DzVx2);
      hTzz1[iptr] = c13_2 * DxVx2 + c33_2 * DzVz2 + c35_2 * (DxVz2 + DzVx2);
      hTxz1[iptr] = c15_3 * DxVx3 + c35_3 * DzVz3 + c55_3 * (DxVz3 + DzVx3);
      
      // second
      // cI_ij I denotes which point
      hTxx2[iptr] = c11_3 * DxVx3 + c13_3 * DzVz3 + c15_3 * (DxVz3 + DzVx3);
      hTzz2[iptr] = c13_3 * DxVx3 + c33_3 * DzVz3 + c35_3 * (DxVz3 + DzVx3);
      hTxz2[iptr] = c15_2 * DxVx2 + c35_2 * DzVz2 + c55_2 * (DxVz2 + DzVx2);
    
      
      iptr += 1;

    }
  }
}



void
sv_eliso1st_curv_lebedev_update_velocity(
float *restrict hVx1, float *restrict hVz1,
float *restrict Vx1,  float *restrict Vz1,
float *restrict hVx2, float *restrict hVz2,
float *restrict Vx2,  float *restrict Vz2,
float dt, int lnx, int lnz,
size_t siz_line_half)
{
  for (size_t k=0; k < lnz; k++)
  {
    size_t iptr_k = k * siz_line_half;
    size_t iptr = iptr_k;

      for (size_t i=0; i < lnx; i++)      
      {
        Vx1[iptr] += dt * hVx1[iptr];
        Vz1[iptr] += dt * hVz1[iptr];

        Vx2[iptr] += dt * hVx2[iptr];
        Vz2[iptr] += dt * hVz2[iptr];

        iptr += 1;
      }
  }
}



void
sv_eliso1st_curv_lebedev_update_stress(
float *restrict hTxx1, float *restrict hTzz1, float *restrict hTxz1,
float *restrict Txx1,  float *restrict Tzz1,  float *restrict Txz1,
float *restrict hTxx2, float *restrict hTzz2, float *restrict hTxz2,
float *restrict Txx2,  float *restrict Tzz2,  float *restrict Txz2,
float dt, int lnx, int lny,
size_t siz_line_half)
{
  for (size_t k=0; k<lny; k++)
  {
    size_t iptr_k = k * siz_line_half;
    size_t iptr = iptr_k;

    for (size_t i=0; i<lnx; i++)
    {
        Txx1[iptr] += dt * hTxx1[iptr];
        Tzz1[iptr] += dt * hTzz1[iptr];
        Txz1[iptr] += dt * hTxz1[iptr];

        Txx2[iptr] += dt * hTxx2[iptr];
        Tzz2[iptr] += dt * hTzz2[iptr];
        Txz2[iptr] += dt * hTxz2[iptr];
    
        
        iptr += 1;

    }
  }
}






