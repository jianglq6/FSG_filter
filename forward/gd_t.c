/*
********************************************************************************
* Curve grid metric calculation                                                *
********************************************************************************
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "netcdf.h"

#include "fdlib_mem.h"
#include "fdlib_math.h"
#include "fd_t.h"
#include "gd_t.h"
#include "md_t.h"
//#include "isPointInHexahedron.h"

// used in read grid file
#define M_gd_INDEX( i, j, k, ni, nj ) ( ( i ) + ( j ) * ( ni ) + ( k ) * ( ni ) * ( nj ) )

void 
gd_curv_init_lebedev(gdinfo_t *gdinfo, gd_t *gdcurv)
{
  /*
   * 0-2: x3d, y3d, z3d
   */

  gdcurv->type = GD_TYPE_CURV;

  gdcurv->ncmp = CONST_NDIM * 2;

  gdcurv->nx   = gdinfo->nx;
  gdcurv->nz   = gdinfo->nz;

  gdcurv->siz_iz   = gdcurv->nx;
  gdcurv->siz_icmp = gdcurv->nx * gdcurv->nz;
  
  // vars
  gdcurv->v4d = (float *) fdlib_mem_calloc_1d_float(
                  gdcurv->siz_icmp * gdcurv->ncmp, 0.0, "gd_curv_init");
  if (gdcurv->v4d == NULL) {
      fprintf(stderr,"Error: failed to alloc coord vars\n");
      fflush(stderr);
  }
  
  // position of each v4d
  size_t *cmp_pos = (size_t *) fdlib_mem_calloc_1d_sizet(gdcurv->ncmp,
                                                         0,
                                                         "gd_curv_init");
  
  // set value
  int icmp = 0;
  cmp_pos[icmp] = icmp * gdcurv->siz_icmp;
  gdcurv->x2d = gdcurv->v4d + cmp_pos[icmp];

  icmp += 1;
  cmp_pos[icmp] = icmp * gdcurv->siz_icmp;
  gdcurv->z2d = gdcurv->v4d + cmp_pos[icmp];

  icmp += 1;
  cmp_pos[icmp] = icmp * gdcurv->siz_icmp;
  gdcurv->x2d_lebedev = gdcurv->v4d + cmp_pos[icmp];

  icmp += 1;
  cmp_pos[icmp] = icmp * gdcurv->siz_icmp;
  gdcurv->z2d_lebedev = gdcurv->v4d + cmp_pos[icmp];
  
  // set pointer
  gdcurv->cmp_pos  = cmp_pos;

  return;
}



/*
 * generate cartesian grid for curv struct
 */
void
gd_curv_gen_cart(
  gdinfo_t *gdinfo,
  gd_t *gdcurv,
  float dx, float x0_glob,
  float dz, float z0_glob,
  float *x1d, float *z1d)

{
  float *x2d = gdcurv->x2d;
  float *z2d = gdcurv->z2d;

  float x0 = x0_glob + (0 - gdinfo->fdx_nghosts) * dx;
  float z0 = z0_glob + (0 - gdinfo->fdz_nghosts) * dz;

  size_t iptr = 0;
  for (size_t k=0; k<gdcurv->nz; k++)
  {
	  z1d[k] = z0 + k * dz;
      for (size_t i=0; i<gdcurv->nx; i++)
      {
        x2d[iptr] = x0 + i * dx;
        z2d[iptr] = z0 + k * dz;

        iptr++;
      }
  }
  
  for (size_t i = 0; i<gdcurv->nx; i++) {
	  x1d[i] = x0 + i * dx;
  }

  return;
}


//
// input/output
//
void
gd_curv_coord_export(
  gdinfo_t *gdinfo,
  gd_t *gdcurv,
  char *output_dir)
{
  size_t *restrict c3d_pos   = gdcurv->cmp_pos;
  char  **restrict c3d_name  = gdcurv->cmp_name;
  int number_of_vars = gdcurv->ncmp;
  int  nx = gdcurv->nx;
  int  nz = gdcurv->nz;
  int  ni1 = gdinfo->ni1;
  int  nk1 = gdinfo->nk1;
  int  ni  = gdinfo->ni;
  int  nk  = gdinfo->nk;

  // construct file name
  char ou_file[CONST_MAX_STRLEN];
  sprintf(ou_file, "%s/coord.nc", output_dir);
  
  // read in nc
  int ncid;
  int varid[gdcurv->ncmp];
  int dimid[CONST_NDIM];

  int ierr = nc_create(ou_file, NC_CLOBBER, &ncid);
  if (ierr != NC_NOERR){
    fprintf(stderr,"creat coord nc error: %s\n", nc_strerror(ierr));
    exit(-1);
  }

  // define dimension
  ierr = nc_def_dim(ncid, "i", nx, &dimid[1]);
  ierr = nc_def_dim(ncid, "k", nz, &dimid[0]);

  // define vars
  for (int ivar=0; ivar<gdcurv->ncmp; ivar++) {
    ierr = nc_def_var(ncid, gdcurv->cmp_name[ivar], NC_FLOAT, CONST_NDIM, dimid, &varid[ivar]);
  }

  // attribute: index in output snapshot, index w ghost in thread
  int l_start[] = { ni1, nk1 };
  nc_put_att_int(ncid,NC_GLOBAL,"local_index_of_first_physical_points",
                   NC_INT,CONST_NDIM,l_start);

  int l_count[] = { ni, nk };
  nc_put_att_int(ncid,NC_GLOBAL,"count_of_physical_points",
                   NC_INT,CONST_NDIM,l_count);

  // end def
  ierr = nc_enddef(ncid);

  // add vars
  for (int ivar=0; ivar<gdcurv->ncmp; ivar++) {
    float *ptr = gdcurv->v4d + gdcurv->cmp_pos[ivar];
    ierr = nc_put_var_float(ncid, varid[ivar],ptr);
  }
  
  // close file
  ierr = nc_close(ncid);
  if (ierr != NC_NOERR){
    fprintf(stderr,"nc error: %s\n", nc_strerror(ierr));
    exit(-1);
  }

  return;
}

void
gd_curv_coord_import(gd_t *gdcurv, char *import_dir)
{
  // construct file name
  char in_file[CONST_MAX_STRLEN];
  sprintf(in_file, "%s/coord.nc", import_dir);
  
  // read in nc
  int ncid;
  int varid;

  int ierr = nc_open(in_file, NC_NOWRITE, &ncid);
  if (ierr != NC_NOERR){
    fprintf(stderr,"open coord nc error: %s\n", nc_strerror(ierr));
    exit(-1);
  }

  // read vars
  for (int ivar=0; ivar<gdcurv->ncmp; ivar++)
  {
    float *ptr = gdcurv->v4d + gdcurv->cmp_pos[ivar];

    ierr = nc_inq_varid(ncid, gdcurv->cmp_name[ivar], &varid);
    if (ierr != NC_NOERR){
      fprintf(stderr,"nc error: %s\n", nc_strerror(ierr));
      exit(-1);
    }

    ierr = nc_get_var(ncid, varid, ptr);
    if (ierr != NC_NOERR){
      fprintf(stderr,"nc error: %s\n", nc_strerror(ierr));
      exit(-1);
    }
  }
  
  // close file
  ierr = nc_close(ncid);
  if (ierr != NC_NOERR){
    fprintf(stderr,"nc error: %s\n", nc_strerror(ierr));
    exit(-1);
  }

  return;
}


void
gd_curv_set_minmax(gd_t *gdcurv)
{
  float xmin = gdcurv->x2d[0], xmax = gdcurv->x2d[0];
  float zmin = gdcurv->z2d[0], zmax = gdcurv->z2d[0];
  
  for (size_t i = 0; i < gdcurv->siz_icmp; i++){
      xmin = xmin < gdcurv->x2d[i] ? xmin : gdcurv->x2d[i];
      xmax = xmax > gdcurv->x2d[i] ? xmax : gdcurv->x2d[i];
      zmin = zmin < gdcurv->z2d[i] ? zmin : gdcurv->z2d[i];
      zmax = zmax > gdcurv->z2d[i] ? zmax : gdcurv->z2d[i];
  }

  gdcurv->xmin = xmin;
  gdcurv->xmax = xmax;
  gdcurv->zmin = zmin;
  gdcurv->zmax = zmax;

  return;
}

/*
 * convert cart coord to global index
 */

int
gd_cart_coord_to_local_indx(gdinfo_t *gdinfo,
                           gd_t *gdcart,
                           float sx,
                           float sz,
                           int   *ou_si, int *ou_sk,
                           float *ou_sx_inc, float *ou_sz_inc)
{
  int ierr = 0;

  int si_glob = (int)( (sx - gdcart->x0_glob) / gdcart->dx + 0.5 );
  int sk_glob = (int)( (sz - gdcart->z0_glob) / gdcart->dz + 0.5 );
  float sx_inc = si_glob * gdcart->dx + gdcart->x0_glob - sx;
  float sz_inc = sk_glob * gdcart->dz + gdcart->z0_glob - sz;

  *ou_si = si_glob + gdinfo->fdx_nghosts;
  *ou_sk = sk_glob + gdinfo->fdz_nghosts;
  *ou_sx_inc = sx_inc;
  *ou_sz_inc = sz_inc;

  return ierr; 
}

/* 
 * if the nearest point in this thread then search its grid index
 *   return value:
 *      1 - in this thread
 *      0 - not in this thread
 */

int
gd_curv_coord_to_local_indx(gdinfo_t *gdinfo,
                        gd_t *gd,
                        float sx, float sz,
                        int *si, int *sk,
                        float *sx_inc, float *sz_inc,
                        float *restrict wrk3d)
{
  int is_here = 0; // default outside

  int nx = gdinfo->nx;
  int nz = gdinfo->nz;
  int ni1 = gdinfo->ni1;
  int ni2 = gdinfo->ni2;
  int nk1 = gdinfo->nk1;
  int nk2 = gdinfo->nk2;
  size_t siz_iz = gdinfo->siz_iz;

  float *restrict x3d = gd->x2d;
  float *restrict z3d = gd->z2d;

  // outside coord range
  if ( sx < gd->xmin || sx > gd->xmax ||
       sz < gd->zmin || sz > gd->zmax)
  {
    is_here = 0;
    return is_here;
  }

  // init closest point
  float min_dist = sqrtf(  (sx - x3d[0]) * (sx - x3d[0])
      + (sz - z3d[0]) * (sz - z3d[0]) );
  int min_dist_i = 0 ;
  int min_dist_k = 0 ;

  // compute distance to each point
  for (int k=0; k<nz; k++) {
      for (int i=0; i<nx; i++)
      {
        size_t iptr = i + k * siz_iz;

        float x = x3d[iptr];
        float z = z3d[iptr];

        float DistInt = sqrtf(  (sx - x) * (sx - x)
            + (sz - z) * (sz - z) );
        wrk3d[iptr] =  DistInt;

        // replace closest point
        if (min_dist > DistInt)
        {
          min_dist = DistInt;
          min_dist_i = i;
          min_dist_k = k;
        }
      }
  }

  // if nearest index is outside phys region, not here
  if ( min_dist_i < ni1 || min_dist_i > ni2 ||
      min_dist_k < nk1 || min_dist_k > nk2 )
  {
    is_here = 0;
    return is_here;
  }

  // in this thread
  is_here = 1;

  float points_x[4];
  float points_z[4];
  float points_i[4];
  float points_k[4];

  for (int kk=0; kk<2; kk++)
  {
      for (int ii=0; ii<2; ii++)
      {
        int cur_i = min_dist_i-1+ii;
        int cur_k = min_dist_k-1+kk;

        for (int n3=0; n3<2; n3++) {
            for (int n1=0; n1<2; n1++) {
              int iptr_cube = n1 + n3 * 2;
              int iptr = (cur_i+n1)  + (cur_k+n3) * siz_iz;
              points_x[iptr_cube] = x3d[iptr];
              points_z[iptr_cube] = z3d[iptr];
              points_i[iptr_cube] = cur_i+n1;
              points_k[iptr_cube] = cur_k+n3;
            }
        }

        if (fdlib_math_isPoint2InQuad(sx,sz,points_x,points_z) == 1)
        {
          float si_curv, sk_curv;

          gd_curv_coord2index_sample(sx,sz,
              points_x,points_z,
              points_i,points_k,
              100,100,
              &si_curv, &sk_curv);

          // convert to return values
          *si = min_dist_i;
          *sk = min_dist_k;
          *sx_inc = si_curv - min_dist_i;
          *sz_inc = sk_curv - min_dist_k;

          return is_here;
        }
      }
  }

  // if not in any cube due to bug, set default value
  //    if everything is right, should be return 10 line before
  *si = min_dist_i;
  *sk = min_dist_k;
  *sx_inc = 0.0;
  *sz_inc = 0.0;

  return is_here;
}

/* 
 * find curv index using sampling
 */

  int
gd_curv_coord2index_sample(float sx, float sz, 
    float *points_x, // x coord of all points
    float *points_z,
    float *points_i, // curv coord of all points
    float *points_k,
    int    nx_sample,
    int    nz_sample,
    float *si_curv, // interped curv coord
    float *sk_curv)
{
  float Lx[2], Lz[2];

  // init closest point
  float min_dist = sqrtf(  (sx - points_x[0]) * (sx - points_x[0])
      + (sz - points_z[0]) * (sz - points_z[0]) );
  int min_dist_i = 0 ;
  int min_dist_k = 0 ;

  // linear interp for all sample
  for (int n3=0; n3<nz_sample+1; n3++)
  {
    Lz[1] = (float)(n3) / (float)(nz_sample);
    Lz[0] = 1.0 - Lz[1];
      for (int n1=0; n1<nx_sample+1; n1++)
      {
        Lx[1] = (float)(n1) / (float)(nx_sample);
        Lx[0] = 1.0 - Lx[1];

        // interp
        float x_pt=0;
        float z_pt=0;
        for (int kk=0; kk<2; kk++) {
            for (int ii=0; ii<2; ii++)
            {
              int iptr_cube = ii + kk * 2;
              x_pt += Lx[ii]*Lz[kk] * points_x[iptr_cube];
              z_pt += Lx[ii]*Lz[kk] * points_z[iptr_cube];
            }
        }

        // find min dist
        float DistInt = sqrtf(  (sx - x_pt) * (sx - x_pt)
            + (sz - z_pt) * (sz - z_pt) );

        // replace closest point
        if (min_dist > DistInt)
        {
          min_dist = DistInt;
          min_dist_i = n1;
          min_dist_k = n3;
        }
      } // n1
  } // n3

  *si_curv = points_i[0] + (float)min_dist_i / (float)nx_sample;
  *sk_curv = points_k[0] + (float)min_dist_k / (float)nz_sample;

  return 0;
}

/* 
 * interp curv coord using inverse distance interp
 */

int
gd_curv_coord2index_rdinterp(float sx, float sz, 
    int num_points,
    float *points_x, // x coord of all points
    float *points_z,
    float *points_i, // curv coord of all points
    float *points_k,
    float *si_curv, // interped curv coord
    float *sk_curv)
{
  float weight[num_points];
  float total_weight = 0.0 ;

  // cal weight
  int at_point_indx = -1;
  for (int i=0; i<num_points; i++)
  {
    float dist = sqrtf ((sx - points_x[i]) * (sx - points_x[i])
        + (sz - points_z[i]) * (sz - points_z[i])
        );
    if (dist < 1e-9) {
      at_point_indx = i;
    } else {
      weight[i]   = 1.0 / dist;
      total_weight += weight[i];
    }
  }
  // if at a point
  if (at_point_indx > 0) {
    total_weight = 1.0;
    // other weight 0
    for (int i=0; i<num_points; i++) {
      weight[i] = 0.0;
    }
    // point weight 1
    weight[at_point_indx] = 1.0;
  }

  // interp

  *si_curv = 0.0;
  *sk_curv = 0.0;

  for (int i=0; i<num_points; i++)
  {
    weight[i] *= 1.0 / total_weight ;

    (*si_curv) += weight[i] * points_i[i];
    (*sk_curv) += weight[i] * points_k[i];  

    fprintf(stdout,"---- i=%d,weight=%f,points_i=%f,points_k=%f\n",
        i,weight[i],points_i[i],points_k[i]);
  }

  return 0;
}

float
gd_coord_get_x(gd_t *gd, int i, int k)
{
  float var = 0.0;

  if (gd->type == GD_TYPE_CART)
  {
    var = gd->x1d[i];
  }
  else if (gd->type == GD_TYPE_CURV)
  {
    size_t iptr = i + k * gd->siz_iz;
    var = gd->x2d[iptr];
  }

  return var;
}

float
gd_coord_get_z(gd_t *gd, int i, int k)
{
  float var = 0.0;

  if (gd->type == GD_TYPE_CART)
  {
    var = gd->z1d[k];
  }
  else if (gd->type == GD_TYPE_CURV)
  {
    size_t iptr = i + k * gd->siz_iz;
    var = gd->z2d[iptr];
  }

  return var;
}


// for Lebedev scheme
void gd_curv_rearrange_c2d(gd_t *gdcurv)
{

  float *x2d = gdcurv->x2d;
  float *z2d = gdcurv->z2d;
  float *x2d_lebedev = gdcurv->x2d_lebedev;
  float *z2d_lebedev = gdcurv->z2d_lebedev;

  size_t siz_icmp_qr = gdcurv->siz_icmp * 0.25;

  int nx = gdcurv->nx;
  int nz = gdcurv->nz;
  size_t siz_iz = gdcurv->siz_iz;


  float *x2d1 = x2d_lebedev + 0 * siz_icmp_qr;
  float *x2d2 = x2d_lebedev + 1 * siz_icmp_qr;
  float *x2d3 = x2d_lebedev + 2 * siz_icmp_qr;
  float *x2d4 = x2d_lebedev + 3 * siz_icmp_qr;

  float *z2d1 = z2d_lebedev + 0 * siz_icmp_qr;
  float *z2d2 = z2d_lebedev + 1 * siz_icmp_qr;
  float *z2d3 = z2d_lebedev + 2 * siz_icmp_qr;
  float *z2d4 = z2d_lebedev + 3 * siz_icmp_qr;
  
  gd_rearrange(x2d, x2d1, x2d2, x2d3, x2d4, nx, nz, siz_iz);
  gd_rearrange(z2d, z2d1, z2d2, z2d3, z2d4, nx, nz, siz_iz);
}
  
void gd_rearrange(float *i3d, float *i3d1, float *i3d2, float *i3d3, float *i3d4,
		int nx, int ny, size_t siz_line)
{
	// just for test
  if ( (nx%2==0) && (ny%2==0) ){
	  int countx;
	  countx = nx*0.5;
	  int iptr1, iptr2;

	  for (size_t i=0; i<nx*0.5; i++){
		  for (size_t j=0; j < ny*0.5; j++){
			  iptr1 = i+j*countx;
			  iptr2 = 2*i+2*j*siz_line;
			  i3d1[iptr1] = i3d[iptr2];
			  i3d2[iptr1] = i3d[iptr2+1];
			  i3d3[iptr1] = i3d[iptr2+siz_line];
			  i3d4[iptr1] = i3d[iptr2+siz_line+1];
		  }
	  }
  }
}

void gd_curv_rearrange_m2d(gd_t *gdcurv, md_t *m2d)
{

    size_t siz_icmp_qr = gdcurv->siz_icmp * 0.25;
    int nx = gdcurv->nx;
    int nz = gdcurv->nz;
    size_t siz_iz = gdcurv->siz_iz;

    float *c_11 = m2d->c11;
    float *c_13 = m2d->c13;
    float *c_15 = m2d->c15;
    float *c_33 = m2d->c33;
    float *c_35 = m2d->c35;
    float *c_55 = m2d->c55;
    float *rho  = m2d->rho;

    float *c11_1 = m2d->c11_lebedev + 0 * siz_icmp_qr;
    float *c11_2 = m2d->c11_lebedev + 1 * siz_icmp_qr;
    float *c11_3 = m2d->c11_lebedev + 2 * siz_icmp_qr;
    float *c11_4 = m2d->c11_lebedev + 3 * siz_icmp_qr;

    float *c13_1 = m2d->c13_lebedev + 0 * siz_icmp_qr;
    float *c13_2 = m2d->c13_lebedev + 1 * siz_icmp_qr;
    float *c13_3 = m2d->c13_lebedev + 2 * siz_icmp_qr;
    float *c13_4 = m2d->c13_lebedev + 3 * siz_icmp_qr;

    float *c15_1 = m2d->c15_lebedev + 0 * siz_icmp_qr;
    float *c15_2 = m2d->c15_lebedev + 1 * siz_icmp_qr;
    float *c15_3 = m2d->c15_lebedev + 2 * siz_icmp_qr;
    float *c15_4 = m2d->c15_lebedev + 3 * siz_icmp_qr;

    float *c33_1 = m2d->c33_lebedev + 0 * siz_icmp_qr;
    float *c33_2 = m2d->c33_lebedev + 1 * siz_icmp_qr;
    float *c33_3 = m2d->c33_lebedev + 2 * siz_icmp_qr;
    float *c33_4 = m2d->c33_lebedev + 3 * siz_icmp_qr;

    float *c35_1 = m2d->c35_lebedev + 0 * siz_icmp_qr;
    float *c35_2 = m2d->c35_lebedev + 1 * siz_icmp_qr;
    float *c35_3 = m2d->c35_lebedev + 2 * siz_icmp_qr;
    float *c35_4 = m2d->c35_lebedev + 3 * siz_icmp_qr;

    float *c55_1 = m2d->c55_lebedev + 0 * siz_icmp_qr;
    float *c55_2 = m2d->c55_lebedev + 1 * siz_icmp_qr;
    float *c55_3 = m2d->c55_lebedev + 2 * siz_icmp_qr;
    float *c55_4 = m2d->c55_lebedev + 3 * siz_icmp_qr;
	
    float *rho1  = m2d->rho_lebedev + 0 * siz_icmp_qr;
    float *rho2  = m2d->rho_lebedev + 1 * siz_icmp_qr;
    float *rho3  = m2d->rho_lebedev + 2 * siz_icmp_qr;
    float *rho4  = m2d->rho_lebedev + 3 * siz_icmp_qr;

    
    gd_rearrange(c_11, c11_1, c11_2, c11_3, c11_4, nx, nz, siz_iz);
    gd_rearrange(c_13, c13_1, c13_2, c13_3, c13_4, nx, nz, siz_iz);
    gd_rearrange(c_15, c15_1, c15_2, c15_3, c15_4, nx, nz, siz_iz);
    gd_rearrange(c_33, c33_1, c33_2, c33_3, c33_4, nx, nz, siz_iz);
    gd_rearrange(c_35, c35_1, c35_2, c35_3, c35_4, nx, nz, siz_iz);
    gd_rearrange(c_55, c55_1, c55_2, c55_3, c55_4, nx, nz, siz_iz);
    gd_rearrange(rho, rho1, rho2, rho3, rho4, nx, nz, siz_iz);
}

