#ifndef GD_CURV_H
#define GD_CURV_H

#include "constants.h"
#include "gd_info.h"
#include "md_t.h"

/*************************************************
 * structure
 *************************************************/

typedef enum {

  GD_TYPE_CART = 1,
  GD_TYPE_VMAP = 2,
  GD_TYPE_CURV = 3

} gd_type_t;

//  grid coordinate for both cart, vmap and curv
//    to reduce duplicated functions
typedef struct {

  gd_type_t type;

  int n1, n2, n3;
  int nx, nz, ncmp;
  float *v4d; // allocated var

  //to avoid ref x3d at different funcs
  float *x2d; // pointer to var
  float *z2d;
  // rearrange for lebedev
  float *x2d_lebedev; 
  float *z2d_lebedev;

  // for cart grid
  float *x1d;
  float *z1d;
  float dx;
  float dz;
  // x0/y0/z0 for grid gen
  float x0_glob;
  float z0_glob;

  // min/max of this thread including ghost points
  float xmin, xmax;
  float zmin, zmax;

  size_t siz_iz;
  size_t siz_icmp;

  size_t *cmp_pos;
  char  **cmp_name;
} gd_t;


//  for metric
typedef struct {
  int n1, n2, n3, n4;
  int nx, nz, ncmp;
  float *v4d; // allocated var

  float *jac; // pointer to var
  float *xi_x;
  float *xi_z;
  float *zeta_x;
  float *zeta_z;

  size_t siz_iz;
  size_t siz_icmp;

  size_t *cmp_pos;
  char  **cmp_name;
} gdcurv_metric_t;


/*************************************************
 * function prototype
 *************************************************/

void 
gd_curv_init(gdinfo_t *gdinfo, gd_t *gdcurv);

void 
gd_curv_metric_init(gdinfo_t        *gdinfo,
                    gdcurv_metric_t *metric);
void
gd_curv_metric_cal(gdinfo_t        *gdinfo,
                   gd_t        *gdcurv,
                   gdcurv_metric_t *metric,
                   int fd_len, int *restrict fd_indx, float *restrict fd_coef);

void
gd_curv_gen_cart(
  gdinfo_t *gdinfo,
  gd_t *gdcurv,
  float dx, float x0_glob,
  float dz, float z0_glob,
  float *x1d, float *z1d);

void
gd_curv_metric_import(gdcurv_metric_t *metric, char *import_dir);

void
gd_curv_coord_import(gd_t *gdcurv, char *import_dir);

void
gd_curv_coord_export(
  gdinfo_t *gdinfo,
  gd_t *gdcurv,
  char *output_dir);

void
gd_cart_coord_export(
  gdinfo_t *gdinfo,
  gd_t *gdcart,
  char *output_dir);

void
gd_curv_metric_export(gdinfo_t        *gdinfo,
                      gdcurv_metric_t *metric,
                      char *output_dir);


int gd_grid_z_interp(float *z3dpart, float *zlayerpart, int *NCellPerlay,
                     int *VmapSpacingIsequal, int nLayers, int nx, int ny);

float gd_seval(int ni, float u,
            int n, float x[], float y[],
            float b[], float c[], float d[],
            int *last);

int gd_SPLine( int n, int end1, int end2,
           float slope1, float slope2,
           float x[], float y[],
           float b[], float c[], float d[],
           int *iflag);

void gd_SPL(int n, float *x, float *y, int ni, float *xi, float *yi);

void
gd_curv_set_minmax(gd_t *gdcurv);

void 
gd_cart_init_set(gdinfo_t *gdinfo, gd_t *gdcart,
  float dx, float x0_glob,
  float dz, float z0_glob);

int
gd_cart_coord_to_local_indx(gdinfo_t *gdinfo,
                           gd_t *gdcart,
                           float sx,
                           float sz,
                           int   *ou_si, int *ou_sk,
                           float *ou_sx_inc, float *ou_sz_inc);

int
gd_curv_coord_to_local_indx(gdinfo_t *gdinfo,
                        gd_t *gd,
                        float sx, float sz,
                        int *si, int *sk,
                        float *sx_inc, float *sz_inc,
                        float *restrict wrk3d);

  int
gd_curv_coord2index_sample(float sx, float sz, 
    float *points_x, // x coord of all points
    float *points_z,
    float *points_i, // curv coord of all points
    float *points_k,
    int    nx_sample,
    int    nz_sample,
    float *si_curv, // interped curv coord
    float *sk_curv);

float
gd_coord_get_x(gd_t *gd, int i, int k);

float
gd_coord_get_z(gd_t *gd, int i, int k);

void 
gd_curv_init_lebedev(gdinfo_t *gdinfo, gd_t *gdcurv);

// for Lebedev scheme
void gd_curv_rearrange_c2d(gd_t *gdcurv);
void gd_rearrange(float *i3d, float *i3d1, float *i3d2, float *i3d3, float *i3d4,
		int nx, int ny, size_t siz_line);

void gd_curv_cal_metric_lebedev(gdinfo_t *gdinfo,
		                        gd_t     *gdcurv,
								gdcurv_metric_t *metric,
								int fd_lebedev_half_len, float *fd_lebedev_coef);
void gd_cal_metric(
		size_t iptr, float *jac2d, 
		float  x_xi, float  x_zt, float  z_xi, float  z_zt,
		float *xi_x, float *xi_z, float *zt_x, float *zt_z);

void gd_curv_rearrange_m2d(gd_t *gdcurv, md_t *m2d);

#endif
