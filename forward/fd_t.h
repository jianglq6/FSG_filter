#ifndef FD_T_H
#define FD_T_H

#define FD_STG_MAX_LEN 4

// use siz_shift to find adjacent point of the stentil for 2d var
#define M_FD_SHIFT(deriv, var, iptr, fd_length, fd_shift, fd_coef, n) \
   deriv = fd_coef[0] * var[iptr + fd_shift[0]]; \
   for (n=1; n<fd_length; n++) { \
       deriv += fd_coef[n] * var[iptr + fd_shift[n]]; \
   }

// 2D lebedev 
//-------------------------------

#define Dx0(var, fdx_coef, n, ix, iy, line)   \
  (fdx_coef[n-1] * ( var[iy*line + ix + (n-1)] - var[iy*line + ix - n] ) )

#define Dx1(var, fdx_coef, n, ix, iy, line)   \
  (fdx_coef[n-1] * ( var[iy*line + ix + n] - var[iy*line + ix - (n-1)] ) )

#define Dz0(var, fdy_coef, n, ix, iy, line)   \
  (fdy_coef[n-1] * ( var[(iy + (n-1))*line + ix  ] - var[(iy -n)*line + ix] ) )

#define Dz1(var, fdy_coef, n, ix, iy, line)   \
  (fdy_coef[n-1] * ( var[(iy + n)*line + ix  ] - var[(iy -(n-1))*line + ix] ) )

/*******************************************************************************
 * structure for different fd schemes
 ******************************************************************************/

/*
 * elementary operator
 */

typedef struct
{
  int total_len;
  int half_len;
  int left_len;
  int right_len;
  int   *indx;  // indx change to cur point as 0 for 1d
  int   *shift; // num of grid points skipped
  float *coef;
} fd_op_t; 

/*
 * fsg scheme
 */

typedef struct {

  float CFL; // 1d cfl value for the scheme
  //----------------------------------------------------------------------------
  // Lebedev scheme
  //----------------------------------------------------------------------------

  // ghost point required 
  int fdx_nghosts;
  int fdz_nghosts;

  int   lebedev_len;
  int 	lebedev_half_len;
  int   filter_coef_len;
  int   filter_coef_half_len;

  int     is_filter;
  float   filter_sigma;
  float  *lebedev_coef;
  float  *filter_coef;

} fd_t;


/*******************************************************************************
 * function prototype
 ******************************************************************************/


void
fd_print(fd_t *fd);

int 
fd_set_lebedev(fd_t *fd, int is_filter, float filter_sigma);

#endif
