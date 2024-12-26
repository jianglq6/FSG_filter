/*********************************************************************
 * setup fd operators
 **********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "constants.h"
#include "fdlib_mem.h"
#include "fd_t.h"

/*
 * set Lebedev scheme
 */
int 
fd_set_lebedev(fd_t *fd, int is_filter, float filter_sigma)
{

  int ierr = 0;
  //----------------------------------------------------------------------------
  // lebedev scheme
  //----------------------------------------------------------------------------
  
  // set len

  fd->fdx_nghosts = 4;
  fd->fdz_nghosts = 4;
  fd->lebedev_len = 8;
  fd->lebedev_half_len = 4;

  fd->is_filter = is_filter;
  fd->filter_sigma = filter_sigma;
  fd->filter_coef_len = 5; // + center points
  fd->filter_coef_half_len = 4; 

#define m_lg_max_len 12
#define m_lg_num_lay 6
#define m_filter_len 7
#define m_filter_lay 3
    
  float lebedev_all_coef[m_lg_num_lay][m_lg_max_len] =
  {
  	{-1.0, 
  	  1.0, 
  	  0.0, 
  	  0.0,
  	  0.0, 
  	  0.0, 
  	  0.0, 
  	  0.0,
  	  0.0, 
  	  0.0, 
  	  0.0, 
  	  0.0},
  
  	{ 0.4301412e-1, 
  	 -0.1129042e+1,
  	  0.1129042e+1,  
  	 -0.4301412e-1, 
  	  0.0,  
  	  0.0, 
  	  0.0,
  	  0.0, 
  	  0.0, 
  	  0.0, 
  	  0.0, 
  	  0.0},
  
  	{-0.6301572e-2, 
  	  0.7249965e-1, 
  	 -0.1185991e+1, 
  	  0.1185991e+1, 
  	 -0.7249965e-1,  
  	  0.6301572e-2, 
  	  0.0, 
  	  0.0, 
  	  0.0,
  	  0.0, 
  	  0.0, 
  	  0.0},

 // 4th-order, for we test
  	{ 0.1700324e-2,
  	 -0.1507536e-1, 
  	  0.9382142e-1, 
  	 -0.1217990e+1,
  	  0.1217990e+1,
  	 -0.9382142e-1, 
  	  0.1507536e-1, 
  	 -0.1700324e-2,
  	  0.0, 
  	  0.0,
  	  0.0, 
  	  0.0}, 
  
  	{-0.6814783e-3,
  	  0.5033546e-2,
  	 -0.2343440e-1,
  	  0.1082265e+0,
  	 -0.1236607e+1,
  	  0.1236607e+1,
  	 -0.1082265e+0,
  	  0.2343440e-1,
  	 -0.5033546e-2,
  	  0.6814783e-3,
  	  0.0,
  	  0.0},
  
  	{ 0.3462075e-3,
  	 -0.2215897e-2,
  	  0.8719078e-2,
  	 -0.2997970e-1,
  	  0.1175538e+0,
  	 -0.1247662e+1,
  	  0.1247662e+1,
  	 -0.1175538e+0,
  	  0.2997970e-1,
  	 -0.8719078e-2,
  	  0.2215897e-2,
  	 -0.3462075e-3}
  };
  
  
  
  float filter_coef[m_filter_lay][m_filter_len] = 
  {
  	{0.243527493120, -0.204788880640, 0.120007591680, -0.045211119360, 0.008228661760, 0.0, 0.0},
  	{0.215044884112, -0.187772883589, 0.123755948787, -0.059227575576, 0.018721609157, -0.002999540835, 0.0},
  	{0.190899511506, -0.171503832236, 0.123632891797, -0.069975429105, 0.029662754736, -0.008520738659, 0.001254597714}
  };

  // sum(|c_m|)^-1, for 4th-order, 1D
  fd->CFL = 0.752;
  fd->lebedev_coef = (float *) fdlib_mem_malloc_1d(fd->lebedev_len*sizeof(float),"fd_set_lebedev");
    
  fd->filter_coef = (float *) fdlib_mem_malloc_1d((fd->filter_coef_len)*sizeof(float),"fd_set_filter");
  
  // fd_print 
  fprintf(stdout, "  fd_lebedev_len %d\n", fd->lebedev_len);
  fprintf(stdout, "  fd_lebedev_half_len %d\n", fd->lebedev_half_len);
  for (int i=0; i<fd->lebedev_len; i++) {
  	fd->lebedev_coef[i] = lebedev_all_coef[fd->lebedev_half_len-1][i];
  	fprintf(stdout, "  fd_coef[%d] = %f\n", i, fd->lebedev_coef[i]);
  }
  for (int i=0; i<fd->filter_coef_len; i++) {
  	fd->filter_coef[i] = filter_coef[fd->filter_coef_len-5][i];
  	printf("  filter_coef[%d] = %f\n", i, filter_coef[0][i]);
  }

  return ierr;
}



