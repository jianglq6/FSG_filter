/*********************************************************************
 * setup fd operators
 **********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fdlib_mem.h"
#include "constants.h"
#include "gd_info.h"

int
gd_info_set_lebedev(
		gdinfo_t *const gdinfo,
        const int number_of_total_grid_points_x,
        const int number_of_total_grid_points_z,
              int abs_num_of_layers[][2],
        const int fdx_nghosts,
        const int fdz_nghosts,
        const int verbose)
{
  int ierr = 0;

  int ni = number_of_total_grid_points_x * 2;
  int nk = number_of_total_grid_points_z * 2;
  
  // add ghost points
  int nx = ni + 2 * fdx_nghosts * 2;
  int nz = nk + 2 * fdz_nghosts * 2;

  gdinfo->ni = ni;
  gdinfo->nk = nk;

  gdinfo->nx = nx;
  gdinfo->nz = nz;

  gdinfo->ni1 = fdx_nghosts * 2;
  gdinfo->ni2 = gdinfo->ni1 + ni - 1;

  gdinfo->nk1 = fdz_nghosts * 2;
  gdinfo->nk2 = gdinfo->nk1 + nk - 1;

  gdinfo->lni = ni * 0.5;
  gdinfo->lnk = nk * 0.5;

  gdinfo->lnx = nx * 0.5;
  gdinfo->lnz = nz * 0.5;

  gdinfo->lni1 = fdx_nghosts;
  gdinfo->lni2 = gdinfo->lni1 + gdinfo->lni - 1;

  gdinfo->lnk1 = fdz_nghosts;
  gdinfo->lnk2 = gdinfo->lnk1 + gdinfo->lnk - 1;

  // new var, will replace above old naming
  gdinfo->siz_iz   = gdinfo->nx;
  gdinfo->siz_icmp = gdinfo->nx * gdinfo->nz;

  // set npoint_ghosts according to fdz_nghosts
  gdinfo->npoint_ghosts = fdz_nghosts * 2;

  gdinfo->fdx_nghosts = fdx_nghosts * 2;
  gdinfo->fdz_nghosts = fdz_nghosts * 2;

  gdinfo->index_name = fdlib_mem_malloc_2l_char(
                        CONST_NDIM, CONST_MAX_STRLEN, "gdinfo name");

  // grid coord name
  sprintf(gdinfo->index_name[0],"%s","i");
  sprintf(gdinfo->index_name[1],"%s","k");

  return ierr;
}
//
// set grid size
//

int
gd_info_set(gdinfo_t *const gdinfo,
            const int number_of_total_grid_points_x,
            const int number_of_total_grid_points_z,
                  int abs_num_of_layers[][2],
            const int fdx_nghosts,
            const int fdz_nghosts,
            const int verbose)
{
  int ierr = 0;

  int ni = number_of_total_grid_points_x;
  int nk = number_of_total_grid_points_z;
  
  // add ghost points
  int nx = ni + 2 * fdx_nghosts;
  int nz = nk + 2 * fdz_nghosts;

  gdinfo->ni = ni;
  gdinfo->nk = nk;

  gdinfo->nx = nx;
  gdinfo->nz = nz;

  gdinfo->ni1 = fdx_nghosts;
  gdinfo->ni2 = gdinfo->ni1 + ni - 1;

  gdinfo->nk1 = fdz_nghosts;
  gdinfo->nk2 = gdinfo->nk1 + nk - 1;

  // new var, will replace above old naming
  gdinfo->siz_iz   = gdinfo->nx;
  gdinfo->siz_icmp = gdinfo->nx * gdinfo->nz;

  // set npoint_ghosts according to fdz_nghosts
  gdinfo->npoint_ghosts = fdz_nghosts;

  gdinfo->fdx_nghosts = fdx_nghosts;
  gdinfo->fdz_nghosts = fdz_nghosts;

  gdinfo->index_name = fdlib_mem_malloc_2l_char(
                        CONST_NDIM, CONST_MAX_STRLEN, "gdinfo name");

  // grid coord name
  sprintf(gdinfo->index_name[0],"%s","i");
  sprintf(gdinfo->index_name[1],"%s","k");

  return ierr;
}

/*
 * give a local index ref, check if in this thread
 */

int
gd_info_lindx_is_inner(int i, int k, gdinfo_t *gdinfo)
{
  int is_in = 0;

  if (   i >= gdinfo->ni1 && i <= gdinfo->ni2
      && k >= gdinfo->nk1 && k <= gdinfo->nk2)
  {
    is_in = 1;
  }

  return is_in;
}  

int
gd_info_pindx_is_inner(int i_phy, int k_phy, gdinfo_t *gdinfo)
{
  int is_in = 0;

  int i = i_phy + gdinfo->fdx_nghosts;
  int k = k_phy + gdinfo->fdz_nghosts;

  if (   i >= gdinfo->ni1 && i <= gdinfo->ni2
      && k >= gdinfo->nk1 && k <= gdinfo->nk2)
  {
    is_in = 1;
  }

  return is_in;
}  

int
gd_info_pindx_is_inner_i(int i_phy, gdinfo_t *gdinfo)
{
  int is_in = 0;

  int i = i_phy + gdinfo->fdx_nghosts;

  if (   i >= gdinfo->ni1 && i <= gdinfo->ni2)
  {
    is_in = 1;
  }

  return is_in;
}  

int
gd_info_pindx_is_inner_k(int k_phy, gdinfo_t *gdinfo)
{
  int is_in = 0;

  int k = k_phy + gdinfo->fdz_nghosts;

  if ( k >= gdinfo->nk1 && k <= gdinfo->nk2)
  {
    is_in = 1;
  }

  return is_in;
}  

/*
 * print for QC
 */

int
gd_info_print(gdinfo_t *gdinfo)
{    
  fprintf(stdout, "-------------------------------------------------------\n");
  fprintf(stdout, "--> grid info:\n");
  fprintf(stdout, "-------------------------------------------------------\n");
  fprintf(stdout, " nx    = %-10d\n", gdinfo->nx);
  fprintf(stdout, " nz    = %-10d\n", gdinfo->nz);
  fprintf(stdout, " ni    = %-10d\n", gdinfo->ni);
  fprintf(stdout, " nk    = %-10d\n", gdinfo->nk);

  fprintf(stdout, " ni1   = %-10d\n", gdinfo->ni1);
  fprintf(stdout, " ni2   = %-10d\n", gdinfo->ni2);
  fprintf(stdout, " nk1   = %-10d\n", gdinfo->nk1);
  fprintf(stdout, " nk2   = %-10d\n", gdinfo->nk2);

  return(0);
}
