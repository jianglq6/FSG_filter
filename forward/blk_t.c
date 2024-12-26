/*********************************************************************
 * setup fd operators
 **********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "fdlib_mem.h"
#include "fdlib_math.h"
#include "blk_t.h"

//
// malloc inner vars
//

int
blk_init(blk_t *blk, const int verbose)
{
  int ierr = 0;

  // alloc struct vars
  blk->fd            = (fd_t *)malloc(sizeof(fd_t));
  blk->gdinfo        = (gdinfo_t *)malloc(sizeof(gdinfo_t));
  blk->gd            = (gd_t        *)malloc(sizeof(gd_t     ));
  blk->gdcurv_metric = (gdcurv_metric_t *)malloc(sizeof(gdcurv_metric_t));
  blk->md            = (md_t      *)malloc(sizeof(md_t     ));
  blk->wav           = (wav_t      *)malloc(sizeof(wav_t     ));
  blk->src           = (src_t      *)malloc(sizeof(src_t     ));

  sprintf(blk->name, "%s", "single");

  return ierr;
}

// set str
int
blk_set_output(blk_t *blk,
               char *output_dir,
               char *grid_export_dir,
               char *media_export_dir,
              const int verbose)
{
  // set name
  //sprintf(blk->name, "%s", name);

  // output
  sprintf(blk->output_dir, "%s", output_dir);
  sprintf(blk->grid_export_dir, "%s", grid_export_dir);
  sprintf(blk->media_export_dir, "%s", media_export_dir);

  return 0;
}


/*********************************************************************
 * estimate dt
 *********************************************************************/

int
blk_dt_esti(float dx, float dz, gdinfo_t *gdinfo, md_t *md,
    float CFL, float *dtmax, float *dtmaxVp, float *dtmaxL,
    int *dtmaxi, int *dtmaxk)
{
  int ierr = 0;

  float dtmax_local = 1.0e10;
  float Vp;

  float dtLe = dx*dz / sqrt(dx*dx+dz*dz);

  for (int k = gdinfo->nk1; k < gdinfo->nk2; k++)
  {
      for (int i = gdinfo->ni1; i < gdinfo->ni2; i++)
      {
        size_t iptr = i + k * gdinfo->siz_iz;

        if (md_is_el_aniso(md)==1) {
          Vp = sqrt( (md->c11[iptr]) / md->rho[iptr] );
        } else {
          fprintf(stderr,"ERROR: medium type is not implemented\n");
          exit(1);
        }

        // convert to dt
        float dt_point = CFL / Vp * dtLe;

        // if smaller
        if (dt_point < dtmax_local) {
          dtmax_local = dt_point;
          *dtmaxi = i;
          *dtmaxk = k;
          *dtmaxVp = Vp;
        }

      } // i
  } //k

  *dtmax  = dtmax_local;
  *dtmaxL = dtLe;

  return ierr;
}

float
blk_keep_two_digi(float dt)
{
  char str[40];
  float dt_2;

  sprintf(str, "%4.2e", dt);

  str[3] = '0';

  sscanf(str, "%f", &dt_2);
  
  return dt_2;
}
