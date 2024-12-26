/*********************************************************************
 * wavefield for 2d elastic 1st-order equations
 **********************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "constants.h"
#include "fdlib_mem.h"
#include "wav_t.h"

int 
wav_init_lebedev(gdinfo_t *gdinfo,
                 wav_t *V)
{
  int ierr = 0;

  // for Lebedev
  V->ncmp = 10;

  V->nx   = gdinfo->nx;
  V->nz   = gdinfo->nz;

  V->siz_iz   = V->nx;
  V->siz_icmp = V->nx * V->nz;

  // vars
  // 3 Vi, 6 Tij, 4 rk stages
  V->v5d = (float *) fdlib_mem_calloc_1d_float(V->siz_icmp * V->ncmp,
                        0.0, "v5d, wf_el3d_1st");
  return ierr;
}


int
wav_check_value(float *restrict w, wav_t *wav)
{
  int ierr = 0;

  for (int icmp=0; icmp < wav->ncmp; icmp++)
  {
    float *ptr = w + icmp * wav->siz_icmp;
    for (size_t iptr=0; iptr < wav->siz_icmp; iptr++)
    {
      if (ptr[iptr] != ptr[iptr])
      {
        fprintf(stderr, "ERROR: NaN occurs at iptr=%d icmp=%d\n", iptr, icmp);
        fflush(stderr);
        exit(-1);
      }
    }
  }

  return ierr;
}

int
wav_zero_edge(gdinfo_t *gdinfo, wav_t *wav,
                                  float *restrict w4d)
{
  int ierr = 0;

  for (int icmp=0; icmp < wav->ncmp; icmp++)
  {
    float *restrict var = w4d + wav->cmp_pos[icmp];

    // z1
    for (int k=0; k < gdinfo->nk1; k++)
    {
      size_t iptr_k = k * gdinfo->siz_iz;
        for (int i=0; i < gdinfo->nx; i++)
        {
          size_t iptr = iptr_k + i;
          var[iptr] = 0.0; 
        }
    }

    // z2
    for (int k=gdinfo->nk2+1; k < gdinfo->nz; k++)
    {
      size_t iptr_k = k * gdinfo->siz_iz;
        for (int i=0; i < gdinfo->nx; i++)
        {
          size_t iptr = iptr_k + i;
          var[iptr] = 0.0; 
        }
    }

    // x1
    for (int k = gdinfo->nk1; k <= gdinfo->nk2; k++)
    {
      size_t iptr_k = k * gdinfo->siz_iz;
        for (int i=0; i < gdinfo->ni1; i++)
        {
          size_t iptr = iptr_k + i;
          var[iptr] = 0.0; 
        }
    } 

    // x2
    for (int k = gdinfo->nk1; k <= gdinfo->nk2; k++)
    {
      size_t iptr_k = k * gdinfo->siz_iz;
        for (int i = gdinfo->ni2+1; i < gdinfo->nx; i++)
        {
          size_t iptr = iptr_k + i;
          var[iptr] = 0.0; 
        }
    } 

  } // icmp

  return ierr;
}
