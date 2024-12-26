#ifndef BLK_T_H
#define BLK_T_H

#include "constants.h"
#include "fd_t.h"
#include "gd_info.h"
#include "gd_t.h"
#include "md_t.h"
#include "wav_t.h"
#include "src_t.h"

/*******************************************************************************
 * structure
 ******************************************************************************/

typedef struct
{
  // name for output file name
  char name[CONST_MAX_STRLEN];

  //// flag of medium
  //int medium_type;

  // fd
  fd_t    *fd;     // collocated grid fd

  // grid index info
  gdinfo_t *gdinfo;
  
  // coordnate: x3d, y3d, z3d
  gd_t *gd;

  // grid metrics: jac, xi_x, etc
  gdcurv_metric_t *gdcurv_metric;

  // media: rho, lambda, mu etc
  md_t *md;

  // wavefield:
  wav_t *wav;
  
  // source term
  src_t *src;
  
  // io
//  iorecv_t  *iorecv;
//  ioline_t  *ioline;
//  iosnap_t  *iosnap;

  // fname and dir
  char output_fname_part[CONST_MAX_STRLEN];
  // wavefield output
  char output_dir[CONST_MAX_STRLEN];
  // seperate grid output to save grid for repeat simulation
  char grid_export_dir[CONST_MAX_STRLEN];
  // seperate medium output to save medium for repeat simulation
  char media_export_dir[CONST_MAX_STRLEN];

  // mem usage
  size_t number_of_float;
  size_t number_of_btye;
} blk_t;

/*******************************************************************************
 * function prototype
 ******************************************************************************/

int
blk_init(blk_t *blk, const int verbose);

// set str
int
blk_set_output(blk_t *blk,
               char *output_dir,
               char *grid_export_dir,
               char *media_export_dir,
               const int verbose);

int
blk_print(blk_t *blk);

int
blk_dt_esti(float dx, float dz, gdinfo_t *gdinfo, md_t *md,
    float CFL, float *dtmax, float *dtmaxVp, float *dtmaxL,
    int *dtmaxi, int *dtmaxk);

float
blk_keep_two_digi(float dt);

#endif
