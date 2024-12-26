/**********************************************************************************
 * Fully staggered Grid Finite Difference Method for 2D Wave Propagation Simulation 
 **********************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stddef.h>
#include <time.h>

#include "constants.h"
#include "par_t.h"
#include "blk_t.h"

#include "media_discrete_model.h"
#include "sv_eq1st_curv_fsg_el_aniso.h"

int main(int argc, char** argv)
{
  int verbose = 1; // default fprint
  char *par_fname;
  char err_message[CONST_MAX_STRLEN];

//-------------------------------------------------------------------------------
// get commond-line argument
//-------------------------------------------------------------------------------

  // argc checking
  if (argc < 2) {
    fprintf(stdout,"usage: main_curv_fsg_2d <par_file> <opt: verbose>\n");
    exit(1);
  }

  par_fname = argv[1];

  if (argc >= 3) {
    verbose = atoi(argv[2]); // verbose number
    fprintf(stdout,"verbose=%d\n", verbose); fflush(stdout);
  }

  fprintf(stdout,"par file =  %s\n", par_fname); 

  // read par

  par_t *par = (par_t *) malloc(sizeof(par_t));

  par_read_from_file(par_fname, par, verbose);

  if (verbose>0) par_print(par);

//-------------------------------------------------------------------------------
// init blk_t
//-------------------------------------------------------------------------------

  if (verbose>0) fprintf(stdout,"create blk ...\n"); 

  // malloc blk
  blk_t *blk = (blk_t *) malloc(sizeof(blk_t));

  // malloc inner vars
  blk_init(blk, verbose);

  fd_t            *fd            = blk->fd    ;
  gdinfo_t        *gdinfo        = blk->gdinfo;
  gd_t            *gdcurv        = blk->gd;
  md_t            *md            = blk->md;
  wav_t           *wav           = blk->wav;
  src_t           *src           = blk->src;

  // set up fd_t
  //    not support selection scheme by par file yet
  if (verbose>0) fprintf(stdout,"set scheme ...\n"); 
  fd_set_lebedev(fd, par->is_filter, par->filter_sigma);

  // set gdinfo
  gd_info_set_lebedev(gdinfo, 
                      par->number_of_total_grid_points_x,
                      par->number_of_total_grid_points_z,
                      par->abs_num_of_layers,
                      fd->fdx_nghosts,
                      fd->fdz_nghosts,
                      verbose);

  // set str in blk
  blk_set_output(blk, 
                 par->output_dir,
                 par->grid_export_dir,
                 par->media_export_dir,
                 verbose);


//-------------------------------------------------------------------------------
//-- grid generation or import
//-------------------------------------------------------------------------------

  float dx = par->cartesian_grid_stepsize[0];
  float dz = par->cartesian_grid_stepsize[1];

  if (verbose>0) fprintf(stdout,"allocate grid vars ...\n"); 

  // malloc var in gdcurv
  gd_curv_init_lebedev(gdinfo, gdcurv);


  float *x1d = (float*) malloc(gdcurv->nx*sizeof(float));
  float *z1d = (float*) malloc(gdcurv->nz*sizeof(float));

  // generate grid coord
  switch (par->grid_generation_itype)
  {
    case PAR_GRID_CARTESIAN :

        fprintf(stdout,"generate cartesian grid in code ...\n"); 


        float x0 = par->cartesian_grid_origin[0];
        float z0 = par->cartesian_grid_origin[1];

        // for lebedev
        gd_curv_gen_cart(gdinfo,gdcurv,0.5*dx,x0,0.5*dz,z0,x1d,z1d);

        break;

    case PAR_GRID_IMPORT:

        fprintf(stdout,"import grid vars ...\n"); 
//       gd_curv_coord_import(gdcurv, par->grid_import_dir);

        break;

  }

  // cal min/max of this thread
  gd_curv_set_minmax(gdcurv);


  // rearrange coord for lebedev
  gd_curv_rearrange_c2d(gdcurv);

 
//-------------------------------------------------------------------------------
//-- media generation or import
//-------------------------------------------------------------------------------

  // allocate media vars
  if (verbose>0) {fprintf(stdout,"allocate media vars ...\n"); fflush(stdout);}
  md_init_lebedev(gdinfo, md);

  // read or discrete velocity model
  switch (par->media_input_itype)
  {

    case PAR_MEDIA_ASCII :
        if (verbose>0) fprintf(stdout,"import discrete medium file ...\n"); 
        media_read_test_bptti(md, par->media_input_file);

        break;

    case PAR_MEDIA_LAY : {
        if (verbose>0) fprintf(stdout,"read and discretize layer medium file ...\n"); 

        media_layer2model_el_aniso(md->rho, md->c11, md->c13,
                                   md->c15, md->c33, md->c35, 
                                   md->c55, x1d, z1d, md->nx, md->nz, 
                                   MEDIA_USE_CART, par->media_input_file,
                                   par->equivalent_medium_method); 
        break;
    }

  }

  // export grid media
  if (par->is_export_media==1)
  {
    if (verbose>0) fprintf(stdout,"export discrete medium to file ...\n"); 

    md_export(gdinfo, md, blk->media_export_dir);
  } else {
    if (verbose>0) fprintf(stdout,"do not export medium\n"); 
  }


  // for lebedev
  gd_curv_rearrange_m2d(gdcurv, md);
  

//-------------------------------------------------------------------------------
//-- estimate/check/set time step
//-------------------------------------------------------------------------------
  float   t0 = par->time_start;
  float   dt = par->size_of_time_step;
  int     nt_total = par->number_of_time_steps+1;

  if (par->time_check_stability == 1) {
    fprintf(stdout, "   eatimate time step...\n");
    float dtmax, dtmaxVp, dtmaxL;
    int   dtmaxi, dtmaxk;
    blk_dt_esti(dx, dz, gdinfo, md, fd->CFL, 
        &dtmax, &dtmaxVp, &dtmaxL, &dtmaxi, &dtmaxk);
    //-- print for QC
    fprintf(stdout, "-> dtmax=%f, Vp=%f, i=%d, k=%d\n",
            dtmax, dtmaxVp, dtmaxi, dtmaxk);
    
    // check valid
    if (dtmax <= 0.0) {
       fprintf(stderr,"ERROR: maximum dt <= 0, stop running\n");
       exit(1);
    }

    //-- auto set stept
    if (dt < 0.0) {
       dt       = blk_keep_two_digi(dtmax);
       nt_total = (int) (par->time_window_length / dt + 0.5);

       fprintf(stdout, "-> Set dt       = %f according to maximum allowed value\n", dt);
       fprintf(stdout, "-> Set nt_total = %d\n", nt_total);
    }

    //-- if input dt, check value
    if (dtmax < dt) {
       fprintf(stdout, "Serious Error: dt=%f > dtmax=%f, stop!\n", dt, dtmax);
       exit(1);
    }	
  }

//-------------------------------------------------------------------------------
//-- initial source
//-------------------------------------------------------------------------------
  src_init_bypar(src, par->source_index, par->source_force_vector,
      par->source_peak_time, par->source_freq);

  if (verbose>0) fprintf(stdout,"source info ...\n"); 
  src_print(src, gdinfo, gdcurv);

//-------------------------------------------------------------------------------
//-- allocate main var
//-------------------------------------------------------------------------------

  if (verbose>0) fprintf(stdout,"allocate solver vars ...\n"); 
  wav_init_lebedev(gdinfo, wav);

  
//-------------------------------------------------------------------------------
//-- slover
//-------------------------------------------------------------------------------
  
  // convert rho to 1 / rho to reduce number of arithmetic cal
  md_rho_to_slow(md->rho_lebedev, md->siz_icmp);

  if (verbose>0) fprintf(stdout,"start solver ...\n"); 
  
  time_t t_start = time(NULL);
  
  sv_eliso1st_curv_lebedev_allstep_intv(fd, gdinfo, gdcurv, md, src,
                                        wav, blk->output_dir, 
                                        dx, dz, dt, nt_total, t0);
  
  time_t t_end = time(NULL);
  
  if (verbose>0) {
    fprintf(stdout,"\n\nRuning Time of time :%f s \n", difftime(t_end,t_start));
  }

  free(x1d);
  free(z1d);

  return 0;
}
