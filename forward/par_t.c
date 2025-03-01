/*
 * 
 */

//#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "par_t.h"

/*
 * read from file
 */

void
par_read_from_file(char *par_fname,  par_t *par, int verbose)
{
  //
  // read whole file inot str
  //
  FILE *fp = fopen(par_fname,"r");
  if (!fp) {
    fprintf(stderr,"Error: can't open par file: %s\n", par_fname);
    exit(1);
  }

  fseek(fp, 0, SEEK_END);
  long len = ftell(fp);

  fseek(fp, 0, SEEK_SET);
  char *str = (char*)malloc(len+1);
  fread(str, 1, len, fp);
  fclose(fp);

  // read from str
  par_read_from_str(str, par);

  return;
}

/*
 * funcs to get par from alread read in str
 */
int 
par_read_from_str(const char *str, par_t *par)
{
  int ierr = 0;

  // allocate
  par->boundary_type_name = (char **)malloc(CONST_NDIM_2 * sizeof(char*));
  for (int i=0; i<CONST_NDIM_2; i++) {
    par->boundary_type_name[i] = (char *)malloc(10*sizeof(char));
  }

  // convert str to json
  cJSON *root = cJSON_Parse(str);
  if (NULL == root) {
    printf("Error at parsing json!\n");
    exit(1);
  }

  cJSON *item;
  cJSON *subitem, *thirditem, *snapitem, *lineitem;

  // no default
  if (item = cJSON_GetObjectItem(root, "number_of_total_grid_points_x")) {
    par->number_of_total_grid_points_x = item->valueint;
  }
  if (item = cJSON_GetObjectItem(root, "number_of_total_grid_points_z")) {
    par->number_of_total_grid_points_z = item->valueint;
  }

  // dis grid
  par->disg_num_level = 0;
  if (item = cJSON_GetObjectItem(root, "dis_grid_at_zindx")) {
    par->disg_num_level = cJSON_GetArraySize(item);
    par->disg_at_zindx = (int *)malloc(par->disg_num_level * sizeof(int));
    for (int n=0; n < par->disg_num_level; n++) {
      par->disg_at_zindx[n] = cJSON_GetArrayItem(item, n)->valueint;
    }
  }
  if (item = cJSON_GetObjectItem(root, "dis_grid_factor")) {
    // should be improved to allow input order change
    if (par->disg_num_level != cJSON_GetArraySize(item)) {
      fprintf(stderr,"ERROR: input size of dis_grid_at_zindx and dis_grid_factor diff\n");
      exit(1);
    }
    par->disg_factor = (int *)malloc(par->disg_num_level * sizeof(int));
    for (int n=0; n < par->disg_num_level; n++) {
      par->disg_factor[n] = cJSON_GetArrayItem(item, n)->valueint;
    }
  }

  // set default values to negative
  par->size_of_time_step = -1.0;
  par->number_of_time_steps = -1;
  par->time_window_length = -1.0;
  par->time_start = 0.0;
  par->time_check_stability = 1;

  if (item = cJSON_GetObjectItem(root, "size_of_time_step")) {
    par->size_of_time_step = item->valuedouble;
  }
  if (item = cJSON_GetObjectItem(root, "number_of_time_steps")) {
    par->number_of_time_steps = item->valueint;
  }
  if (item = cJSON_GetObjectItem(root, "time_window_length")) {
    par->time_window_length = item->valuedouble;
  }
  if (item = cJSON_GetObjectItem(root, "time_start")) {
    par->time_start = item->valuedouble;
  }
  if (item = cJSON_GetObjectItem(root, "time_check_stability")) {
    par->time_check_stability = item->valueint;
  }

  // check 
  //int num_time_minus = 0;
  //if (par->size_of_time_step    < 0.0) num_time_minus += 1;
  //if (par->time_window_length   < 0.0) num_time_minus += 1;
  //if (par->number_of_time_steps < 0  ) num_time_minus += 1;
  //if (num_time_minus >= 2)
  //{
  //  fprintf(stderr," --> size_of_time_step   =%f\n", par->size_of_time_step);
  //  fprintf(stderr," --> number_of_time_steps=%d\n", par->number_of_time_steps);
  //  fprintf(stderr," --> time_window_length  =%f\n", par->time_window_length);
  //  fprintf(stderr,"Error: at lest two of above three paras should > 0\n");
  //  exit(-1);
  //}

  if (par->size_of_time_step < 0.0 && par->time_window_length < 0)
  {
    fprintf(stderr," --> size_of_time_step   =%f\n", par->size_of_time_step);
    fprintf(stderr," --> time_window_length  =%f\n", par->time_window_length);
    fprintf(stderr,"Error: at lest one of above paras should > 0\n");
    exit(1);
  }

  if (par->size_of_time_step > 0.0)
  {
    if (par->number_of_time_steps < 0 && par->time_window_length < 0.0)
    {
      fprintf(stderr,"Error: both size_of_time_step=%f, time_window_length=%f less 0\n",
             par->size_of_time_step, par->time_window_length);
      exit(1);
    }

    if (par->number_of_time_steps < 0) 
    {
      par->number_of_time_steps = (int)(par->time_window_length / par->size_of_time_step + 0.5);
    }

    par->time_window_length = par->size_of_time_step * par->number_of_time_steps;
  }

  //par->time_end   = par->time_start + 
  //        par->number_of_time_steps * par->size_of_time_step;
  //int     nt_total = (int) ((par->time_end - par->time_start) / dt+0.5);

  //
  // boundary default values
  //
  for (int idim=0; idim < CONST_NDIM; idim++) {
    for (int iside=0; iside < 2; iside++) {
      par->abs_num_of_layers[idim][iside] = 0;
      par->cfspml_is_sides[idim][iside] = 0;
      par->free_is_sides  [idim][iside] = 0;
    }
  }
  par->bdry_has_cfspml = 0;
  par->bdry_has_free   = 0;

  //
  //-- grid
  //

  // default output grid

  par->grid_generation_itype = PAR_GRID_IMPORT;
  if (item = cJSON_GetObjectItem(root, "grid_generation_method")) {
    // import grid
    if (subitem = cJSON_GetObjectItem(item, "import")) {
       par->grid_generation_itype = PAR_GRID_IMPORT;
       sprintf(par->grid_import_dir, "%s", subitem->valuestring);
    }
    // generate cartesian grid
    if (subitem = cJSON_GetObjectItem(item, "cartesian")) {
       par->grid_generation_itype = PAR_GRID_CARTESIAN;
       if (thirditem = cJSON_GetObjectItem(subitem, "origin")) {
         for (int i = 0; i < CONST_NDIM; i++) {
           par->cartesian_grid_origin[i] = cJSON_GetArrayItem(thirditem, i)->valuedouble;
         }
       }
       if (thirditem = cJSON_GetObjectItem(subitem, "inteval")) {
         for (int i = 0; i < CONST_NDIM; i++) {
           par->cartesian_grid_stepsize[i] = cJSON_GetArrayItem(thirditem, i)->valuedouble;
         }
       }
    }
    // layer interp
    if (subitem = cJSON_GetObjectItem(item, "layer_interp")) {
       par->grid_generation_itype = PAR_GRID_LAYER_INTERP;
       if (thirditem = cJSON_GetObjectItem(subitem, "in_grid_layer_file")) {
          sprintf(par->in_grid_layer_file, "%s", thirditem->valuestring);
       }
       if (thirditem = cJSON_GetObjectItem(subitem, "refine_factor")) {
         for (int i = 0; i < CONST_NDIM; i++) {
           par->grid_layer_resample_factor[i] = cJSON_GetArrayItem(thirditem, i)->valueint;
         }
       }
       if (thirditem = cJSON_GetObjectItem(subitem, "horizontal_start_index")) {
         for (int i = 0; i < CONST_NDIM-1; i++) {
           par->grid_layer_start[i] = cJSON_GetArrayItem(thirditem, i)->valueint;
         }
       }
       if (thirditem = cJSON_GetObjectItem(subitem, "vertical_ToFreeSurf_resample_index")) {
         par->grid_layer_start[CONST_NDIM-1] = thirditem->valueint;
       }
    }
  }

  par->is_export_grid = 1;
  if (item = cJSON_GetObjectItem(root, "is_export_grid")) {
     par->is_export_grid = item->valueint;
  }
  if (item = cJSON_GetObjectItem(root, "grid_export_dir")) {
      sprintf(par->grid_export_dir,"%s",item->valuestring);
  }

  //
  //-- metric
  //

  par->metric_method_itype = PAR_METRIC_CALCULATE;
  if (item = cJSON_GetObjectItem(root, "metric_calculation_method")) {
    if (subitem = cJSON_GetObjectItem(item, "import")) {
        par->metric_method_itype = PAR_METRIC_IMPORT;
        sprintf(par->metric_import_dir, "%s", subitem->valuestring);
    }
    if (subitem = cJSON_GetObjectItem(item, "calculate")) {
        par->metric_method_itype = PAR_METRIC_CALCULATE;
    }
  }

  par->is_export_metric = 1;
  if (item = cJSON_GetObjectItem(root, "is_export_metric")) {
     par->is_export_metric = item->valueint;
  }

  //
  //-- medium
  //

  par->media_input_itype = PAR_MEDIA_IMPORT;
  if (item = cJSON_GetObjectItem(root, "medium")) {
    if (subitem = cJSON_GetObjectItem(item, "input_way")) {
        sprintf(par->media_input_way, "%s", subitem->valuestring);
    }
    // if input by import
    if (strcmp(par->media_input_way,"import") == 0) {
        par->media_input_itype = PAR_MEDIA_IMPORT;
        sprintf(par->media_import_dir, "%s", subitem->valuestring);
    }
    // if input by generate in side
    if (strcmp(par->media_input_way,"infile_layer") == 0) {
        par->media_input_itype = PAR_MEDIA_LAY;
        if (subitem = cJSON_GetObjectItem(item, "input_file")) {
          sprintf(par->media_input_file, "%s", subitem->valuestring);
        } else {
          fprintf(stderr, "Please give the medium file! \n");
          exit(1);
        }
        if (subitem = cJSON_GetObjectItem(item, "equivalent_medium_method")) {
          sprintf(par->equivalent_medium_method, "%s", subitem->valuestring);
        } else {
          fprintf(stderr, "Please give the type of equivalent medium methods\n");
          exit(1);
        }
    }
    if (strcmp(par->media_input_way,"infile_ascii") == 0) {
        par->media_input_itype = PAR_MEDIA_ASCII;
        if (subitem = cJSON_GetObjectItem(item, "input_file")) {
          sprintf(par->media_input_file, "%s", subitem->valuestring);
        } else {
          fprintf(stderr, "Please give the medium file! \n");
          exit(1);
        }
    }
  }

  par->is_export_media = 1;
  if (item = cJSON_GetObjectItem(root, "is_export_media")) {
     par->is_export_media = item->valueint;
  }
  if (item = cJSON_GetObjectItem(root, "media_export_dir")) {
      sprintf(par->media_export_dir,"%s",item->valuestring);
  }

  // filtering, for fsg
  par->is_filter = 1;
  par->filter_sigma = 0.2;
  if (item = cJSON_GetObjectItem(root, "is_filter")) {
     par->is_filter = item->valueint;
  }
  if (item = cJSON_GetObjectItem(root, "filter_sigma")) {
     par->filter_sigma = item->valuedouble;
  }
  

  //
  //-- source
  //

  if (item = cJSON_GetObjectItem(root, "source"))
  {
    if (subitem = cJSON_GetObjectItem(item, "index")) {
      for (int i = 0; i < CONST_NDIM; i++) {
        int indx = cJSON_GetArrayItem(subitem, i)->valueint;
        par->source_index[i] = indx;
      }
    }
    if (subitem = cJSON_GetObjectItem(item, "force_vector")) {
      for (int i = 0; i < 2; i++) {
        float F = cJSON_GetArrayItem(subitem, i)->valuedouble;
        par->source_force_vector[i] = F;
      }
    }
    if (subitem = cJSON_GetObjectItem(item, "ricker_center_frequency")) {
      par->source_freq = subitem->valuedouble;
    }
    if (subitem = cJSON_GetObjectItem(item, "ricker_peak_time")) {
      par->source_peak_time = subitem->valuedouble;
    }
  }

  //-- output dir
  if (item = cJSON_GetObjectItem(root, "output_dir")) {
      sprintf(par->output_dir,"%s",item->valuestring);
  }

  //-- receiver
  if (item = cJSON_GetObjectItem(root, "in_station_file")) {
    sprintf(par->in_station_file, "%s", item->valuestring);
  }

  //-- receiver line
  if (item = cJSON_GetObjectItem(root, "receiver_line"))
  {
    par->number_of_receiver_line = cJSON_GetArraySize(item);
    par->receiver_line_index_start  = (int *)malloc(par->number_of_receiver_line*sizeof(int)*CONST_NDIM);
    par->receiver_line_index_incre  = (int *)malloc(par->number_of_receiver_line*sizeof(int)*CONST_NDIM);
    par->receiver_line_count  = (int *)malloc(par->number_of_receiver_line*sizeof(int));
    //par->receiver_line_time_interval  = (int *)malloc(par->number_of_receiver_line*sizeof(int));
    par->receiver_line_name = (char **)malloc(par->number_of_receiver_line*sizeof(char*));
    for (int n=0; n<par->number_of_receiver_line; n++) {
      par->receiver_line_name[n] = (char *)malloc(PAR_MAX_STRLEN*sizeof(char));
    }
    // each line
    for (int i=0; i < cJSON_GetArraySize(item) ; i++)
    {
      lineitem = cJSON_GetArrayItem(item, i);

      if (subitem = cJSON_GetObjectItem(lineitem, "name"))
      {
        sprintf(par->receiver_line_name[i],"%s",subitem->valuestring);
      }

      if (subitem = cJSON_GetObjectItem(lineitem, "grid_index_start"))
      {
        for (int j = 0; j < CONST_NDIM; j++) {
          par->receiver_line_index_start[i*CONST_NDIM+j] = cJSON_GetArrayItem(subitem, j)->valueint;
        }
      }

      if (subitem = cJSON_GetObjectItem(lineitem, "grid_index_incre"))
      {
        for (int j = 0; j < CONST_NDIM; j++) {
          par->receiver_line_index_incre[i*CONST_NDIM+j] = cJSON_GetArrayItem(subitem, j)->valueint;
        }
      }

      if (subitem = cJSON_GetObjectItem(lineitem, "grid_index_count"))
      {
         par->receiver_line_count[i] = subitem->valueint;
      }

      //if (subitem = cJSON_GetObjectItem(lineitem, "t_index_interval"))
      //{
      //   par->receiver_line_tinterval[i] = cJSON_GetArrayItem(subitem, j)->valueint;
      //}
    }
  }

  // snapshot
  if (item = cJSON_GetObjectItem(root, "snapshot"))
  {
    par->number_of_snapshot = cJSON_GetArraySize(item);
    //fprintf(stdout,"size=%d, %d, %d\n", par->number_of_snapshot, sizeof(int), CONST_NDIM);
    //fflush(stdout);
    par->snapshot_index_start  = (int *)malloc(par->number_of_snapshot*sizeof(int)*CONST_NDIM);
    //if (par->snapshot_index_start == NULL) {
    //  fprintf(stdout,"eror\n");
    //  fflush(stdout);
    //}
    par->snapshot_index_count  = (int *)malloc(par->number_of_snapshot*sizeof(int)*CONST_NDIM);
    par->snapshot_index_incre = (int *)malloc(par->number_of_snapshot*sizeof(int)*CONST_NDIM);
    par->snapshot_time_start  = (int *)malloc(par->number_of_snapshot*sizeof(int));
    par->snapshot_time_incre = (int *)malloc(par->number_of_snapshot*sizeof(int));
    par->snapshot_save_velocity = (int *)malloc(par->number_of_snapshot*sizeof(int));
    par->snapshot_save_stress  = (int *)malloc(par->number_of_snapshot*sizeof(int));
    par->snapshot_save_strain = (int *)malloc(par->number_of_snapshot*sizeof(int));
    // name of snapshot
    par->snapshot_name = (char **)malloc(par->number_of_snapshot*sizeof(char*));
    for (int n=0; n<par->number_of_snapshot; n++) {
      par->snapshot_name[n] = (char *)malloc(PAR_MAX_STRLEN*sizeof(char));
    }

    // each snapshot
    for (int i=0; i < cJSON_GetArraySize(item) ; i++)
    {
      snapitem = cJSON_GetArrayItem(item, i);

      if (subitem = cJSON_GetObjectItem(snapitem, "name"))
      {
        sprintf(par->snapshot_name[i],"%s",subitem->valuestring);
      }

      if (subitem = cJSON_GetObjectItem(snapitem, "grid_index_start"))
      {
        for (int j = 0; j < CONST_NDIM; j++) {
          par->snapshot_index_start[i*CONST_NDIM+j] = cJSON_GetArrayItem(subitem, j)->valueint;
        }
      }

      if (subitem = cJSON_GetObjectItem(snapitem, "grid_index_count"))
      {
        for (int j = 0; j < CONST_NDIM; j++) {
          par->snapshot_index_count[i*CONST_NDIM+j] = cJSON_GetArrayItem(subitem, j)->valueint;
        }
      }

      if (subitem = cJSON_GetObjectItem(snapitem, "grid_index_incre"))
      {
        for (int j = 0; j < CONST_NDIM; j++) {
          par->snapshot_index_incre[i*CONST_NDIM+j] = cJSON_GetArrayItem(subitem, j)->valueint;
        }
      }

      if (subitem = cJSON_GetObjectItem(snapitem, "time_index_start")) {
        par->snapshot_time_start[i]  = subitem->valueint;
      }
      if (subitem = cJSON_GetObjectItem(snapitem, "time_index_incre")) {
        par->snapshot_time_incre[i]  = subitem->valueint;
      }
      if (subitem = cJSON_GetObjectItem(snapitem, "save_velocity")) {
        par->snapshot_save_velocity[i]  = subitem->valueint;
      }
      if (subitem = cJSON_GetObjectItem(snapitem, "save_stress")) {
        par->snapshot_save_stress[i]  = subitem->valueint;
      }
      if (subitem = cJSON_GetObjectItem(snapitem, "save_strain")) {
        par->snapshot_save_strain[i] = subitem->valueint;
      }
    }
  }

  //-- misc
  if (item = cJSON_GetObjectItem(root, "check_nan_every_nummber_of_steps")) {
      par->check_nan_every_nummber_of_steps = item->valueint;
  }
  if (item = cJSON_GetObjectItem(root, "output_all")) {
      par->output_all = item->valueint;
  }

  //if (item = cJSON_GetObjectItem(root, "grid_name")) {
  //    sprintf(par->grid_name,"%s",item->valuestring);
  //}

  cJSON_Delete(root);

  // set values to default ones if no input

  return ierr;
}

/*
 * funcs to read cfspml para from json str
 */

void 
par_read_json_cfspml(cJSON *item,
      int *nlay, float *amax, float *bmax, float *vel)
{
  cJSON *subitem;

  if (subitem = cJSON_GetObjectItem(item, "number_of_layers"))
  {
    *nlay = subitem->valueint;
  }
  if (subitem = cJSON_GetObjectItem(item, "alpha_max"))
  {
    *amax = subitem->valuedouble;
  }
  if (subitem = cJSON_GetObjectItem(item, "beta_max"))
  {
    *bmax = subitem->valuedouble;
  }
  if (subitem = cJSON_GetObjectItem(item, "ref_vel"))
  {
    *vel = subitem->valuedouble;
  }
}

/*
 * funcs to read index/wavelet para from json str
 */
void 
par_read_json_source(cJSON *item,
      float *src_coord, int *grid_index, float *grid_inc,
      float *force_vector,  int *force_actived,
      float *moment_tensor, int *moment_actived,
      char *wavelet_name, float *wavelet_coefs, float *t_start, float *t_end)
{
  cJSON *subitem;

  // default values
  for (int idim=0; idim < CONST_NDIM; idim++) {
    grid_index[idim] = -1; // compare -1 if coord input
    grid_inc[idim]   = 0.0;
    src_coord[idim]  = 0.0;
    force_vector[idim]  = 0.0;
  }
  for (int idim=0; idim < 3; idim++) {
    moment_tensor[idim  ]  = 0.0;
  }

  *force_actived = 0;
  *moment_actived = 0;

  // if coordinate
  if (subitem = cJSON_GetObjectItem(item, "coord")) {
    for (int i = 0; i < CONST_NDIM; i++) {
      src_coord[i] = cJSON_GetArrayItem(subitem, i)->valuedouble;
    }
  }

  // if grid index
  if (subitem = cJSON_GetObjectItem(item, "index")) {
    for (int i = 0; i < CONST_NDIM; i++) {
      float indx_w_inc = cJSON_GetArrayItem(subitem, i)->valuedouble;
      grid_index[i] = (int) (indx_w_inc + 0.5);
      grid_inc  [i] = indx_w_inc - grid_index[i];
    }
  }

  // stf acts from start time
  if (subitem = cJSON_GetObjectItem(item, "start_time")) {
     *t_start = subitem->valuedouble;
  }
  // stf ends at end time
  if (subitem = cJSON_GetObjectItem(item, "end_time")) {
     *t_end = subitem->valuedouble;
  }

  // if force vector
  if (subitem = cJSON_GetObjectItem(item, "force_vector")) {
    *force_actived = 1;
    for (int i = 0; i < CONST_NDIM; i++) {
      force_vector[i] = cJSON_GetArrayItem(subitem, i)->valuedouble;
    }
  }

  // if moment vector
  if (subitem = cJSON_GetObjectItem(item, "moment_tensor")) {
    *moment_actived = 1;
    for (int i = 0; i < 3; i++) {
      moment_tensor[i] = cJSON_GetArrayItem(subitem, i)->valuedouble;
    }
  }

  // wavelet name
  if (subitem = cJSON_GetObjectItem(item, "wavelet_name"))
  {
    sprintf(wavelet_name,"%s",subitem->valuestring);
  }

  // coefs

  // ricker
  if (strcmp(wavelet_name, "ricker")==0 || strcmp(wavelet_name, "ricker_deriv")==0) {
    if (subitem = cJSON_GetObjectItem(item, "ricker_center_frequency")) {
      wavelet_coefs[0] = subitem->valuedouble;
    }
    if (subitem = cJSON_GetObjectItem(item, "ricker_peak_time")) {
      wavelet_coefs[1] = subitem->valuedouble;
    }
  }

  // gaussian
  if (strcmp(wavelet_name, "gaussian")==0 || strcmp(wavelet_name, "gaussian_deriv")==0) {
    if (subitem = cJSON_GetObjectItem(item, "gaussian_rms_width")) {
      wavelet_coefs[0] = subitem->valuedouble;
    }
    if (subitem = cJSON_GetObjectItem(item, "gaussian_peak_time"))
    {
      wavelet_coefs[1] = subitem->valuedouble;
    }
  }

  return;
}

int
par_print(par_t *par)
{    
  int ierr = 0;

  fprintf(stdout, "-------------------------------------------------------\n");

  fprintf(stdout, "-------------------------------------------------------\n");
  fprintf(stdout, "--> Time Integration information:\n");
  fprintf(stdout, "-------------------------------------------------------\n");
  fprintf(stdout, " size_of_time_step = %10.4e\n", par->size_of_time_step);
  fprintf(stdout, " number_of_time_steps = %-10d\n", par->number_of_time_steps);

  fprintf(stdout, "-------------------------------------------------------\n");
  fprintf(stdout, "--> GRID information:\n");
  fprintf(stdout, "-------------------------------------------------------\n");
  fprintf(stdout, " grid_export_dir = %s\n", par->grid_export_dir);
  fprintf(stdout, " number_of_total_grid_points_x = %-10d\n", par->number_of_total_grid_points_x);
  fprintf(stdout, " number_of_total_grid_points_z = %-10d\n", par->number_of_total_grid_points_z);

  fprintf(stdout, " disg_num_level = %-10d\n", par->disg_num_level);
  for (int n; n < par->disg_num_level; n++) {
    fprintf(stdout, "    #%d: at %d, factor=%d\n",
          n, par->disg_at_zindx[n], par->disg_factor[n]);
  }

  fprintf(stdout, " grid_generation_itype = %d\n", par->grid_generation_itype);
  if (par->grid_generation_itype==PAR_GRID_CARTESIAN) {
    fprintf(stdout, " cartesian_grid_x0 = %10.4e\n", par->cartesian_grid_origin[0]);
    fprintf(stdout, " cartesian_grid_z0 = %10.4e\n", par->cartesian_grid_origin[2]);
    fprintf(stdout, " cartesian_grid_dx = %10.4e\n", par->cartesian_grid_stepsize[0]);
    fprintf(stdout, " cartesian_grid_dz = %10.4e\n", par->cartesian_grid_stepsize[2]);
  }

  fprintf(stdout, " metric_method_itype = %d\n", par->metric_method_itype);

  fprintf(stdout, "-------------------------------------------------------\n");
  fprintf(stdout, "--> media info.\n");
  fprintf(stdout, "-------------------------------------------------------\n");
  fprintf(stdout, " media_type = %s\n", par->media_type);
  fprintf(stdout, " media_export_dir = %s\n", par->media_export_dir);
  fprintf(stdout, " media_input_itype = %d\n", par->media_input_itype);


  //if (par->media_input_itype == PAR_MEDIA_3LAY) {
  //  fprintf()
  //}
  //fprintf(stdout, "\n --> the media filename is:\n");
  //fprintf(stdout, " velp_file  = %s\n", PSV->fnm_velp);
  //fprintf(stdout, " vels_file  = %s\n", PSV->fnm_vels);
  //fprintf(stdout, " rho_file   = %s\n", PSV->fnm_rho);

  fprintf(stdout, "-------------------------------------------------------\n");
  fprintf(stdout, "--> source info.\n");
  fprintf(stdout, "-------------------------------------------------------\n");
  // grid index 
  fprintf(stdout, " source index = [%d,  %d]\n", par->source_index[0], par->source_index[1]);
  fprintf(stdout, " source freq = %g\n", par->source_freq);
  fprintf(stdout, " source peak_time = %g\n", par->source_peak_time);
  fprintf(stdout, " force vector = [%g, %g]\n", par->source_force_vector[0],
               par->source_force_vector[1]);

  fprintf(stdout, "\n");
  fprintf(stdout, "-------------------------------------------------------\n");
  fprintf(stdout, "--> output information.\n");
  fprintf(stdout, "-------------------------------------------------------\n");
  fprintf(stdout, "--> output_dir = %s\n", par->output_dir);

  fprintf(stdout, "--> snapshot information.\n");
  fprintf(stdout, "number_of_snapshot=%d\n", par->number_of_snapshot);
  if (par->number_of_snapshot > 0)
  {
      fprintf(stdout, "#  name  i0  k0  ni  nk  di  dk  it0  nt  dit\n");
      for(int n=0; n<par->number_of_snapshot; n++)
      {
         fprintf(stdout, "%6d %s %6d %6d %6d %6d %6d %6d %6d %6d\n",
             n,
             par->snapshot_name[n],
             par->snapshot_index_start[n*2+0],
             par->snapshot_index_start[n*2+1],
             par->snapshot_index_count[n*2+0],
             par->snapshot_index_count[n*2+1],
             par->snapshot_index_incre[n*2+0],
             par->snapshot_index_incre[n*2+1],
             par->snapshot_time_start[n],
             par->snapshot_time_incre[n]);
      }
  }

  fprintf(stdout, "--> qc parameters:\n");
  fprintf(stdout, "check_nan_every_nummber_of_steps=%d\n", par->check_nan_every_nummber_of_steps);
  fprintf(stdout, "output_all=%d\n", par->output_all);

  return ierr;
}
