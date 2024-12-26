/*
 *
 */

#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "netcdf.h"

#include "constants.h"
#include "fdlib_mem.h"
#include "md_t.h"

int
md_init_lebedev(gdinfo_t *gdinfo, md_t *md)
{
  int ierr = 0;

  md->medium_type = CONST_MEDIUM_ELASTIC;
  md->aniso_type  = CONST_ANISO_ALL; 
  md->visco_type  = CONST_VISCO_NONE; 
  

  md->nx   = gdinfo->nx;
  md->nz   = gdinfo->nz;

  md->siz_iz   = md->nx;
  md->siz_icmp = md->nx * md->nz;

  md->ncmp = 7*2; // cij-6 + rho; for lebedev *2

  /*
   * 0: rho
   * 1: lambda
   * 2: mu
   */
  
  // vars
  md->v4d = (float *) fdlib_mem_calloc_1d_float(
                          md->siz_icmp * md -> ncmp,
                          0.0, "md_init");

  if (md->v4d == NULL) {
      fprintf(stderr,"Error: failed to alloc medium_el_iso\n");
      fflush(stderr);
  }

  // position of each var
  size_t *cmp_pos = (size_t *) fdlib_mem_calloc_1d_sizet(md->ncmp,
                                                         0,
                                                         "medium_init");

  // name of each var
  char **cmp_name = (char **) fdlib_mem_malloc_2l_char(md->ncmp,
                                                       CONST_MAX_STRLEN,
                                                       "medium_init");

  // set pos
  for (int icmp=0; icmp < md->ncmp; icmp++)
  {
    cmp_pos[icmp] = icmp * md->siz_icmp;
  }

  // init
  int icmp = 0;
  sprintf(cmp_name[icmp],"%s","rho");
  md->rho = md->v4d + cmp_pos[icmp];


  // aniso
  icmp += 1;
  sprintf(cmp_name[icmp],"%s","c11");
  md->c11 = md->v4d + cmp_pos[icmp];

  icmp += 1;
  sprintf(cmp_name[icmp],"%s","c13");
  md->c13 = md->v4d + cmp_pos[icmp];

  icmp += 1;
  sprintf(cmp_name[icmp],"%s","c15");
  md->c15 = md->v4d + cmp_pos[icmp];

  icmp += 1;
  sprintf(cmp_name[icmp],"%s","c33");
  md->c33 = md->v4d + cmp_pos[icmp];

  icmp += 1;
  sprintf(cmp_name[icmp],"%s","c35");
  md->c35 = md->v4d + cmp_pos[icmp];

  icmp += 1;
  sprintf(cmp_name[icmp],"%s","c55");
  md->c55 = md->v4d + cmp_pos[icmp];


	// for lebedev
  icmp += 1;
  sprintf(cmp_name[icmp],"%s","rho_lg");
  md->rho_lebedev = md->v4d + cmp_pos[icmp];

  icmp += 1;
  sprintf(cmp_name[icmp],"%s","c11_lg");
  md->c11_lebedev = md->v4d + cmp_pos[icmp];

  icmp += 1;
  sprintf(cmp_name[icmp],"%s","c13_lg");
  md->c13_lebedev = md->v4d + cmp_pos[icmp];

  icmp += 1;
  sprintf(cmp_name[icmp],"%s","c15_lg");
  md->c15_lebedev = md->v4d + cmp_pos[icmp];

  icmp += 1;
  sprintf(cmp_name[icmp],"%s","c33_lg");
  md->c33_lebedev = md->v4d + cmp_pos[icmp];

  icmp += 1;
  sprintf(cmp_name[icmp],"%s","c35_lg");
  md->c35_lebedev = md->v4d + cmp_pos[icmp];

  icmp += 1;
  sprintf(cmp_name[icmp],"%s","c55_lg");
  md->c55_lebedev = md->v4d + cmp_pos[icmp];


  // set pointer
  md->cmp_pos  = cmp_pos;
  md->cmp_name = cmp_name;

  return ierr;
}



int
md_import(md_t *md, char *in_dir)
{
  int ierr = 0;

  char in_file[CONST_MAX_STRLEN];
  
  int ncid, varid;
  
  // construct file name
  sprintf(in_file, "%s/media.nc", in_dir);
  
  // read in nc
  ierr = nc_open(in_file, NC_NOWRITE, &ncid);
  if (ierr != NC_NOERR){
    fprintf(stderr,"nc error: %s\n", nc_strerror(ierr));
    exit(-1);
  }
  
  for (int icmp=0; icmp < md->ncmp; icmp++) {
      ierr = nc_inq_varid(ncid, md->cmp_name[icmp], &varid);
      if (ierr != NC_NOERR){
        fprintf(stderr,"nc error: %s\n", nc_strerror(ierr));
        exit(-1);
      }
  
      ierr = nc_get_var_float(ncid,varid,md->v4d + md->cmp_pos[icmp]);
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

  return ierr;
}

int
md_export(gdinfo_t  *gdinfo,
                 md_t *md,
                 char *output_dir)
{
  int ierr = 0;

  size_t *restrict m3d_pos   = md->cmp_pos;
  char  **restrict m3d_name  = md->cmp_name;
  int  number_of_vars = md->ncmp;
  int  nx = md->nx;
  int  nz = md->nz;
  int  ni1 = gdinfo->ni1;
  int  nk1 = gdinfo->nk1;
  int  ni  = gdinfo->ni;
  int  nk  = gdinfo->nk;

  // construct file name
  char ou_file[CONST_MAX_STRLEN];
  sprintf(ou_file, "%s/media.nc", output_dir);
  
  // read in nc
  int ncid;
  int varid[number_of_vars];
  int dimid[CONST_NDIM];

  ierr = nc_create(ou_file, NC_CLOBBER, &ncid);
  if (ierr != NC_NOERR){
    fprintf(stderr,"creat coord nc error: %s\n", nc_strerror(ierr));
    exit(-1);
  }

  // define dimension
  ierr = nc_def_dim(ncid, "i", nx, &dimid[1]);
  ierr = nc_def_dim(ncid, "k", nz, &dimid[0]);

  // define vars
  for (int ivar=0; ivar<number_of_vars; ivar++) {
    ierr = nc_def_var(ncid, m3d_name[ivar], NC_FLOAT, CONST_NDIM, dimid, &varid[ivar]);
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
  for (int ivar=0; ivar<number_of_vars; ivar++) {
    float *ptr = md->v4d + m3d_pos[ivar];
    ierr = nc_put_var_float(ncid, varid[ivar],ptr);
  }
  
  // close file
  ierr = nc_close(ncid);
  if (ierr != NC_NOERR){
    fprintf(stderr,"nc error: %s\n", nc_strerror(ierr));
    exit(-1);
  }

  return ierr;
}

/*
 * test
 */

int
md_gen_test_ac_iso(md_t *md)
{
  int ierr = 0;

  int nx = md->nx;
  int nz = md->nz;
  int siz_iz = md->siz_iz;

  float *kappa3d = md->kappa;
  float *rho3d = md->rho;

  for (size_t k=0; k<nz; k++)
  {
      for (size_t i=0; i<nx; i++)
      {
        size_t iptr = i + k * siz_iz;
        float Vp=3000.0;
        float rho=1500.0;
        float kappa = Vp*Vp*rho;
        kappa3d[iptr] = kappa;
        rho3d[iptr] = rho;
      }
  }

  return ierr;
}

int
md_gen_test_el_iso(md_t *md)
{
  int ierr = 0;

  int nx = md->nx;
  int nz = md->nz;
  int siz_iz = md->siz_iz;

  float *lam3d = md->lambda;
  float  *mu3d = md->mu;
  float *rho3d = md->rho;

  for (size_t k=0; k<nz; k++)
  {
      for (size_t i=0; i<nx; i++)
      {
        size_t iptr = i + k * siz_iz;
        float Vp=3000.0;
        float Vs=2000.0;
        float rho=1500.0;
        if (k < nz - 30) {
          Vp=5000.0;
          Vs=2500.0;
          rho=2000.0;
        }
        float mu = Vs*Vs*rho;
        float lam = Vp*Vp*rho - 2.0*mu;
        lam3d[iptr] = lam;
         mu3d[iptr] = mu;
        rho3d[iptr] = rho;
      }
  }

  return ierr;
}

int
md_gen_test_Qs(md_t *md, float Qs_freq)
{
  int ierr = 0;

  int nx = md->nx;
  int nz = md->nz;
  int siz_iz = md->siz_iz;

  md->visco_Qs_freq = Qs_freq;

  float *Qs = md->Qs;

  for (size_t k=0; k<nz; k++)
  {
      for (size_t i=0; i<nx; i++)
      {
        size_t iptr = i + k * siz_iz;
        Qs[iptr] = 20;
      }
  }

  return ierr;
}


void
media_read_test_bptti(md_t *md, const char *media_file)
{
	int NX;
	int NZ;
	int nx = md->nx;
	int nz = md->nz;
	int siz_iz = md->siz_iz;

    FILE *fp;
    fp = fopen(media_file, "r");

    fscanf(fp, "%d", &NX);
    fscanf(fp, "%d", &NZ);
    if(NX != nx || NZ != nz) {
        fprintf(stdout, "Error: >>>>  Nx or Nz Mismatch\n");
		fprintf(stdout, "NX = %d, NZ = %d;  nx = %d, nz = %d\n", NX, NZ, nx, nz);
	
        exit(-1);
    }

    for (size_t j=0; j<nz; j++) {
        for (size_t i=0; i<nx; i++) {
      	  size_t iptr = i + j*siz_iz;

            fscanf(fp, "%f", &md->rho[iptr]);
            fscanf(fp, "%f", &md->c11[iptr]);
            fscanf(fp, "%f", &md->c13[iptr]);
            fscanf(fp, "%f", &md->c15[iptr]);
            fscanf(fp, "%f", &md->c33[iptr]);
            fscanf(fp, "%f", &md->c35[iptr]);
            fscanf(fp, "%f", &md->c55[iptr]);
        }
    }
    fclose(fp);
	
}


int
md_gen_test_el_aniso(md_t *md)
{
  int ierr = 0;

  int nx = md->nx;
  int nz = md->nz;
  int siz_iz = md->siz_iz;

  // N0.1 layer
  float rho_1, vp_1, vs_1, epsilon_1, gama_1, delta_1, theta_1;

  float C11_1, C12_1, C13_1, C14_1, C15_1, C16_1;
  float        C22_1, C23_1, C24_1, C25_1, C26_1;
  float               C33_1, C34_1, C35_1, C36_1;
  float                      C44_1, C45_1, C46_1;
  float                             C55_1, C56_1;
  float                                    C66_1;




    rho_1 = 1000.0;
     vp_1 = 3500.0;
     vs_1 = 2400.0;
epsilon_1 = 0.095;
  delta_1 = 0.150;
  theta_1 = PI/6.0;
   gama_1 = 0.0;// no need for 2d

  md_Thomsen_to_Cij_ForTTI(
  		rho_1, vp_1, vs_1, epsilon_1, gama_1, delta_1, theta_1,
  		&C11_1, &C12_1, &C13_1, &C14_1, &C15_1, &C16_1,
  		        &C22_1, &C23_1, &C24_1, &C25_1, &C26_1,
  		                &C33_1, &C34_1, &C35_1, &C36_1,
  		                        &C44_1, &C45_1, &C46_1,
  		                                &C55_1, &C56_1,
  		                                        &C66_1);
  printf("------------------------------\n");
  printf("---------- cij (homo) --------\n");
  printf("------------------------------\n");
  printf("C11 = %f\n", C11_1);
  printf("C13 = %f\n", C13_1);
  printf("C15 = %f\n", C15_1);
  printf("C33 = %f\n", C33_1);
  printf("C35 = %f\n", C35_1);
  printf("C55 = %f\n", C55_1);
  printf("rho = %f\n", rho_1);



  for (size_t k=0; k<nz; k++)
  {
      for (size_t i=0; i<nx; i++)
      {
        size_t iptr = i + k * siz_iz;
		
	      md->c11[iptr] = C11_1; 
	      md->c13[iptr] = C13_1; 
	      md->c15[iptr] = C15_1; 
	      md->c33[iptr] = C33_1; 
	      md->c35[iptr] = C35_1; 
	      md->c55[iptr] = C55_1; 
          md->rho[iptr] = rho_1;
      }
  }

  return ierr;
}

/*
 * convert rho to slowness to reduce number of arithmetic cal
 */

int
md_rho_to_slow(float *restrict rho, size_t siz_volume)
{
  int ierr = 0;

  for (size_t iptr=0; iptr<siz_volume; iptr++) {
    if (rho[iptr] > 1e-10) {
      rho[iptr] = 1.0 / rho[iptr];
    } else {
      rho[iptr] = 0.0;
    }
  }

  return ierr;
}

int
md_is_el_iso(md_t *md)
{
  int is = 0;
  if (md->medium_type == CONST_MEDIUM_ELASTIC &&
      md->aniso_type == CONST_ANISO_ISO)
  {
    is = 1;
  }

  return is;
}

int
md_is_el_aniso(md_t *md)
{
  int is = 0;
  if (md->medium_type == CONST_MEDIUM_ELASTIC &&
      md->aniso_type == CONST_ANISO_ALL)
  {
    is = 1;
  }

  return is;
}

int
md_is_ac_iso(md_t *md)
{
  int is = 0;
  if (md->medium_type == CONST_MEDIUM_ACOUSTIC &&
      md->aniso_type == CONST_ANISO_ISO)
  {
    is = 1;
  }

  return is;
}


void 
md_Thomsen_to_Cij_ForTTI(
		float rho, float vp, float vs, float sigma, float gama, float delta, float theta,
	    float *C11, float *C12, float *C13, float *C14, float *C15, float *C16,
	                float *C22, float *C23, float *C24, float *C25, float *C26,
	                            float *C33, float *C34, float *C35, float *C36,
	                                        float *C44, float *C45, float *C46,
	                                                    float *C55, float *C56,
	                                                                float *C66)

{
	float C_33, C_44, C_11, C_66, C_13, C_12;

    C_33 = vp*vp*rho;
    C_44 = vs*vs*rho;
    C_11 = 2*C_33*sigma + C_33;
    C_66 = 2*C_44*gama + C_44;
    C_13 = sqrt(( 2*C_33*C_33*delta + (C_33-C_44)*(C_11+C_33-2*C_44) ) * 0.5) - C_44;
    C_12 = C_11-2*C_66;

    *C11 = ( cos(theta)*cos(theta)*C_11 + sin(theta)*sin(theta)*C_13 ) * cos(theta)*cos(theta) + ( cos(theta)*cos(theta)*C_13 + sin(theta)*sin(theta)*C_33 ) * sin(theta)*sin(theta)
        + sin(2*theta)*sin(2*theta)*C_44;
 
    *C12 = cos(theta)*cos(theta)*C_12 + sin(theta)*sin(theta)*C_13;

    *C13 = ( cos(theta)*cos(theta)*C_11 + sin(theta)*sin(theta)*C_13 ) * sin(theta)*sin(theta)
        + ( cos(theta)*cos(theta)*C_13 + sin(theta)*sin(theta)*C_33 ) * cos(theta)*cos(theta)
        - sin(2*theta)*sin(2*theta)*C_44;

	*C14 = 0.0;

    *C15 = 0.5*(cos(theta)*cos(theta)*C_11 + sin(theta)*sin(theta)*C_13) * sin(2*theta)
        - 0.5*(cos(theta)*cos(theta)*C_13 + sin(theta)*sin(theta)*C_33) * sin(2*theta)
        - sin(2*theta)*C_44 * cos(2*theta);

	*C16 = 0.0;


    *C22 = *C11;
    
    *C23 = sin(theta)*sin(theta)*C_12 + cos(theta)*cos(theta)*C_13;

	*C24 = 0.0;

    *C25 = 0.5 * sin(2*theta)*C_12 - 0.5*sin(2*theta)*C_13;

	*C26 = 0.0;


    *C33 = (sin(theta)*sin(theta)*C_11 + cos(theta)*cos(theta)*C_13) * sin(theta)*sin(theta)
        + (sin(theta)*sin(theta)*C_13 + cos(theta)*cos(theta)*C_33) * cos(theta)*cos(theta)
        + sin(2*theta)*sin(2*theta)*C_44;
	 
	*C34 = 0.0;

    *C35 = 0.5 * (sin(theta)*sin(theta)*C_11 + cos(theta)*cos(theta)*C_13) * sin(2*theta)
        - 0.5 * (sin(theta)*sin(theta)*C_13 + cos(theta)*cos(theta)*C_33) * sin(2*theta)
        + sin(2*theta)*C_44*cos(2*theta);
	
	*C36 = 0.0;
	

    *C44 = cos(theta)*cos(theta)*C_44 + sin(theta)*sin(theta)*C_66;
	 
	*C45 = 0.0;

    *C46 = -cos(theta)*sin(theta)*C_44 + cos(theta)*sin(theta)*C_66;


    *C55 = 0.25*(sin(2*theta)*C_11 - sin(2*theta)*C_13) * sin(2*theta)
        - 0.25*(sin(2*theta)*C_13 - sin(2*theta)*C_33) * sin(2*theta)
        + pow(cos(2*theta), 2) * C_44;

	*C56 = 0.0;

    *C66 = sin(theta)*sin(theta)*C_44 + cos(theta)*cos(theta)*C_66;
}

