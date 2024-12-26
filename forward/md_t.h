#ifndef MD_EL_ISO_H
#define MD_EL_ISO_H

#include "gd_info.h"

/*************************************************
 * structure
 *************************************************/

typedef struct {
  int n1, n2, n3, n4;
  int nx, nz, ncmp;
  float *v4d; // allocated var

  size_t siz_iz;
  size_t siz_icmp;

  size_t *cmp_pos;
  char  **cmp_name;

  // flag to determine medium type
  int medium_type;
  int aniso_type;
  int visco_type;

  // rho for all media
  float *rho;

  // for acustic
  float *kappa; // pointer to var

  // for isotropic media
  float *lambda; // pointer to var
  float *mu;

  // for visco attenuation
  float *Qs;

  // for anisotropic media
  float *c11;
  float *c13;
  float *c15;
  float *c33;
  float *c35;
  float *c55;

  // for lebedev
  float *c11_lebedev;
  float *c13_lebedev;
  float *c15_lebedev;
  float *c33_lebedev;
  float *c35_lebedev;
  float *c55_lebedev;
  float *rho_lebedev;

  float visco_Qs_freq;

} md_t;

/*************************************************
 * function prototype
 *************************************************/


int
md_import(md_t *md, char *in_dir);

int
md_export(gdinfo_t  *gdinfo,
                 md_t *md,
                 char *output_dir);

int
md_gen_test_ac_iso(md_t *md);

int
md_gen_test_el_iso(md_t *md);

int
md_gen_test_Qs(md_t *md, float Qs_freq);

void
media_read_test_bptti(md_t *md, const char *media_file);

int
md_gen_test_el_aniso(md_t *md);

int
md_rho_to_slow(float *restrict rho, size_t siz_volume);

int
md_is_el_iso(md_t *md);

int
md_is_el_aniso(md_t *md);

int
md_is_ac_iso(md_t *md);

void 
md_Thomsen_to_Cij_ForTTI(
		float rho, float vp, float vs, float sigma, float gama, float delta, float theta,
	    float *C11, float *C12, float *C13, float *C14, float *C15, float *C16,
	                float *C22, float *C23, float *C24, float *C25, float *C26,
	                            float *C33, float *C34, float *C35, float *C36,
	                                        float *C44, float *C45, float *C46,
	                                                    float *C55, float *C56,
	                                                                float *C66);

// for lebedev
int md_init_lebedev(gdinfo_t *gdinfo, md_t *md);

#endif
