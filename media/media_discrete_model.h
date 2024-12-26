#ifndef MEDIA_DISCRETE_MODEL_H
#define MEDIA_DISCRETE_MODEL_H

// for C code call
#define MEDIA_USE_CART 1
#define MEDIA_USE_VMAP 2
#define MEDIA_USE_CURV 3

/*--------------------------- layer2model --------------------- */

//--- 4. elastic anisotropic/TTI
int media_layer2model_el_aniso(
        float *rho,
        float *c11, float *c13, float *c15,
        float *c33, float *c35, 
        float *c55,
        const float *x2d,
        const float *z2d,
        size_t nx,
        size_t nz,
        int grid_type, 
        const char *in_lay_file,
        const char *equivalent_medium_method) ; 


#endif
