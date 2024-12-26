#ifndef _MEDIA_LAYER2MODEL_
#define _MEDIA_LAYER2MODEL_

#include "media_geometry2d.hpp"
#include "media_utility.hpp"


#ifdef __cplusplus
extern "C" {
#endif

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
        const char *equivalent_medium_method); 

#ifdef __cplusplus
}
#endif


int AssignLayerMediaPara2Point(
    size_t ix, size_t iz,         /* To print Error messages */ 
    Point2 A,  
    inter_t *interfaces,
    int media_type,                /* the type can be found in media_utility.hpp */ 
    std::vector<float> &var); 

//- Calculate the value of the point for different media type (to avoid multiple geometric calculations) 
//   for layer2model
void CalPointValue_layer(int media_type, 
                   inter_t *interfaces,
                   Point2 &A,
                   std::vector<float> &elevation,  /*the elevation of point A at the projection position of the interface mesh. */
                   std::vector<int> &internum4elev,
                   int mi,
                   std::vector<float> &var);

void MarkInterfaceNumber(
        int grid_type,
        float *Hx, float *Hz,
        size_t nx, size_t nz,
        int *MaterNum, // nx*ny*nz
        inter_t *interfaces);


//- 4. Assign the parameter directly (use the local values): elastic tti
void parametrization_layer_el_aniso_loc(
    size_t nx, 
    size_t nz,
    const float *Gridx, 
    const float *Gridz,
    int grid_type, 
    int media_type,
    inter_t *interfaces,
    float *c11,
    float *c13,
    float *c15,
    float *c33,
    float *c35,
    float *c55,
    float *rho);

//- 4.1 Assign the parameter by volume arithmetic and harmonic averaging method 
//- elasic aniso
void parametrization_layer_el_aniso_har(
    size_t nx, 
    size_t nz,
    const float *Gridx, 
    const float *Gridz,
    int grid_type, 
    int media_type,
    inter_t *interfaces,
    float *c11,
    float *c13,
    float *c15,
    float *c33,
    float *c35,
    float *c55,
    float *rho);

//- 4.2 Assign the parameter by volume arithmetic averaging method 
//- elasic tti
void parametrization_layer_el_aniso_ari(
    size_t nx, 
    size_t nz,
    const float *Gridx, 
    const float *Gridz,
    int grid_type, 
    int media_type,
    inter_t *interfaces,
    float *c11,
    float *c13,
    float *c15,
    float *c33,
    float *c35,
    float *c55,
    float *rho);


void parametrization_layer_el_aniso_tti(
    size_t nx, 
    size_t nz,
    const float *Gridx, 
    const float *Gridz,
    int grid_type, 
    int media_type,
    inter_t *interfaces,
    float *c11,
    float *c13,
    float *c15,
    float *c33,
    float *c35,
    float *c55,
    float *rho);

#endif /* __MEDID_LAYER2MODEL__ */
