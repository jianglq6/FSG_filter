/***************************************************************************
 *
 * This function is used for medium parameterization.
 *
 * Authors: Luqian Jiang <jianglq@mail.ustc.edu.cn>
 *          Wei Zhang <zhangwei@sustech.edu.cn>
 *
 * Copyright (c) 2021 zwlab
 *
 * ChangeLog: 
 *    09/2021: Created by Luqian Jiang 
 *
 ***************************************************************************/
#include <iostream>
#include <vector>
#include <string.h>
#include <cmath>
#include <numeric>
#include <set>
#include "media_geometry2d.hpp"
#include "media_layer2model.hpp"
#include "media_read_file.hpp"
#include "media_utility.hpp"

// for reporting error, media_type int -> string
std::map<int, std::string> md_map = create_md2str_map();
extern int edgeTable[16];

/*============================= for C call =================================*/

//--- 4. elastic anisotropic/TTI
//  the TTI medium or the other medium but used TTI equivalent medium parametrization
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
        const char *equivalent_medium_method) 
{

    inter_t *interfaces;
    int md_type = -1;
    size_t siz_slice = nx*nz;

    /* Read interface file */
    read_interface_file(in_lay_file, &interfaces, md_type);   

    // the function does not support one component and acoustic wave
    // if the 
    if (md_type == ONE_COMPONENT){
        fprintf(stderr, "Error: media_type=one_component is not supported in media_layer2model_el_aniso(),\n"\
                        "       please check the media_type of %s! \n", in_lay_file);        
        fflush(stderr);
        exit(1);
    } else if (md_type == ACOUSTIC_ISOTROPIC){
        fprintf(stderr, "Error: media_type=acoustic_isotropic is not supported in media_layer2model_el_aniso(),\n"\
                        "       please check the media_type of %s! \n", in_lay_file);        
        fflush(stderr);
        exit(1);
    }

    //- isotropic: loc, har, tti
    if (md_type == ELASTIC_ISOTROPIC) {
        fprintf(stderr, "Error: no such equivalent_medium_method: %s!\n", equivalent_medium_method);        
        fprintf(stderr, "Contact the author for this code. \n");        
        fflush(stderr);
        exit(1);        
    // vti: loc, har, ari    
    } else if (md_type == ELASTIC_VTI_PREM || md_type == ELASTIC_VTI_THOMSEN || md_type == ELASTIC_VTI_CIJ) {
        fprintf(stderr, "Error: no such equivalent_medium_method: %s! \n", equivalent_medium_method);        
        fprintf(stderr, "Contact the author for this code. \n");        
        fflush(stderr);
        exit(1);        
    } else { // for anisotropic
        if (strcmp(equivalent_medium_method, "loc") == 0) {
            parametrization_layer_el_aniso_loc(nx, nz, x2d, z2d, grid_type, 
                md_type, interfaces, c11, c13, c15, c33, c35, c55, rho); 
        } else if (strcmp(equivalent_medium_method, "har") == 0) {
            parametrization_layer_el_aniso_har(nx, nz, x2d, z2d, grid_type, 
                md_type, interfaces, c11, c13, c15, c33, c35, c55, rho); 
        } else if (strcmp(equivalent_medium_method, "ari") == 0) {
            parametrization_layer_el_aniso_ari(nx, nz, x2d, z2d, grid_type, 
                md_type, interfaces, c11, c13, c15, c33, c35, c55, rho); 
        } else if(strcmp(equivalent_medium_method,"tti") == 0) {
            parametrization_layer_el_aniso_tti(nx, nz, x2d, z2d, grid_type, 
                md_type, interfaces, c11, c13, c15, c33, c35, c55, rho); 
        } else { //default
            fprintf(stderr, "Error: no such equivalent_medium_method: %s! \n", equivalent_medium_method);        
            fflush(stderr);
            exit(1);        
        }
    }

    delete []interfaces;
    fprintf(stdout, " - Done\n"); 

    return 0; 
}

/*=================================================================================*/

int AssignLayerMediaPara2Point(
    size_t ix, size_t iz,         /* To print Error messages */ 
    Point2 A,  
    inter_t *interfaces,
    int media_type,                /* the type can be found in media_utility.hpp */ 
    std::vector<float> &var)
{
    size_t  NI = interfaces[0].NI;
    int num_out_range = 0;
    int isPointOnInter = 0; 

    std::vector<float> elevation;
    std::vector<int> internum4elev;  // interface number for the elevation, for find the number of interfaces to used
    for (int ni = 0; ni < NI; ni++) {
        int npoint = interfaces[ni].NX;
        float MINX = interfaces[ni].xloc[0];
        float MAXX = interfaces[ni].xloc[npoint-1];  
        if (A.x < MINX || A.x > MAXX) {
            num_out_range++;
        } else {
            float elev = LinearInterp(npoint, interfaces[ni].xloc, interfaces[ni].elevation, A.x);
            elevation.emplace_back(elev); 
            internum4elev.emplace_back(ni);
        }
    }

    // out every interface range, exit!
    if (num_out_range == NI) {
        fprintf(stderr, "Error: Grid(%d, %d)=(%f,%f) is out of every interface,\n"\
                        "       please check the in_lay_file and give a larger interface range in x!\n", 
                        (int)ix, (int)iz, A.x, A.z);
        fflush(stderr);
        exit(1);
    }

    /* Find which material_index to use */
    int mi = findLastGreaterEqualIndex(A.z, elevation);

    if (mi > -1 && isEqual(A.z, elevation[mi])) 
        isPointOnInter = 1;

    CalPointValue_layer(media_type, interfaces, A, elevation, internum4elev, mi, var);

    return isPointOnInter;
}

//- Calculate the value of the point for different media type (to avoid multiple geometric calculations) 
//   for layer2model
void CalPointValue_layer(int media_type, 
                   inter_t *interfaces,
                   Point2 &A,
                   std::vector<float> &elevation,  /*the elevation of point A at the projection position of the interface mesh. */
                   std::vector<int> &internum4elev,
                   int mi,
                   std::vector<float> &var)
{
    float dz = elevation[mi] - A.z;

    switch(media_type)
    {
    case ONE_COMPONENT: /* 0. var */
        /* If grid_z > elevation of top_interface, it given by the medium of top non-zero thickness layer */
        if (mi == -1) {
            mi = findLastGreaterEqualIndex(elevation[0], elevation);
            int k = internum4elev[mi]; 
            var[0] = LinearInterp(interfaces[k].NX, interfaces[k].xloc, interfaces[k].var, A.x); 
        } else {
            int k = internum4elev[mi]; 
            int n = interfaces[k].NX;
            float var0      = LinearInterp(n, interfaces[k].xloc, interfaces[k].var     , A.x);
            float var_grad  = LinearInterp(n, interfaces[k].xloc, interfaces[k].var_grad, A.x);
            float var_pow   = LinearInterp(n, interfaces[k].xloc, interfaces[k].var_pow , A.x);
            var[0]  = var0  + pow(dz, var_pow)* var_grad;
            if (isEqual(dz, 0.0)) {
                int fi = findFirstGreaterEqualIndex(A.z, elevation)-1;
                if (fi >= 0) {
                    k = internum4elev[fi];
                    n = interfaces[k].NX;
                    var0      = LinearInterp(n, interfaces[k].xloc, interfaces[k].var     , A.x);
                    var_grad  = LinearInterp(n, interfaces[k].xloc, interfaces[k].var_grad, A.x);
                    var_pow   = LinearInterp(n, interfaces[k].xloc, interfaces[k].var_pow , A.x);
                    dz = elevation[fi] - A.z;
                    float var_top = var0  + pow(dz, var_pow)* var_grad;
                    var[0] = (var[0] + var_top)/2.0;
                }
            }
        }
    break;

    case ELASTIC_ISOTROPIC: /*1. rho, vp, vs*/
        if (mi == -1) {
            mi = findLastGreaterEqualIndex(elevation[0], elevation);
            int k = internum4elev[mi]; 
            int n = interfaces[k].NX;
            var[0] = LinearInterp(n, interfaces[k].xloc, interfaces[k].rho, A.x);
            var[1] = LinearInterp(n, interfaces[k].xloc, interfaces[k].vp , A.x);
            var[2] = LinearInterp(n, interfaces[k].xloc, interfaces[k].vs , A.x);
        } else {
            int k = internum4elev[mi]; 
            int n = interfaces[k].NX;
            float vp       = LinearInterp(n, interfaces[k].xloc, interfaces[k].vp      , A.x);
            float vp_grad  = LinearInterp(n, interfaces[k].xloc, interfaces[k].vp_grad , A.x);
            float vp_pow   = LinearInterp(n, interfaces[k].xloc, interfaces[k].vp_pow  , A.x);
            float vs       = LinearInterp(n, interfaces[k].xloc, interfaces[k].vs      , A.x);
            float vs_grad  = LinearInterp(n, interfaces[k].xloc, interfaces[k].vs_grad , A.x);
            float vs_pow   = LinearInterp(n, interfaces[k].xloc, interfaces[k].vs_pow  , A.x);
            float rho      = LinearInterp(n, interfaces[k].xloc, interfaces[k].rho     , A.x);
            float rho_grad = LinearInterp(n, interfaces[k].xloc, interfaces[k].rho_grad, A.x);
            float rho_pow  = LinearInterp(n, interfaces[k].xloc, interfaces[k].rho_pow , A.x);
            var[1] = vp  + pow(dz, vp_pow)* vp_grad;
            var[2] = vs  + pow(dz, vs_pow)* vs_grad;
            var[0] = rho + pow(dz,rho_pow)*rho_grad;
            // if the point on the interface, average for loc
            if (isEqual(dz, 0.0)) {
                int fi = findFirstGreaterEqualIndex(A.z, elevation)-1;
                if (fi >= 0) {
                    k = internum4elev[fi];
                    n = interfaces[k].NX;
                    vp       = LinearInterp(n, interfaces[k].xloc, interfaces[k].vp      , A.x);
                    vp_grad  = LinearInterp(n, interfaces[k].xloc, interfaces[k].vp_grad , A.x);
                    vp_pow   = LinearInterp(n, interfaces[k].xloc, interfaces[k].vp_pow  , A.x);
                    vs       = LinearInterp(n, interfaces[k].xloc, interfaces[k].vs      , A.x);
                    vs_grad  = LinearInterp(n, interfaces[k].xloc, interfaces[k].vs_grad , A.x);
                    vs_pow   = LinearInterp(n, interfaces[k].xloc, interfaces[k].vs_pow  , A.x);
                    rho      = LinearInterp(n, interfaces[k].xloc, interfaces[k].rho     , A.x);
                    rho_grad = LinearInterp(n, interfaces[k].xloc, interfaces[k].rho_grad, A.x);
                    rho_pow  = LinearInterp(n, interfaces[k].xloc, interfaces[k].rho_pow , A.x);
                    dz = elevation[fi] - A.z;
                    float vp_top  = vp  + pow(dz,  vp_pow)* vp_grad;
                    float vs_top  = vs  + pow(dz,  vs_pow)* vs_grad;
                    float rho_top = rho + pow(dz, rho_pow)* rho_grad;
                    var[0] = (var[0] + rho_top)/2.0;
                    var[1] = (var[1] +  vp_top)/2.0;
                    var[2] = (var[2] +  vs_top)/2.0;
                }
            }
        }
    break;

    case ELASTIC_VTI_PREM: /*2. rho, vph, vpv, vsv, eta */
        if (mi == -1) {
            mi = findLastGreaterEqualIndex(elevation[0], elevation);
            int k = internum4elev[mi]; 
            int n = interfaces[k].NX;
            var[0] = LinearInterp(n, interfaces[k].xloc, interfaces[k].rho, A.x);
            var[1] = LinearInterp(n, interfaces[k].xloc, interfaces[k].vph, A.x);
            var[2] = LinearInterp(n, interfaces[k].xloc, interfaces[k].vpv, A.x);
            var[3] = LinearInterp(n, interfaces[k].xloc, interfaces[k].vsv, A.x);
            var[4] = LinearInterp(n, interfaces[k].xloc, interfaces[k].eta, A.x);  
        } else {
            int k = internum4elev[mi]; 
            int n = interfaces[k].NX;
            float vph = LinearInterp(n, interfaces[k].xloc, interfaces[k].vph, A.x);
            float vpv = LinearInterp(n, interfaces[k].xloc, interfaces[k].vpv, A.x);
            float vsv = LinearInterp(n, interfaces[k].xloc, interfaces[k].vsv, A.x);
            float rho = LinearInterp(n, interfaces[k].xloc, interfaces[k].rho, A.x);
            float eta = LinearInterp(n, interfaces[k].xloc, interfaces[k].eta, A.x);
            float vph_grad = LinearInterp(n, interfaces[k].xloc, interfaces[k].vph_grad, A.x);
            float vpv_grad = LinearInterp(n, interfaces[k].xloc, interfaces[k].vpv_grad, A.x);
            float vsv_grad = LinearInterp(n, interfaces[k].xloc, interfaces[k].vsv_grad, A.x);
            float rho_grad = LinearInterp(n, interfaces[k].xloc, interfaces[k].rho_grad, A.x);
            float eta_grad = LinearInterp(n, interfaces[k].xloc, interfaces[k].eta_grad, A.x);
            float vph_pow  = LinearInterp(n, interfaces[k].xloc, interfaces[k].vph_pow , A.x);
            float vpv_pow  = LinearInterp(n, interfaces[k].xloc, interfaces[k].vpv_pow , A.x);
            float vsv_pow  = LinearInterp(n, interfaces[k].xloc, interfaces[k].vsv_pow , A.x);
            float rho_pow  = LinearInterp(n, interfaces[k].xloc, interfaces[k].rho_pow , A.x);
            float eta_pow  = LinearInterp(n, interfaces[k].xloc, interfaces[k].eta_pow , A.x);
            var[0] = rho + pow(dz, rho_pow)* rho_grad;
            var[1] = vph + pow(dz, vph_pow)* vph_grad;
            var[2] = vpv + pow(dz, vpv_pow)* vpv_grad;
            var[3] = vsv + pow(dz, vsv_pow)* vsv_grad;
            var[4] = eta + pow(dz, eta_pow)* eta_grad;
            if (isEqual(dz, 0.0)) {
                int fi = findFirstGreaterEqualIndex(A.z, elevation)-1;
                if (fi >= 0) {
                    k = internum4elev[fi];
                    n = interfaces[k].NX;
                    vph = LinearInterp(n, interfaces[k].xloc, interfaces[k].vph, A.x);
                    vpv = LinearInterp(n, interfaces[k].xloc, interfaces[k].vpv, A.x);
                    vsv = LinearInterp(n, interfaces[k].xloc, interfaces[k].vsv, A.x);
                    rho = LinearInterp(n, interfaces[k].xloc, interfaces[k].rho, A.x);
                    eta = LinearInterp(n, interfaces[k].xloc, interfaces[k].eta, A.x);
                    vph_grad = LinearInterp(n, interfaces[k].xloc, interfaces[k].vph_grad, A.x);
                    vpv_grad = LinearInterp(n, interfaces[k].xloc, interfaces[k].vpv_grad, A.x);
                    vsv_grad = LinearInterp(n, interfaces[k].xloc, interfaces[k].vsv_grad, A.x);
                    rho_grad = LinearInterp(n, interfaces[k].xloc, interfaces[k].rho_grad, A.x);
                    eta_grad = LinearInterp(n, interfaces[k].xloc, interfaces[k].eta_grad, A.x);
                    vph_pow  = LinearInterp(n, interfaces[k].xloc, interfaces[k].vph_pow , A.x);
                    vpv_pow  = LinearInterp(n, interfaces[k].xloc, interfaces[k].vpv_pow , A.x);
                    vsv_pow  = LinearInterp(n, interfaces[k].xloc, interfaces[k].vsv_pow , A.x);
                    rho_pow  = LinearInterp(n, interfaces[k].xloc, interfaces[k].rho_pow , A.x);
                    eta_pow  = LinearInterp(n, interfaces[k].xloc, interfaces[k].eta_pow , A.x);
                    dz = elevation[fi] - A.z;
                    float rho_top = rho + pow(dz, rho_pow)* rho_grad;
                    float vph_top = vph + pow(dz, vph_pow)* vph_grad;
                    float vpv_top = vpv + pow(dz, vpv_pow)* vpv_grad;
                    float vsv_top = vsv + pow(dz, vsv_pow)* vsv_grad;
                    float eta_top = eta + pow(dz, eta_pow)* eta_grad;
                    var[0] = (var[0] + rho_top)/2.0;
                    var[1] = (var[1] + vph_top)/2.0;
                    var[2] = (var[2] + vpv_top)/2.0;
                    var[3] = (var[3] + vsv_top)/2.0;
                    var[4] = (var[4] + eta_top)/2.0;
                }
            }
        }   
    break;

    case ELASTIC_VTI_THOMSEN: /*3. rho, vp0, vs0, epsilon, delta */
        if (mi == -1) {
            mi = findLastGreaterEqualIndex(elevation[0], elevation);
            int k = internum4elev[mi]; 
            int n = interfaces[k].NX;
            var[0] = LinearInterp(n, interfaces[k].xloc, interfaces[k].rho    , A.x);
            var[1] = LinearInterp(n, interfaces[k].xloc, interfaces[k].vp0    , A.x);
            var[2] = LinearInterp(n, interfaces[k].xloc, interfaces[k].vs0    , A.x);
            var[3] = LinearInterp(n, interfaces[k].xloc, interfaces[k].epsilon, A.x);
            var[4] = LinearInterp(n, interfaces[k].xloc, interfaces[k].delta  , A.x);
        } else {
            int k = internum4elev[mi]; 
            int n = interfaces[k].NX;
            float rho     = LinearInterp(n, interfaces[k].xloc, interfaces[k].rho    , A.x);
            float vp0     = LinearInterp(n, interfaces[k].xloc, interfaces[k].vp0    , A.x);
            float vs0     = LinearInterp(n, interfaces[k].xloc, interfaces[k].vs0    , A.x);
            float epsil   = LinearInterp(n, interfaces[k].xloc, interfaces[k].epsilon, A.x);  // epsilon
            float delta   = LinearInterp(n, interfaces[k].xloc, interfaces[k].delta  , A.x);
            float rho_grad     = LinearInterp(n, interfaces[k].xloc, interfaces[k].rho_grad    , A.x);
            float vp0_grad     = LinearInterp(n, interfaces[k].xloc, interfaces[k].vp0_grad    , A.x);
            float vs0_grad     = LinearInterp(n, interfaces[k].xloc, interfaces[k].vs0_grad    , A.x);
            float epsilon_grad = LinearInterp(n, interfaces[k].xloc, interfaces[k].epsilon_grad, A.x);
            float delta_grad   = LinearInterp(n, interfaces[k].xloc, interfaces[k].delta_grad  , A.x);
            float rho_pow      = LinearInterp(n, interfaces[k].xloc, interfaces[k].rho_pow     , A.x);
            float vp0_pow      = LinearInterp(n, interfaces[k].xloc, interfaces[k].vp0_pow     , A.x);
            float vs0_pow      = LinearInterp(n, interfaces[k].xloc, interfaces[k].vs0_pow     , A.x);
            float epsilon_pow  = LinearInterp(n, interfaces[k].xloc, interfaces[k].epsilon_pow , A.x);
            float delta_pow    = LinearInterp(n, interfaces[k].xloc, interfaces[k].delta_pow   , A.x);
            var[0] = rho     + pow(dz, rho_pow)    * rho_grad;
            var[1] = vp0     + pow(dz, vp0_pow)    * vp0_grad;
            var[2] = vs0     + pow(dz, vs0_pow)    * vs0_grad;
            var[3] = epsil   + pow(dz, epsilon_pow)* epsilon_grad;
            var[4] = delta   + pow(dz, delta_pow)  * delta_grad;
            if (isEqual(dz, 0.0)) {
                int fi = findFirstGreaterEqualIndex(A.z, elevation)-1;
                if (fi >= 0) {
                    k = internum4elev[fi];
                    n = interfaces[k].NX;
                    rho     = LinearInterp(n, interfaces[k].xloc, interfaces[k].rho    , A.x);
                    vp0     = LinearInterp(n, interfaces[k].xloc, interfaces[k].vp0    , A.x);
                    vs0     = LinearInterp(n, interfaces[k].xloc, interfaces[k].vs0    , A.x);
                    epsil   = LinearInterp(n, interfaces[k].xloc, interfaces[k].epsilon, A.x);  // epsilon
                    delta   = LinearInterp(n, interfaces[k].xloc, interfaces[k].delta  , A.x);
                    rho_grad     = LinearInterp(n, interfaces[k].xloc, interfaces[k].rho_grad    , A.x);
                    vp0_grad     = LinearInterp(n, interfaces[k].xloc, interfaces[k].vp0_grad    , A.x);
                    vs0_grad     = LinearInterp(n, interfaces[k].xloc, interfaces[k].vs0_grad    , A.x);
                    epsilon_grad = LinearInterp(n, interfaces[k].xloc, interfaces[k].epsilon_grad, A.x);
                    delta_grad   = LinearInterp(n, interfaces[k].xloc, interfaces[k].delta_grad  , A.x);
                    rho_pow      = LinearInterp(n, interfaces[k].xloc, interfaces[k].rho_pow     , A.x);
                    vp0_pow      = LinearInterp(n, interfaces[k].xloc, interfaces[k].vp0_pow     , A.x);
                    vs0_pow      = LinearInterp(n, interfaces[k].xloc, interfaces[k].vs0_pow     , A.x);
                    epsilon_pow  = LinearInterp(n, interfaces[k].xloc, interfaces[k].epsilon_pow , A.x);
                    delta_pow    = LinearInterp(n, interfaces[k].xloc, interfaces[k].delta_pow   , A.x);
                    dz = elevation[fi] - A.z;
                    float var0_top = rho     + pow(dz, rho_pow)    * rho_grad;
                    float var1_top = vp0     + pow(dz, vp0_pow)    * vp0_grad;
                    float var2_top = vs0     + pow(dz, vs0_pow)    * vs0_grad;
                    float var3_top = epsil   + pow(dz, epsilon_pow)* epsilon_grad;
                    float var4_top = delta   + pow(dz, delta_pow)  * delta_grad;
                    var[0] = (var[0] + var0_top)/2.0;
                    var[1] = (var[1] + var1_top)/2.0;
                    var[2] = (var[2] + var2_top)/2.0;
                    var[3] = (var[3] + var3_top)/2.0;
                    var[4] = (var[4] + var4_top)/2.0;
                }
            }

        }   
    break;

    case ELASTIC_VTI_CIJ: /*4. rho c11 c33 c55 c66 c13 */
         if (mi == -1) {
            mi = findLastGreaterEqualIndex(elevation[0], elevation);
            int k = internum4elev[mi]; 
            int n = interfaces[k].NX;
            var[0] = LinearInterp(n, interfaces[k].xloc, interfaces[k].rho, A.x);
            var[1] = LinearInterp(n, interfaces[k].xloc, interfaces[k].c11, A.x);
            var[2] = LinearInterp(n, interfaces[k].xloc, interfaces[k].c33, A.x);
            var[3] = LinearInterp(n, interfaces[k].xloc, interfaces[k].c55, A.x);
            var[4] = LinearInterp(n, interfaces[k].xloc, interfaces[k].c13, A.x);   
        } else {
            int k = internum4elev[mi]; 
            int n = interfaces[k].NX;
            float rho = LinearInterp(n, interfaces[k].xloc, interfaces[k].rho, A.x);
            float c11 = LinearInterp(n, interfaces[k].xloc, interfaces[k].c11, A.x);
            float c33 = LinearInterp(n, interfaces[k].xloc, interfaces[k].c33, A.x);
            float c55 = LinearInterp(n, interfaces[k].xloc, interfaces[k].c55, A.x);
            float c13 = LinearInterp(n, interfaces[k].xloc, interfaces[k].c13, A.x);
            float c11_grad = LinearInterp(n, interfaces[k].xloc, interfaces[k].c11_grad, A.x);
            float c33_grad = LinearInterp(n, interfaces[k].xloc, interfaces[k].c33_grad, A.x);
            float c55_grad = LinearInterp(n, interfaces[k].xloc, interfaces[k].c55_grad, A.x);
            float c13_grad = LinearInterp(n, interfaces[k].xloc, interfaces[k].c13_grad, A.x);
            float rho_grad = LinearInterp(n, interfaces[k].xloc, interfaces[k].rho_grad, A.x);
            float c11_pow  = LinearInterp(n, interfaces[k].xloc, interfaces[k].c11_pow , A.x);
            float c33_pow  = LinearInterp(n, interfaces[k].xloc, interfaces[k].c33_pow , A.x);
            float c55_pow  = LinearInterp(n, interfaces[k].xloc, interfaces[k].c55_pow , A.x);
            float c13_pow  = LinearInterp(n, interfaces[k].xloc, interfaces[k].c13_pow , A.x);
            float rho_pow  = LinearInterp(n, interfaces[k].xloc, interfaces[k].rho_pow , A.x);
            var[0] = rho + pow(dz, rho_pow)* rho_grad;
            var[1] = c11 + pow(dz, c11_pow)* c11_grad;
            var[2] = c33 + pow(dz, c33_pow)* c33_grad;
            var[3] = c55 + pow(dz, c55_pow)* c55_grad;
            var[4] = c13 + pow(dz, c13_pow)* c13_grad;
            if (isEqual(dz, 0.0)) {
                int fi = findFirstGreaterEqualIndex(A.z, elevation)-1;
                if (fi >= 0) {
                    k = internum4elev[fi];
                    n = interfaces[k].NX;
                    rho = LinearInterp(n, interfaces[k].xloc, interfaces[k].rho, A.x);
                    c11 = LinearInterp(n, interfaces[k].xloc, interfaces[k].c11, A.x);
                    c33 = LinearInterp(n, interfaces[k].xloc, interfaces[k].c33, A.x);
                    c55 = LinearInterp(n, interfaces[k].xloc, interfaces[k].c55, A.x);
                    c13 = LinearInterp(n, interfaces[k].xloc, interfaces[k].c13, A.x);
                    c11_grad = LinearInterp(n, interfaces[k].xloc, interfaces[k].c11_grad, A.x);
                    c33_grad = LinearInterp(n, interfaces[k].xloc, interfaces[k].c33_grad, A.x);
                    c55_grad = LinearInterp(n, interfaces[k].xloc, interfaces[k].c55_grad, A.x);
                    c13_grad = LinearInterp(n, interfaces[k].xloc, interfaces[k].c13_grad, A.x);
                    rho_grad = LinearInterp(n, interfaces[k].xloc, interfaces[k].rho_grad, A.x);
                    c11_pow  = LinearInterp(n, interfaces[k].xloc, interfaces[k].c11_pow , A.x);
                    c33_pow  = LinearInterp(n, interfaces[k].xloc, interfaces[k].c33_pow , A.x);
                    c55_pow  = LinearInterp(n, interfaces[k].xloc, interfaces[k].c55_pow , A.x);
                    c13_pow  = LinearInterp(n, interfaces[k].xloc, interfaces[k].c13_pow , A.x);
                    rho_pow  = LinearInterp(n, interfaces[k].xloc, interfaces[k].rho_pow , A.x);
                    dz = elevation[fi] - A.z;
                    float rho_top = rho + pow(dz, rho_pow)* rho_grad;
                    float c11_top = c11 + pow(dz, c11_pow)* c11_grad;
                    float c33_top = c33 + pow(dz, c33_pow)* c33_grad;
                    float c55_top = c55 + pow(dz, c55_pow)* c55_grad;
                    float c13_top = c13 + pow(dz, c13_pow)* c13_grad;
                    var[0] = (var[0] + rho_top)/2.0;
                    var[1] = (var[1] + c11_top)/2.0;
                    var[2] = (var[2] + c33_top)/2.0;
                    var[3] = (var[3] + c55_top)/2.0;
                    var[4] = (var[4] + c13_top)/2.0;
                }
            }
        }   
    break;

    case ELASTIC_TTI_THOMSEN: /*5. rho, vp0, vs0, epsilon, delta, dip */
        if (mi == -1) {
            mi = findLastGreaterEqualIndex(elevation[0], elevation);
            int k = internum4elev[mi]; 
            int n = interfaces[k].NX;
            var[0] = LinearInterp(n, interfaces[k].xloc, interfaces[k].rho    , A.x);
            var[1] = LinearInterp(n, interfaces[k].xloc, interfaces[k].vp0    , A.x);
            var[2] = LinearInterp(n, interfaces[k].xloc, interfaces[k].vs0    , A.x);
            var[3] = LinearInterp(n, interfaces[k].xloc, interfaces[k].epsilon, A.x);
            var[4] = LinearInterp(n, interfaces[k].xloc, interfaces[k].delta  , A.x);
            var[5] = LinearInterp(n, interfaces[k].xloc, interfaces[k].dip    , A.x);
        } else {
            int k = internum4elev[mi]; 
            int n = interfaces[k].NX;
            float vp0     = LinearInterp(n, interfaces[k].xloc, interfaces[k].vp0    , A.x);
            float vs0     = LinearInterp(n, interfaces[k].xloc, interfaces[k].vs0    , A.x);
            float epsil   = LinearInterp(n, interfaces[k].xloc, interfaces[k].epsilon, A.x);
            float delta   = LinearInterp(n, interfaces[k].xloc, interfaces[k].delta  , A.x);
            float rho     = LinearInterp(n, interfaces[k].xloc, interfaces[k].rho    , A.x);
            float dip     = LinearInterp(n, interfaces[k].xloc, interfaces[k].dip    , A.x);
            float vp0_grad     = LinearInterp(n, interfaces[k].xloc, interfaces[k].vp0_grad    , A.x);
            float vs0_grad     = LinearInterp(n, interfaces[k].xloc, interfaces[k].vs0_grad    , A.x);
            float epsilon_grad = LinearInterp(n, interfaces[k].xloc, interfaces[k].epsilon_grad, A.x);
            float delta_grad   = LinearInterp(n, interfaces[k].xloc, interfaces[k].delta_grad  , A.x);
            float rho_grad     = LinearInterp(n, interfaces[k].xloc, interfaces[k].rho_grad    , A.x);
            float dip_grad     = LinearInterp(n, interfaces[k].xloc, interfaces[k].dip_grad    , A.x);
            float vp0_pow      = LinearInterp(n, interfaces[k].xloc, interfaces[k].vp0_pow     , A.x);
            float vs0_pow      = LinearInterp(n, interfaces[k].xloc, interfaces[k].vs0_pow     , A.x);
            float epsilon_pow  = LinearInterp(n, interfaces[k].xloc, interfaces[k].epsilon_pow , A.x);
            float delta_pow    = LinearInterp(n, interfaces[k].xloc, interfaces[k].delta_pow   , A.x);
            float rho_pow      = LinearInterp(n, interfaces[k].xloc, interfaces[k].rho_pow     , A.x);
            float dip_pow      = LinearInterp(n, interfaces[k].xloc, interfaces[k].dip_pow     , A.x);

            var[0] = rho     + pow(dz, rho_pow)    * rho_grad;
            var[1] = vp0     + pow(dz, vp0_pow)    * vp0_grad;
            var[2] = vs0     + pow(dz, vs0_pow)    * vs0_grad;
            var[3] = epsil   + pow(dz, epsilon_pow)* epsilon_grad;
            var[4] = delta   + pow(dz, delta_pow)  * delta_grad;
            var[5] = dip     + pow(dz, dip_pow)    * dip_grad;

            if (isEqual(dz, 0.0)) {
                int fi = findFirstGreaterEqualIndex(A.z, elevation)-1;
                if (fi >= 0) {
                    k = internum4elev[fi];
                    n = interfaces[k].NX;
                    vp0     = LinearInterp(n, interfaces[k].xloc, interfaces[k].vp0    , A.x);
                    vs0     = LinearInterp(n, interfaces[k].xloc, interfaces[k].vs0    , A.x);
                    epsil   = LinearInterp(n, interfaces[k].xloc, interfaces[k].epsilon, A.x);
                    delta   = LinearInterp(n, interfaces[k].xloc, interfaces[k].delta  , A.x);
                    rho     = LinearInterp(n, interfaces[k].xloc, interfaces[k].rho    , A.x);
                    dip     = LinearInterp(n, interfaces[k].xloc, interfaces[k].dip    , A.x);
                    vp0_grad     = LinearInterp(n, interfaces[k].xloc, interfaces[k].vp0_grad    , A.x);
                    vs0_grad     = LinearInterp(n, interfaces[k].xloc, interfaces[k].vs0_grad    , A.x);
                    epsilon_grad = LinearInterp(n, interfaces[k].xloc, interfaces[k].epsilon_grad, A.x);
                    delta_grad   = LinearInterp(n, interfaces[k].xloc, interfaces[k].delta_grad  , A.x);
                    rho_grad     = LinearInterp(n, interfaces[k].xloc, interfaces[k].rho_grad    , A.x);
                    dip_grad     = LinearInterp(n, interfaces[k].xloc, interfaces[k].dip_grad    , A.x);
                    vp0_pow      = LinearInterp(n, interfaces[k].xloc, interfaces[k].vp0_pow     , A.x);
                    vs0_pow      = LinearInterp(n, interfaces[k].xloc, interfaces[k].vs0_pow     , A.x);
                    epsilon_pow  = LinearInterp(n, interfaces[k].xloc, interfaces[k].epsilon_pow , A.x);
                    delta_pow    = LinearInterp(n, interfaces[k].xloc, interfaces[k].delta_pow   , A.x);
                    rho_pow      = LinearInterp(n, interfaces[k].xloc, interfaces[k].rho_pow     , A.x);
                    dip_pow      = LinearInterp(n, interfaces[k].xloc, interfaces[k].dip_pow     , A.x);

                    dz = elevation[fi] - A.z;
                    float var0_top = rho   + pow(dz, rho_pow)    * rho_grad;
                    float var1_top = vp0   + pow(dz, vp0_pow)    * vp0_grad;
                    float var2_top = vs0   + pow(dz, vs0_pow)    * vs0_grad;
                    float var3_top = epsil + pow(dz, epsilon_pow)* epsilon_grad;
                    float var4_top = delta + pow(dz, delta_pow)  * delta_grad;
                    float var5_top = dip   + pow(dz, dip_pow)    * dip_grad;

                    var[0] = (var[0] + var0_top)/2.0;
                    var[1] = (var[1] + var1_top)/2.0;
                    var[2] = (var[2] + var2_top)/2.0;
                    var[3] = (var[3] + var3_top)/2.0;
                    var[4] = (var[4] + var4_top)/2.0;
                    var[5] = (var[5] + var5_top)/2.0;

                }
            }
        }   
    break;

    case ELASTIC_ANISO_CIJ: /* 7. rho c11 c13 c15 c33 c35 c55 */
         if (mi == -1) {
            mi = findLastGreaterEqualIndex(elevation[0], elevation);
            int k = internum4elev[mi]; 
            int n = interfaces[k].NX;
            var[0] = LinearInterp(n, interfaces[k].xloc, interfaces[k].rho, A.x);
            var[1] = LinearInterp(n, interfaces[k].xloc, interfaces[k].c11, A.x);
            var[2] = LinearInterp(n, interfaces[k].xloc, interfaces[k].c13, A.x);
            var[3] = LinearInterp(n, interfaces[k].xloc, interfaces[k].c15, A.x);
            var[4] = LinearInterp(n, interfaces[k].xloc, interfaces[k].c33, A.x);
            var[5] = LinearInterp(n, interfaces[k].xloc, interfaces[k].c35, A.x);
            var[6] = LinearInterp(n, interfaces[k].xloc, interfaces[k].c55, A.x);
        } else {
            int k = internum4elev[mi]; 
            int n = interfaces[k].NX;
            float c11 = LinearInterp(n, interfaces[k].xloc, interfaces[k].c11, A.x);
            float c13 = LinearInterp(n, interfaces[k].xloc, interfaces[k].c13, A.x);
            float c15 = LinearInterp(n, interfaces[k].xloc, interfaces[k].c15, A.x);
            float c33 = LinearInterp(n, interfaces[k].xloc, interfaces[k].c33, A.x);
            float c35 = LinearInterp(n, interfaces[k].xloc, interfaces[k].c35, A.x);
            float c55 = LinearInterp(n, interfaces[k].xloc, interfaces[k].c55, A.x);
            float rho = LinearInterp(n, interfaces[k].xloc, interfaces[k].rho, A.x);
            float c11_pow = LinearInterp( n, interfaces[k].xloc, interfaces[k].c11_pow,  A.x);
            float c13_pow = LinearInterp( n, interfaces[k].xloc, interfaces[k].c13_pow,  A.x);
            float c15_pow = LinearInterp( n, interfaces[k].xloc, interfaces[k].c15_pow,  A.x);
            float c33_pow = LinearInterp( n, interfaces[k].xloc, interfaces[k].c33_pow,  A.x);
            float c35_pow = LinearInterp( n, interfaces[k].xloc, interfaces[k].c35_pow,  A.x);
            float c55_pow = LinearInterp( n, interfaces[k].xloc, interfaces[k].c55_pow,  A.x);
            float rho_pow = LinearInterp( n, interfaces[k].xloc, interfaces[k].rho_pow,  A.x);
            float c11_grad = LinearInterp(n, interfaces[k].xloc, interfaces[k].c11_grad, A.x);
            float c13_grad = LinearInterp(n, interfaces[k].xloc, interfaces[k].c13_grad, A.x);
            float c15_grad = LinearInterp(n, interfaces[k].xloc, interfaces[k].c15_grad, A.x);
            float c33_grad = LinearInterp(n, interfaces[k].xloc, interfaces[k].c33_grad, A.x);
            float c35_grad = LinearInterp(n, interfaces[k].xloc, interfaces[k].c35_grad, A.x);
            float c55_grad = LinearInterp(n, interfaces[k].xloc, interfaces[k].c55_grad, A.x);
            float rho_grad = LinearInterp(n, interfaces[k].xloc, interfaces[k].rho_grad, A.x);
            var[0] = rho + pow(dz, rho_pow) * rho_grad;
            var[1] = c11 + pow(dz, c11_pow) * c11_grad; 
            var[2] = c13 + pow(dz, c13_pow) * c13_grad;
            var[3] = c15 + pow(dz, c15_pow) * c15_grad;
            var[4] = c33 + pow(dz, c33_pow) * c33_grad;
            var[5] = c35 + pow(dz, c35_pow) * c35_grad;
            var[6] = c55 + pow(dz, c55_pow) * c55_grad;

            if (isEqual(dz, 0.0)) {
                int fi = findFirstGreaterEqualIndex(A.z, elevation)-1;
                if (fi >= 0) {
                    k = internum4elev[fi];
                    n = interfaces[k].NX;
                    c11 = LinearInterp(n, interfaces[k].xloc, interfaces[k].c11, A.x);
                    c13 = LinearInterp(n, interfaces[k].xloc, interfaces[k].c13, A.x);
                    c15 = LinearInterp(n, interfaces[k].xloc, interfaces[k].c15, A.x);
                    c33 = LinearInterp(n, interfaces[k].xloc, interfaces[k].c33, A.x);
                    c35 = LinearInterp(n, interfaces[k].xloc, interfaces[k].c35, A.x);
                    c55 = LinearInterp(n, interfaces[k].xloc, interfaces[k].c55, A.x);
                    rho = LinearInterp(n, interfaces[k].xloc, interfaces[k].rho, A.x);
                    c11_pow = LinearInterp( n, interfaces[k].xloc, interfaces[k].c11_pow,  A.x);
                    c13_pow = LinearInterp( n, interfaces[k].xloc, interfaces[k].c13_pow,  A.x);
                    c15_pow = LinearInterp( n, interfaces[k].xloc, interfaces[k].c15_pow,  A.x);
                    c33_pow = LinearInterp( n, interfaces[k].xloc, interfaces[k].c33_pow,  A.x);
                    c35_pow = LinearInterp( n, interfaces[k].xloc, interfaces[k].c35_pow,  A.x);
                    c55_pow = LinearInterp( n, interfaces[k].xloc, interfaces[k].c55_pow,  A.x);
                    rho_pow = LinearInterp( n, interfaces[k].xloc, interfaces[k].rho_pow,  A.x);
                    c11_grad = LinearInterp(n, interfaces[k].xloc, interfaces[k].c11_grad, A.x);
                    c13_grad = LinearInterp(n, interfaces[k].xloc, interfaces[k].c13_grad, A.x);
                    c15_grad = LinearInterp(n, interfaces[k].xloc, interfaces[k].c15_grad, A.x);
                    c33_grad = LinearInterp(n, interfaces[k].xloc, interfaces[k].c33_grad, A.x);
                    c35_grad = LinearInterp(n, interfaces[k].xloc, interfaces[k].c35_grad, A.x);
                    c55_grad = LinearInterp(n, interfaces[k].xloc, interfaces[k].c55_grad, A.x);
                    rho_grad = LinearInterp(n, interfaces[k].xloc, interfaces[k].rho_grad, A.x);
                    dz = elevation[fi] - A.z;
                    float rho_top = rho + pow(dz, rho_pow) * rho_grad;
                    float c11_top = c11 + pow(dz, c11_pow) * c11_grad;
                    float c33_top = c33 + pow(dz, c33_pow) * c33_grad;
                    float c55_top = c55 + pow(dz, c55_pow) * c55_grad;
                    float c13_top = c13 + pow(dz, c13_pow) * c13_grad;
                    float c35_top = c35 + pow(dz, c35_pow) * c35_grad;
                    float c15_top = c15 + pow(dz, c15_pow) * c15_grad;
                    var[0] = (rho_top + var[0])/2.0;
                    var[1] = (c11_top + var[1])/2.0; 
                    var[2] = (c13_top + var[2])/2.0;
                    var[3] = (c15_top + var[3])/2.0;
                    var[4] = (c33_top + var[4])/2.0;
                    var[5] = (c35_top + var[5])/2.0;
                    var[6] = (c55_top + var[6])/2.0;

                }
            }
        }     
    break;

    case ELASTIC_TTI_BOND: /* 6. rho c11 c33 c55 c13 dip */
        if (mi == -1) {
            mi = findLastGreaterEqualIndex(elevation[0], elevation);
            int k = internum4elev[mi]; 
            int n = interfaces[k].NX;
            var[0] = LinearInterp(n, interfaces[k].xloc, interfaces[k].rho, A.x);
            var[1] = LinearInterp(n, interfaces[k].xloc, interfaces[k].c11, A.x);
            var[2] = LinearInterp(n, interfaces[k].xloc, interfaces[k].c33, A.x);
            var[3] = LinearInterp(n, interfaces[k].xloc, interfaces[k].c55, A.x);
            var[4] = LinearInterp(n, interfaces[k].xloc, interfaces[k].c13, A.x);
            var[5] = LinearInterp(n, interfaces[k].xloc, interfaces[k].dip, A.x);
        } else {
            int k = internum4elev[mi]; 
            int n = interfaces[k].NX;
            float rho = LinearInterp(n, interfaces[k].xloc, interfaces[k].rho, A.x);
            float c11 = LinearInterp(n, interfaces[k].xloc, interfaces[k].c11, A.x);
            float c33 = LinearInterp(n, interfaces[k].xloc, interfaces[k].c33, A.x);
            float c55 = LinearInterp(n, interfaces[k].xloc, interfaces[k].c55, A.x);
            float c13 = LinearInterp(n, interfaces[k].xloc, interfaces[k].c13, A.x);
            float dip = LinearInterp(n, interfaces[k].xloc, interfaces[k].dip, A.x);
            float c11_grad = LinearInterp(n, interfaces[k].xloc, interfaces[k].c11_grad, A.x);
            float c33_grad = LinearInterp(n, interfaces[k].xloc, interfaces[k].c33_grad, A.x);
            float c55_grad = LinearInterp(n, interfaces[k].xloc, interfaces[k].c55_grad, A.x);
            float c13_grad = LinearInterp(n, interfaces[k].xloc, interfaces[k].c13_grad, A.x);
            float rho_grad = LinearInterp(n, interfaces[k].xloc, interfaces[k].rho_grad, A.x);
            float dip_grad = LinearInterp(n, interfaces[k].xloc, interfaces[k].dip_grad, A.x);
            float c11_pow  = LinearInterp(n, interfaces[k].xloc, interfaces[k].c11_pow , A.x);
            float c33_pow  = LinearInterp(n, interfaces[k].xloc, interfaces[k].c33_pow , A.x);
            float c55_pow  = LinearInterp(n, interfaces[k].xloc, interfaces[k].c55_pow , A.x);
            float c13_pow  = LinearInterp(n, interfaces[k].xloc, interfaces[k].c13_pow , A.x);
            float rho_pow  = LinearInterp(n, interfaces[k].xloc, interfaces[k].rho_pow , A.x);
            float dip_pow  = LinearInterp(n, interfaces[k].xloc, interfaces[k].dip_pow , A.x);
            var[0] = rho + pow(dz, rho_pow)* rho_grad;
            var[1] = c11 + pow(dz, c11_pow)* c11_grad;
            var[2] = c33 + pow(dz, c33_pow)* c33_grad;
            var[3] = c55 + pow(dz, c55_pow)* c55_grad;
            var[4] = c13 + pow(dz, c13_pow)* c13_grad;
            var[5] = dip + pow(dz, dip_pow)* dip_grad;

            if (isEqual(dz, 0.0)) {
                int fi = findFirstGreaterEqualIndex(A.z, elevation)-1;
                if (fi >= 0) {
                    k = internum4elev[fi];
                    n = interfaces[k].NX;
                    rho = LinearInterp(n, interfaces[k].xloc, interfaces[k].rho, A.x);
                    c11 = LinearInterp(n, interfaces[k].xloc, interfaces[k].c11, A.x);
                    c33 = LinearInterp(n, interfaces[k].xloc, interfaces[k].c33, A.x);
                    c55 = LinearInterp(n, interfaces[k].xloc, interfaces[k].c55, A.x);
                    c13 = LinearInterp(n, interfaces[k].xloc, interfaces[k].c13, A.x);
                    dip = LinearInterp(n, interfaces[k].xloc, interfaces[k].dip, A.x);
                    dip_grad = LinearInterp(n, interfaces[k].xloc, interfaces[k].dip_grad, A.x);
                    dip_pow  = LinearInterp(n, interfaces[k].xloc, interfaces[k].dip_pow , A.x);
                    c11_grad = LinearInterp(n, interfaces[k].xloc, interfaces[k].c11_grad, A.x);
                    c33_grad = LinearInterp(n, interfaces[k].xloc, interfaces[k].c33_grad, A.x);
                    c55_grad = LinearInterp(n, interfaces[k].xloc, interfaces[k].c55_grad, A.x);
                    c13_grad = LinearInterp(n, interfaces[k].xloc, interfaces[k].c13_grad, A.x);
                    rho_grad = LinearInterp(n, interfaces[k].xloc, interfaces[k].rho_grad, A.x);
                    c11_pow  = LinearInterp(n, interfaces[k].xloc, interfaces[k].c11_pow , A.x);
                    c33_pow  = LinearInterp(n, interfaces[k].xloc, interfaces[k].c33_pow , A.x);
                    c55_pow  = LinearInterp(n, interfaces[k].xloc, interfaces[k].c55_pow , A.x);
                    c13_pow  = LinearInterp(n, interfaces[k].xloc, interfaces[k].c13_pow , A.x);
                    rho_pow  = LinearInterp(n, interfaces[k].xloc, interfaces[k].rho_pow , A.x);
                    dz = elevation[fi] - A.z;
                    float rho_top = rho + pow(dz, rho_pow)* rho_grad;
                    float c11_top = c11 + pow(dz, c11_pow)* c11_grad;
                    float c33_top = c33 + pow(dz, c33_pow)* c33_grad;
                    float c55_top = c55 + pow(dz, c55_pow)* c55_grad;
                    float c13_top = c13 + pow(dz, c13_pow)* c13_grad;
                    float dip_top = dip + pow(dz, dip_pow)* dip_grad;
                    var[0] = (var[0] + rho_top)/2.0;
                    var[1] = (var[1] + c11_top)/2.0;
                    var[2] = (var[2] + c33_top)/2.0;
                    var[3] = (var[3] + c55_top)/2.0;
                    var[4] = (var[4] + c13_top)/2.0;
                    var[5] = (var[5] + dip_top)/2.0;
                }
            }
        }   
    break;

    case ACOUSTIC_ISOTROPIC: /* 7. rho vp */
        if (mi == -1) {
            mi = findLastGreaterEqualIndex(elevation[0], elevation);
            int k = internum4elev[mi]; 
            int n = interfaces[k].NX;
            var[0] = LinearInterp(n, interfaces[k].xloc, interfaces[k].rho, A.x);
            var[1] = LinearInterp(n, interfaces[k].xloc, interfaces[k].vp , A.x);
        } else {
            int k = internum4elev[mi]; 
            int n = interfaces[k].NX;
            float vp_grad  = LinearInterp(n, interfaces[k].xloc, interfaces[k].vp_grad , A.x);
            float vp       = LinearInterp(n, interfaces[k].xloc, interfaces[k].vp      , A.x);
            float vp_pow   = LinearInterp(n, interfaces[k].xloc, interfaces[k].vp_pow  , A.x);
            float rho      = LinearInterp(n, interfaces[k].xloc, interfaces[k].rho     , A.x);
            float rho_grad = LinearInterp(n, interfaces[k].xloc, interfaces[k].rho_grad, A.x);
            float rho_pow  = LinearInterp(n, interfaces[k].xloc, interfaces[k].rho_pow , A.x);
            var[0] = rho + pow(dz,rho_pow)*rho_grad;
            var[1] = vp  + pow(dz, vp_pow)* vp_grad;
            // if the point on the interface, average for loc
            if (isEqual(dz, 0.0)) {
                int fi = findFirstGreaterEqualIndex(A.z, elevation)-1;
                if (fi >= 0) {
                    k = internum4elev[fi];
                    n = interfaces[k].NX;
                    vp       = LinearInterp(n, interfaces[k].xloc, interfaces[k].vp      , A.x);
                    vp_grad  = LinearInterp(n, interfaces[k].xloc, interfaces[k].vp_grad , A.x);
                    vp_pow   = LinearInterp(n, interfaces[k].xloc, interfaces[k].vp_pow  , A.x);
                    rho      = LinearInterp(n, interfaces[k].xloc, interfaces[k].rho     , A.x);
                    rho_grad = LinearInterp(n, interfaces[k].xloc, interfaces[k].rho_grad, A.x);
                    rho_pow  = LinearInterp(n, interfaces[k].xloc, interfaces[k].rho_pow , A.x);
                    dz = elevation[fi] - A.z;
                    float vp_top  = vp  + pow(dz,  vp_pow)* vp_grad;
                    float rho_top = rho + pow(dz, rho_pow)* rho_grad;
                    var[0] = (var[0] + rho_top)/2.0;
                    var[1] = (var[1] +  vp_top)/2.0;
                }
            }
        }
    break;

    default: // for self-check
        fprintf(stderr,"Error: Unknow media\n");
        fflush(stderr);
        exit(1);

    }     
} 


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
    float *rho)
{
    for (size_t k = 0; k < nz; ++k) {
        for (size_t i = 0; i < nx; ++i) {
            size_t indx =  i + k*nx;

            size_t indx_z = indx;          // for vmap and curv: z
            size_t indx_x = i; // for vmap and cart: x, y

            if (grid_type == GRID_CART) {
                indx_z = k;                // for cart: z
            } else if (grid_type == GRID_CURV) {
                indx_x = indx;
            }

            std::vector<float> var(7, 0.0); 

            AssignLayerMediaPara2Point(i, k,
                        Point2(Gridx[indx_x], Gridz[indx_z]), 
                        interfaces, media_type, var);

            para2tti(var, media_type, // return cij
                c11[indx], c13[indx], c15[indx],
                c33[indx], c35[indx],
                c55[indx], rho[indx]);
        }
    }
}

//======================== averaging/equivalent medium method ================================

/* 
 * For half-grid point, marked the materials number.  
 */
void MarkInterfaceNumber(
        int grid_type,
        float *Hx, float *Hz,
        size_t nx, size_t nz,
        int *MaterNum, // nx*ny*nz
        inter_t *interfaces) 
{
    size_t NI = interfaces[0].NI;
    size_t siz_slice = nx*nz;
    
    for (int k = 0; k < nz; k++) {
        for (int i = 0; i < nx; i++) {
            int indx = i + k*nx;
            int indx_x = indx, indx_z = indx;
            if (grid_type != GRID_CURV)
                indx_x = i;
            if (grid_type == GRID_CART)
                indx_z = k;

            // for every point, calculation the elevation 
            std::vector<float> elevation;
            std::vector<int> internum4elev;
            for (int ni = 0; ni < NI; ni++) {
                int n = interfaces[ni].NX;
                if (Hx[indx_x] >= interfaces[ni].xloc[0] && Hx[indx_x] <= interfaces[ni].xloc[n-1]) {
                    float elev = LinearInterp(n, interfaces[ni].xloc, interfaces[ni].elevation, Hx[indx_x]);
                    elevation.emplace_back(elev);
                    internum4elev.emplace_back(ni);
                }
            }

            // find the which elevation used
            int mi = findLastGreaterEqualIndex(Hz[indx_z], elevation);
            if (mi == -1) 
                mi = findLastGreaterEqualIndex(elevation[0], elevation);

            MaterNum[indx] = internum4elev[mi]; 
        }
    }
    
}


//- 4.0 Assign the parameter by volume arithmetic and harmonic averaging method 
//- elasic tti
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
    float *rho) 
{

    size_t siz_slice = nz * nx;
    size_t NI = interfaces[0].NI;
   
    // assign the local value first.
    parametrization_layer_el_aniso_loc(nx, nz, Gridx, Gridz,
        grid_type, media_type, interfaces, c11, c13, c15, c33, c35, c55, rho);

    float *Hx = nullptr, *Hz = nullptr; 
    GenerateHalfGrid(nx, nz, grid_type, Gridx, Gridz, &Hx, &Hz); 

    /* mark the interface number at the half grid. */
    int *MaterNum = new int[siz_slice];

    MarkInterfaceNumber(grid_type, Hx, Hz, nx, nz, MaterNum, interfaces);

    for (size_t k = 1; k < nz-1; k++) {
        for (size_t i = 1; i < nx-1; i++) {
            size_t indx =  i + k * nx; 
            /* Check if the corresponding the half-grid mesh have different values */
            std::vector<int> v(4);
            v[0] = MaterNum[indx-1-nx];
            v[1] = MaterNum[indx  -nx];
            v[2] = MaterNum[indx     ];
            v[3] = MaterNum[indx-1   ];

            // There is more than one medium value in the half-grid mesh, 
            if ( NumOfValues(v, NI) > 1) {
                
                Mesh2 M = GenerateHalfMesh(grid_type, i, k, indx, 
                            nx, siz_slice, Hx, Hz);

                // recalculate the material value of the point
                float ari_rho = 0.0;
                float har_c11 = 0.0;
                float har_c13 = 0.0;
                float har_c15 = 0.0;
                float har_c33 = 0.0;
                float har_c35 = 0.0;
                float har_c55 = 0.0;

                Point2 *SubGrid = MeshSubdivide(M);

                int nsg = (NG+1)*(NG+1);
                int num_dis = nsg;
                for (int isg = 0; isg < nsg; isg++) {

                    std::vector<float> var(7, 0.0);
                    int isPointOnInter = AssignLayerMediaPara2Point(i, k,
                        SubGrid[isg], interfaces, media_type, var);

                    if (isPointOnInter == 1) {
                        num_dis--;
                    } else{
                        // for every sub-point, transfer para to cij
                        float c11_p = 0.0, c13_p = 0.0;
                        float c15_p = 0.0, c33_p = 0.0;
                        float c35_p = 0.0, c55_p = 0.0;
                        float rho_p = 0.0;
    
                        para2tti(var, media_type, // return cij
                            c11_p, c13_p, c15_p, 
                            c33_p, c35_p, c55_p, rho_p);
    
                        har_c11 += (1.0/c11_p);
                        har_c13 += (1.0/c13_p);
                        har_c15 += (1.0/c15_p);
                        har_c33 += (1.0/c33_p);
                        har_c35 += (1.0/c35_p);
                        har_c55 += (1.0/c55_p);
                        ari_rho += rho_p;
                    }
                }

                c11[indx] = (1.0*num_dis)/har_c11;
                c13[indx] = (1.0*num_dis)/har_c13;
                c15[indx] = (1.0*num_dis)/har_c15;
                c33[indx] = (1.0*num_dis)/har_c33;
                c35[indx] = (1.0*num_dis)/har_c35;
                c55[indx] = (1.0*num_dis)/har_c55;
                rho[indx] = ari_rho/num_dis;

                if(SubGrid != nullptr) delete []SubGrid;
            }

        }
    }
    

    if (Hx != nullptr) delete [] Hx;
    if (Hz != nullptr) delete [] Hz;
    if (MaterNum != nullptr) delete [] MaterNum;    
}

//- 4.1 Assign the parameter by volume arithmetic averaging method 
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
    float *rho) 
{
    size_t siz_slice = nz * nx;
    size_t NI = interfaces[0].NI;
    
    // assign the local value first.
    parametrization_layer_el_aniso_loc(nx, nz, Gridx, Gridz,
        grid_type, media_type, interfaces, c11, c13, c15, c33, c35, c55, rho);

    float *Hx = nullptr, *Hz = nullptr; 
    GenerateHalfGrid(nx, nz, grid_type, Gridx, Gridz, &Hx, &Hz); 

    /* mark the interface number at the half grid. */
    int *MaterNum = new int[siz_slice];

    MarkInterfaceNumber(grid_type, Hx, Hz, nx, nz, MaterNum, interfaces);

    for (size_t k = 1; k < nz-1; k++) {
        for (size_t i = 1; i < nx-1; i++) {
            size_t indx =  i + k*nx; 

            /* Check if the corresponding the half-grid mesh have different values */
            std::vector<int> v(4);
            v[0] = MaterNum[indx-1-nx];
            v[1] = MaterNum[indx  -nx];
            v[2] = MaterNum[indx     ];
            v[3] = MaterNum[indx-1   ];

            // There is more than one medium value in the half-grid mesh, 
            if ( NumOfValues(v, NI) > 1) {
                
                Mesh2 M = GenerateHalfMesh(grid_type, i, k, indx, 
                        nx, siz_slice, Hx, Hz);

                // recalculate the material value of the point
                float ari_rho = 0.0;
                float ari_c11 = 0.0;
                float ari_c13 = 0.0;
                float ari_c15 = 0.0;
                float ari_c33 = 0.0;
                float ari_c35 = 0.0;
                float ari_c55 = 0.0;

                Point2 *SubGrid = MeshSubdivide(M);

                int nsg = (NG+1)*(NG+1);
                int num_dis = nsg;
                for (int isg = 0; isg < nsg; isg++) {

                    std::vector<float> var(7, 0.0);
                    int isPointOnInter = AssignLayerMediaPara2Point(i, k,
                        SubGrid[isg], interfaces, media_type, var);

                    if (isPointOnInter == 1) {
                        num_dis--;
                    } else {
                        // for every sub-point, transfer para to cij
                        float c11_p = 0.0, c13_p = 0.0;
                        float c15_p = 0.0, c33_p = 0.0;
                        float c35_p = 0.0, c55_p = 0.0;
                        float rho_p = 0.0;
    
                        para2tti(var, media_type, // return cij
                            c11_p, c13_p, c15_p, 
                            c33_p, c35_p, c55_p, rho_p);
    
                        ari_c11 += c11_p;
                        ari_c13 += c13_p;
                        ari_c15 += c15_p;
                        ari_c33 += c33_p;
                        ari_c35 += c35_p;
                        ari_c55 += c55_p;
                        ari_rho += rho_p;
                    }
                }

                c11[indx] = ari_c11/num_dis;
                c13[indx] = ari_c13/num_dis;
                c15[indx] = ari_c15/num_dis;
                c33[indx] = ari_c33/num_dis;
                c35[indx] = ari_c35/num_dis;
                c55[indx] = ari_c55/num_dis;
                
                rho[indx] = ari_rho/num_dis;

                if(SubGrid != nullptr) delete []SubGrid;
            }

        }
    }
    

    if (Hx != nullptr) delete [] Hx;
    if (Hz != nullptr) delete [] Hz;
    if (MaterNum != nullptr) delete [] MaterNum;    
}



void getInterfaceIntersection_layer(
    int npoint, float *x, float *z,        // for interface
    const Point2 &v1, const Point2 &v2,    // the edge
    std::set<Point2> &intersectionList)    // segment(vertex1, vertex2)
{
    for (int i = 0; i < npoint-1; i++) {
        Point2 b1(x[i], z[i]);
        Point2 b2(x[i+1], z[i+1]);
        if (SegmentProperIntersection(v1, v2, b1, b2)) {
            Point2 intersectionP = GetLineIntersection(v1, v2-v1, b1, b2-b1);
            intersectionList.insert(intersectionP);
        } else{
            if (b1 == v1 || b1 == v2) {
                intersectionList.insert(b1);
            } else if (b2 == v1 || b2 == v2)  {
                intersectionList.insert(b2); 
            } else{
                if (isPointOnSegment(b1, v1, v2)) intersectionList.insert(b1);
                if (isPointOnSegment(b2, v1, v2)) intersectionList.insert(b2);
                if (isPointOnSegment(v1, b1, b2)) intersectionList.insert(v1);
                if (isPointOnSegment(v2, b1, b2)) intersectionList.insert(v2);
            }
        }
    }
}


/*
 * The S-M method for horizontally layered anisotropic medium,
 * Ref:
 *  1. Schoenberg and Muir, 1989, A calculus for finely layered anisotropic media
 *  2. Jiang, L. and Zhang, W., 2021. TTI equivalent medium parametrization method for the seismic waveform modelling 
 *     of heterogeneous media with coarse grids. GJI, 227(3), pp.2016-2043. DOI:10.1093/gji/ggab310	
 */
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
    float *rho)
{
    size_t siz_slice = nz * nx;
    size_t NI = interfaces[0].NI;

    // assign the local value by har first 
    // (can just consider the numofvalue == 2).
    parametrization_layer_el_aniso_ari(nx, nz, Gridx, Gridz,
        grid_type, media_type, interfaces, c11, c13, c15, c33, c35, c55, rho);

    float *Hx = nullptr, *Hz = nullptr; 
    GenerateHalfGrid(nx, nz, grid_type, Gridx, Gridz, &Hx, &Hz); 

    /* mark the interface number at the half grid. */
    int *MaterNum = new int[siz_slice];
    MarkInterfaceNumber(grid_type, Hx, Hz, nx, nz, MaterNum, interfaces);

    for (size_t k = 1; k < nz-1; k++) {
      for (size_t i = 1; i < nx-1; i++) {
        size_t indx = i+k*nx;
        std::vector<int> v(4);
        v[0] = MaterNum[indx-1-nx];
        v[1] = MaterNum[indx  -nx];
        v[2] = MaterNum[indx     ];
        v[3] = MaterNum[indx-1   ];

        if (NumOfValues(v, NI) == 2) {
          Mesh2 M = GenerateHalfMesh(grid_type, i, k, indx, 
                  nx, siz_slice, Hx, Hz);

          //- convert value to sign +(>= 0), -(< 0)
          int sum = std::accumulate(v.begin(), v.end(), 0);
          int interNum = -1;
          sum = sum/4+1; 
          for (auto &sign:v) {
              interNum = std::max(interNum, sign);
              sign -= sum;
          }

          int cubeState = 0;
          if (v[0] < 0) cubeState |= 1;
          if (v[1] < 0) cubeState |= 2;
          if (v[2] < 0) cubeState |= 4;
          if (v[3] < 0) cubeState |= 8;

          // get intersection points, use set to avoid duplication
          std::set<Point2> intersectionList;

          if (edgeTable[cubeState] & 1)  // v[0] v[1]
            getInterfaceIntersection_layer(interfaces[interNum].NX,
                interfaces[interNum].xloc, interfaces[interNum].elevation,
                M.v[0], M.v[1], intersectionList);

          if (edgeTable[cubeState] & 2) // v[1] v[2]
            getInterfaceIntersection_layer(interfaces[interNum].NX,
                interfaces[interNum].xloc, interfaces[interNum].elevation,
                M.v[1], M.v[2], intersectionList);

          if (edgeTable[cubeState] & 4) // v[2] v[3]
            getInterfaceIntersection_layer(interfaces[interNum].NX,
                interfaces[interNum].xloc, interfaces[interNum].elevation,
                M.v[2], M.v[3], intersectionList);

          if (edgeTable[cubeState] & 8) // v[3] v[0]
            getInterfaceIntersection_layer(interfaces[interNum].NX,
                interfaces[interNum].xloc, interfaces[interNum].elevation,
                M.v[3], M.v[0], intersectionList);

          // ignore coincidence of the point and multiple intersections 
          // tti equivalent
          if (intersectionList.size() == 2) {
            std::set<Point2>::iterator it = intersectionList.begin();
            Point2 P1 = *it;
            it++;
            Point2 P2 = *it;
            Vector2 x_new = P2 - P1;
            if (P1.x > P2.x) 
                x_new = P1-P2;

            // -pi - pi
            double theta = atan2(x_new.z,x_new.x); 

            // A longer array for initializing matrix 
            float f[4] = {0.0, 0.0, 0.0, 0.0};
            Matrix<float> har_Cnn(2,2,f);
            Matrix<float> Ctn_1st(1,2,f);   // ctn * cnn^-1
            Matrix<float> Ctt_1st(1,1,f);   // Ctt_1
            Matrix<float> Ctt_mid(1,1,f);   // Ctn*Cnn^-1*Ctn^T
            Matrix<float> Ctt_lst(2,1,f);   // Cnn^-1 * Ctn^T

            Point2 *SubGrid = MeshSubdivide(M);

            int nsg = (NG+1)*(NG+1);
            int num_dis = nsg;
            for (int isg = 0; isg < nsg; isg++) {

                std::vector<float> var(7, 0.0); 

                int isPointOnInter =AssignLayerMediaPara2Point(i, k,
                    SubGrid[isg], interfaces, media_type, var);

                if (isPointOnInter == 1)
                    num_dis--;
                else {
                    float c11_p = 0.0, c13_p = 0.0;
                    float c15_p = 0.0, c33_p = 0.0;
                    float c35_p = 0.0, c55_p = 0.0;
                    float rho_p = 0.0;
    
                    para2tti(var, media_type, 
                        c11_p, c13_p, c15_p,
                        c33_p, c35_p, c55_p, rho_p);

                    float cnn_p[4] = {c33_p, c35_p, c35_p, c55_p};
                    float ctn_p[2] = {c13_p, c15_p};
                    Matrix<float> cnn(2,2,cnn_p);
                    Matrix<float> ctn(1,2,ctn_p);
                    Matrix<float> ctt(1,1,&c11_p);
                    
                    har_Cnn = har_Cnn + cnn.inverse2x2();
                    Ctn_1st = Ctn_1st + ctn*cnn.inverse2x2();
                    Ctt_1st = Ctt_1st + ctt;
                    Ctt_mid = Ctt_mid + ctn*cnn.inverse2x2()*ctn.transpose();
                    Ctt_lst = Ctt_lst + cnn.inverse2x2()*ctn.transpose();
                }
            }

            har_Cnn = har_Cnn/(num_dis*1.0);
            Ctn_1st = Ctn_1st/(num_dis*1.0);
            Ctt_1st = Ctt_1st/(num_dis*1.0);
            Ctt_mid = Ctt_mid/(num_dis*1.0);
            Ctt_lst = Ctt_lst/(num_dis*1.0);

            // after averaging
            Matrix<float> CNN(2,2), CTN(1,2), CTT(1,1);
            CNN = har_Cnn.inverse2x2();
            CTN = Ctn_1st * CNN;
            CTT = Ctt_1st - Ctt_mid + CTN*Ctt_lst;

            BondTransform2d(CTT(0,0), CTN(0,0), CTN(0,1), CNN(0,0), CNN(0,1), CNN(1,1), theta,
                c11[indx], c13[indx], c15[indx], c33[indx], c35[indx], c55[indx]);
          }

        } // numofvalue == 2
      }
    }
    
    if (Hx != nullptr) delete [] Hx;
    if (Hz != nullptr) delete [] Hz;
    if (MaterNum != nullptr) delete [] MaterNum;    
}
