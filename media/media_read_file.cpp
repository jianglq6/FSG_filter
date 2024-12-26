#include <iostream>
#include <string.h>
#include <math.h>
#include "media_read_file.hpp"
#include "media_utility.hpp"

FILE *gfopen(const char *filename, const char *mode)
{
    fprintf(stdout," - Reading media data from file: %s\n", filename);
    FILE *fp;
    if ((fp = fopen(filename,mode)) == NULL) {
        fprintf(stderr, "Error: Cannot open %s, " \
            "please check your file path and run-directory.\n",filename);
        exit(1);
    }
    return fp;
}


void read_interface_file(
    const char *interface_file,
    inter_t **interfaces,
    int &md_type)   // interfaces[ni].
{
    char line[MAX_BUF_LEN];
    FILE *file = gfopen(interface_file, "r");
    FILE *tmp_file = tmpfile();
    char media_type[MAX_BUF_LEN];
    int   NI;

    /* Read every line in the file. */
    while(fgets(line, MAX_BUF_LEN, file) != NULL)
    {
        if (line[0] == '#' || line[0] == '\n')
            continue;
        fputs(line,tmp_file);
    } 

    /* Set the file pointer at the beginning of the stream */
    rewind(tmp_file);

    /* Read the temporary data file, Can not use the feof() */
    while(feof(tmp_file) != EOF)
    {
    //- header
        if (fscanf(tmp_file, "%s", media_type) < 1) {
            fprintf(stderr,"Error: Please give a media_type in %s!\n", interface_file);
            fflush(stderr);
            exit(1);
        } 

        if (strcmp(media_type, "one_component") == 0) {
            md_type = ONE_COMPONENT; 
        } else if (strcmp(media_type, "acoustic_isotropic") == 0) {
            md_type = ACOUSTIC_ISOTROPIC; 
        } else if (strcmp(media_type, "elastic_isotropic") == 0) {
            md_type = ELASTIC_ISOTROPIC; 
        } else if (strcmp(media_type, "elastic_vti_prem") == 0) {
            md_type = ELASTIC_VTI_PREM; 
        } else if (strcmp(media_type, "elastic_vti_thomsen") == 0) {
            md_type = ELASTIC_VTI_THOMSEN; 
        } else if (strcmp(media_type, "elastic_vti_cij") == 0) {
            md_type = ELASTIC_VTI_CIJ; 
        } else if (strcmp(media_type, "elastic_tti_thomsen") == 0) {
            md_type = ELASTIC_TTI_THOMSEN; 
        } else if (strcmp(media_type, "elastic_tti_bond") == 0) {
            md_type = ELASTIC_TTI_BOND;
        } else if (strcmp(media_type, "elastic_aniso_cij") == 0) {
            md_type = ELASTIC_ANISO_CIJ; 
        } else {
            fprintf(stderr,"Error: media_type=%s is not supported, \n"\
                           "       please check %s!\n", media_type, interface_file);
            fflush(stderr);
            exit(1);
        }

    //- how many interfaces
        if (fscanf(tmp_file, "%d", &NI) < 1) {
            fprintf(stderr,"Error: please give a number of layers in %s!\n", interface_file);
            fflush(stderr);
            exit(1);
        }

        if (NI < 1) {
            fprintf(stderr, "Error: No enough interfaes (minimum is 1), please check file %s!\n", interface_file);
            fflush(stderr);
            exit(1);
        }
        
        (*interfaces) = new inter_t [NI];
        
        for (int ni = 0; ni < NI; ni++) {

            (*interfaces)[ni].NI = NI;
    //- how many points in ni-layer
            int npoint = 0;
            if (fscanf(tmp_file, "%d", &npoint) < 1) {
                fprintf(stderr,"Error: please check the given interfaces mesh in %s!\n", interface_file);
                fflush(stderr);
                exit(1);
            }

            if (npoint < 2) {
                fprintf(stderr, "Error: No enough point (minimum is 2), please check file %s!\n", interface_file);
                fflush(stderr);   
                exit(1);         
            }

            (*interfaces)[ni].NX = npoint;
            // for every point, there is a x_loc and the corresponding z_loc (elevation)
            (*interfaces)[ni].xloc = new float [npoint];
            (*interfaces)[ni].elevation = new float [npoint];

    // - read data of every point in #ni-interface
            switch (md_type) 
            {
            case ONE_COMPONENT:
                (*interfaces)[ni].var      = new float[npoint];
                (*interfaces)[ni].var_grad = new float[npoint];
                (*interfaces)[ni].var_pow  = new float[npoint];
                for (int ip = 0; ip < npoint; ip++) {
                    int num_read = 0;
                    num_read = fscanf( tmp_file, "%f %f %f %f %f",
                                       &(*interfaces)[ni].xloc[ip], 
                                       &(*interfaces)[ni].elevation[ip], 
                                       &(*interfaces)[ni].var[ip]      ,
                                       &(*interfaces)[ni].var_grad[ip] ,
                                       &(*interfaces)[ni].var_pow[ip]   );
                    if (num_read < 5) {
                        fprintf(stderr, "Error: Insufficient data in %s.\n", interface_file);
                        fflush(stderr);
                        exit(1);
                    }
                    if (ip > 0 && (*interfaces)[ni].xloc[ip] < (*interfaces)[ni].xloc[ip-1]) {
                        fprintf(stderr, "Error: interface point must be sorted in increasing X, \n" \
                                        "       please check #%d point of #%d interface in %s",
                                        ip, ni, interface_file);
                        fflush(stderr);
                        exit(1);
                    }
                }
            break;

            case ACOUSTIC_ISOTROPIC:
                (*interfaces)[ni].rho      = new float[npoint];
                (*interfaces)[ni].rho_grad = new float[npoint];
                (*interfaces)[ni].rho_pow  = new float[npoint];
                (*interfaces)[ni].vp       = new float[npoint];
                (*interfaces)[ni].vp_grad  = new float[npoint];
                (*interfaces)[ni].vp_pow   = new float[npoint];
                for (int ip = 0; ip < npoint; ip++) {
                    int num_read = 0;
                    num_read = fscanf( tmp_file, "%f %f %f %f %f %f %f %f",
                                       &(*interfaces)[ni].xloc[ip]     , 
                                       &(*interfaces)[ni].elevation[ip], 
                                       &(*interfaces)[ni].rho[ip]      ,
                                       &(*interfaces)[ni].rho_grad[ip] ,
                                       &(*interfaces)[ni].rho_pow[ip]  , 
                                       &(*interfaces)[ni].vp[ip]       ,
                                       &(*interfaces)[ni].vp_grad[ip]  ,
                                       &(*interfaces)[ni].vp_pow[ip]    );
                    if (num_read < 8) {
                        fprintf(stderr, "Error: Insufficient data in %s.\n", interface_file);
                        fflush(stderr);
                        exit(1);
                    }
                    if (ip > 0 && (*interfaces)[ni].xloc[ip] <= (*interfaces)[ni].xloc[ip-1]) {
                        fprintf(stderr, "Error: interface point must be sorted in increasing X, \n" \
                                        "       please check #%d point of #%d interface in %s",
                                        ip, ni, interface_file);
                        fflush(stderr);
                        exit(1);
                    }
                }
            break;

            case ELASTIC_ISOTROPIC:
                (*interfaces)[ni].rho      = new float[npoint];
                (*interfaces)[ni].rho_grad = new float[npoint];
                (*interfaces)[ni].rho_pow  = new float[npoint];
                (*interfaces)[ni].vp       = new float[npoint];
                (*interfaces)[ni].vp_grad  = new float[npoint];
                (*interfaces)[ni].vp_pow   = new float[npoint];
                (*interfaces)[ni].vs       = new float[npoint];
                (*interfaces)[ni].vs_grad  = new float[npoint];
                (*interfaces)[ni].vs_pow   = new float[npoint];
                for (int ip = 0; ip < npoint; ip++) {
                    int num_read = 0;
                    num_read = fscanf( tmp_file, "%f %f %f %f %f %f %f %f %f %f %f",
                                       &(*interfaces)[ni].xloc[ip]     , 
                                       &(*interfaces)[ni].elevation[ip], 
                                       &(*interfaces)[ni].rho[ip]      ,
                                       &(*interfaces)[ni].rho_grad[ip] ,
                                       &(*interfaces)[ni].rho_pow[ip]  , 
                                       &(*interfaces)[ni].vp[ip]       ,
                                       &(*interfaces)[ni].vp_grad[ip]  ,
                                       &(*interfaces)[ni].vp_pow[ip]   ,
                                       &(*interfaces)[ni].vs[ip]       ,
                                       &(*interfaces)[ni].vs_grad[ip]  ,
                                       &(*interfaces)[ni].vs_pow[ip]    );
                    if (num_read < 11) {
                        fprintf(stderr, "Error: Insufficient data in %s.\n", interface_file);
                        fflush(stderr);
                        exit(1);
                    }
                    if (ip > 0 && (*interfaces)[ni].xloc[ip] < (*interfaces)[ni].xloc[ip-1]) {
                        fprintf(stderr, "Error: interface point must be sorted in increasing X, \n" \
                                        "       please check #%d point of #%d interface in %s",
                                        ip, ni, interface_file);
                        fflush(stderr);
                        exit(1);
                    }
                }
            break;

            case ELASTIC_VTI_PREM:
                (*interfaces)[ni].rho      = new float[npoint];
                (*interfaces)[ni].rho_grad = new float[npoint];
                (*interfaces)[ni].rho_pow  = new float[npoint];
                (*interfaces)[ni].vph      = new float[npoint];
                (*interfaces)[ni].vph_grad = new float[npoint];
                (*interfaces)[ni].vph_pow  = new float[npoint];
                (*interfaces)[ni].vpv      = new float[npoint];
                (*interfaces)[ni].vpv_grad = new float[npoint];
                (*interfaces)[ni].vpv_pow  = new float[npoint];
                (*interfaces)[ni].vsv      = new float[npoint];
                (*interfaces)[ni].vsv_grad = new float[npoint];
                (*interfaces)[ni].vsv_pow  = new float[npoint];
                (*interfaces)[ni].eta      = new float[npoint];
                (*interfaces)[ni].eta_grad = new float[npoint];
                (*interfaces)[ni].eta_pow  = new float[npoint];

                for (int ip = 0; ip < npoint; ip++) {
                    int num_read = 0;
                    num_read = fscanf(tmp_file,
                        "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f ",
                        &(*interfaces)[ni].xloc[ip]     , 
                        &(*interfaces)[ni].elevation[ip], 
                        &(*interfaces)[ni].rho[ip]      ,
                        &(*interfaces)[ni].rho_grad[ip] ,
                        &(*interfaces)[ni].rho_pow[ip]  ,
                        &(*interfaces)[ni].vph[ip]      ,
                        &(*interfaces)[ni].vph_grad[ip] ,
                        &(*interfaces)[ni].vph_pow[ip]  ,
                        &(*interfaces)[ni].vpv[ip]      ,
                        &(*interfaces)[ni].vpv_grad[ip] ,
                        &(*interfaces)[ni].vpv_pow[ip]  ,
                        &(*interfaces)[ni].vsv[ip]      ,
                        &(*interfaces)[ni].vsv_grad[ip] ,
                        &(*interfaces)[ni].vsv_pow[ip]  ,
                        &(*interfaces)[ni].eta[ip]      ,
                        &(*interfaces)[ni].eta_grad[ip] ,
                        &(*interfaces)[ni].eta_pow[ip]   );

                    if (num_read < 17) {
                        fprintf(stderr, "Error: Insufficient data in %s.\n", interface_file);
                        fflush(stderr);
                        exit(1);
                    }
                    if (ip > 0 && (*interfaces)[ni].xloc[ip] <= (*interfaces)[ni].xloc[ip-1]) {
                        fprintf(stderr, "Error: interface point must be sorted in increasing X, \n" \
                                        "       please check #%d point of #%d interface in %s",
                                        ip, ni, interface_file);
                        fflush(stderr);
                        exit(1);
                    }
                }
            break;

            case ELASTIC_VTI_THOMSEN:
                (*interfaces)[ni].rho          = new float[npoint];
                (*interfaces)[ni].rho_grad     = new float[npoint];
                (*interfaces)[ni].rho_pow      = new float[npoint];
                (*interfaces)[ni].vp0          = new float[npoint];
                (*interfaces)[ni].vp0_grad     = new float[npoint];
                (*interfaces)[ni].vp0_pow      = new float[npoint];
                (*interfaces)[ni].vs0          = new float[npoint];
                (*interfaces)[ni].vs0_grad     = new float[npoint];
                (*interfaces)[ni].vs0_pow      = new float[npoint];
                (*interfaces)[ni].epsilon      = new float[npoint];
                (*interfaces)[ni].epsilon_grad = new float[npoint];
                (*interfaces)[ni].epsilon_pow  = new float[npoint];
                (*interfaces)[ni].delta        = new float[npoint];
                (*interfaces)[ni].delta_grad   = new float[npoint];
                (*interfaces)[ni].delta_pow    = new float[npoint];

                for (int ip = 0; ip < npoint; ip++) {
                    int num_read = 0;
                    num_read = fscanf(tmp_file, 
                        "%f %f %f %f %f %f %f %f %f %f "\
                        "%f %f %f %f %f %f %f ",
                        &(*interfaces)[ni].xloc[ip]        , 
                        &(*interfaces)[ni].elevation[ip]   , 
                        &(*interfaces)[ni].rho[ip]         ,
                        &(*interfaces)[ni].rho_grad[ip]    ,
                        &(*interfaces)[ni].rho_pow[ip]     ,
                        &(*interfaces)[ni].vp0[ip]         ,
                        &(*interfaces)[ni].vp0_grad[ip]    ,
                        &(*interfaces)[ni].vp0_pow[ip]     ,
                        &(*interfaces)[ni].vs0[ip]         ,
                        &(*interfaces)[ni].vs0_grad[ip]    ,
                        &(*interfaces)[ni].vs0_pow[ip]     ,
                        &(*interfaces)[ni].epsilon[ip]     ,
                        &(*interfaces)[ni].epsilon_grad[ip],
                        &(*interfaces)[ni].epsilon_pow[ip] ,
                        &(*interfaces)[ni].delta[ip]       ,
                        &(*interfaces)[ni].delta_grad[ip]  ,
                        &(*interfaces)[ni].delta_pow[ip]   );

                    if (num_read < 17) {
                        fprintf(stderr, "Error: Insufficient data in %s.\n", interface_file);
                        fflush(stderr);
                        exit(1);
                    }
                    if (ip > 0 && (*interfaces)[ni].xloc[ip] <= (*interfaces)[ni].xloc[ip-1]) {
                        fprintf(stderr, "Error: interface point must be sorted in increasing X, \n" \
                                        "       please check #%d point of #%d interface in %s",
                                        ip, ni, interface_file);
                        fflush(stderr);
                        exit(1);
                    }
                }
            break;

            case ELASTIC_VTI_CIJ:
                (*interfaces)[ni].rho      = new float[npoint];
                (*interfaces)[ni].rho_grad = new float[npoint];
                (*interfaces)[ni].rho_pow  = new float[npoint];
                (*interfaces)[ni].c11      = new float[npoint];
                (*interfaces)[ni].c11_grad = new float[npoint];
                (*interfaces)[ni].c11_pow  = new float[npoint];
                (*interfaces)[ni].c33      = new float[npoint];
                (*interfaces)[ni].c33_grad = new float[npoint];
                (*interfaces)[ni].c33_pow  = new float[npoint];
                (*interfaces)[ni].c55      = new float[npoint];
                (*interfaces)[ni].c55_grad = new float[npoint];
                (*interfaces)[ni].c55_pow  = new float[npoint];
                (*interfaces)[ni].c13      = new float[npoint];
                (*interfaces)[ni].c13_grad = new float[npoint];
                (*interfaces)[ni].c13_pow  = new float[npoint];

                 for (int ip = 0; ip < npoint; ip++) {
                    int num_read = 0;
                    num_read = fscanf( tmp_file, 
                        "%f %f %f %f %f %f %f %f %f %f "\
                        "%f %f %f %f %f %f %f",
                         &(*interfaces)[ni].xloc[ip], 
                         &(*interfaces)[ni].elevation[ip], 
                         &(*interfaces)[ni].rho[ip]      ,
                         &(*interfaces)[ni].rho_grad[ip] ,
                         &(*interfaces)[ni].rho_pow[ip]  ,
                         &(*interfaces)[ni].c11[ip]      ,
                         &(*interfaces)[ni].c11_grad[ip] ,
                         &(*interfaces)[ni].c11_pow[ip]  ,
                         &(*interfaces)[ni].c33[ip]      ,
                         &(*interfaces)[ni].c33_grad[ip] ,
                         &(*interfaces)[ni].c33_pow[ip]  ,
                         &(*interfaces)[ni].c55[ip]      ,
                         &(*interfaces)[ni].c55_grad[ip] ,
                         &(*interfaces)[ni].c55_pow[ip]  ,
                         &(*interfaces)[ni].c13[ip]      ,
                         &(*interfaces)[ni].c13_grad[ip] ,
                         &(*interfaces)[ni].c13_pow[ip]   );

                    if (num_read < 17) {
                        fprintf(stderr, "Error: Insufficient data in %s.\n", interface_file);
                        fflush(stderr);
                        exit(1);
                    }
                    if (ip > 0 && (*interfaces)[ni].xloc[ip] <= (*interfaces)[ni].xloc[ip-1]) {
                        fprintf(stderr, "Error: interface point must be sorted in increasing X, \n" \
                                        "       please check #%d point of #%d interface in %s",
                                        ip, ni, interface_file);
                        fflush(stderr);
                        exit(1);
                    }                    
                }
            break;

            case ELASTIC_TTI_THOMSEN:
                (*interfaces)[ni].rho          = new float[npoint];
                (*interfaces)[ni].rho_grad     = new float[npoint];
                (*interfaces)[ni].rho_pow      = new float[npoint];
                (*interfaces)[ni].vp0          = new float[npoint];
                (*interfaces)[ni].vp0_grad     = new float[npoint];
                (*interfaces)[ni].vp0_pow      = new float[npoint];
                (*interfaces)[ni].vs0          = new float[npoint];
                (*interfaces)[ni].vs0_grad     = new float[npoint];
                (*interfaces)[ni].vs0_pow      = new float[npoint];
                (*interfaces)[ni].epsilon      = new float[npoint];
                (*interfaces)[ni].epsilon_grad = new float[npoint];
                (*interfaces)[ni].epsilon_pow  = new float[npoint];
                (*interfaces)[ni].delta        = new float[npoint];
                (*interfaces)[ni].delta_grad   = new float[npoint];
                (*interfaces)[ni].delta_pow    = new float[npoint];
                (*interfaces)[ni].dip          = new float[npoint];
                (*interfaces)[ni].dip_grad     = new float[npoint];
                (*interfaces)[ni].dip_pow      = new float[npoint];

                for (int ip = 0; ip < npoint; ip++) {
                    int num_read = 0;
                    num_read = fscanf( tmp_file, 
                        "%f %f %f %f %f %f %f %f %f %f " \
                        "%f %f %f %f %f %f %f %f %f %f ",
                        &(*interfaces)[ni].xloc[ip]        , 
                        &(*interfaces)[ni].elevation[ip]   , 
                        &(*interfaces)[ni].rho[ip]         ,
                        &(*interfaces)[ni].rho_grad[ip]    ,
                        &(*interfaces)[ni].rho_pow[ip]     ,
                        &(*interfaces)[ni].vp0[ip]         ,
                        &(*interfaces)[ni].vp0_grad[ip]    ,
                        &(*interfaces)[ni].vp0_pow[ip]     ,
                        &(*interfaces)[ni].vs0[ip]         ,
                        &(*interfaces)[ni].vs0_grad[ip]    ,
                        &(*interfaces)[ni].vs0_pow[ip]     ,  
                        &(*interfaces)[ni].epsilon[ip]     ,
                        &(*interfaces)[ni].epsilon_grad[ip],
                        &(*interfaces)[ni].epsilon_pow[ip] ,
                        &(*interfaces)[ni].delta[ip]       ,
                        &(*interfaces)[ni].delta_grad[ip]  ,
                        &(*interfaces)[ni].delta_pow[ip]   ,
                        &(*interfaces)[ni].dip[ip]         ,
                        &(*interfaces)[ni].dip_grad[ip]    ,
                        &(*interfaces)[ni].dip_pow[ip]      );

                    if (num_read < 20) {
                        fprintf(stderr, "Error: Insufficient data in %s.\n", interface_file);
                        fflush(stderr);
                        exit(1);
                    }
                    
                    if (ip > 0 && (*interfaces)[ni].xloc[ip] <= (*interfaces)[ni].xloc[ip-1]) {
                        fprintf(stderr, "Error: interface point must be sorted in increasing X, \n" \
                                        "       please check #%d point of #%d interface in %s",
                                        ip, ni, interface_file);
                        fflush(stderr);
                        exit(1);
                    }
                }
            break;

            case ELASTIC_TTI_BOND:
                (*interfaces)[ni].rho      = new float[npoint];
                (*interfaces)[ni].rho_grad = new float[npoint];
                (*interfaces)[ni].rho_pow  = new float[npoint];
                (*interfaces)[ni].c11      = new float[npoint];
                (*interfaces)[ni].c11_grad = new float[npoint];
                (*interfaces)[ni].c11_pow  = new float[npoint];
                (*interfaces)[ni].c33      = new float[npoint];
                (*interfaces)[ni].c33_grad = new float[npoint];
                (*interfaces)[ni].c33_pow  = new float[npoint];
                (*interfaces)[ni].c55      = new float[npoint];
                (*interfaces)[ni].c55_grad = new float[npoint];
                (*interfaces)[ni].c55_pow  = new float[npoint];
                (*interfaces)[ni].c13      = new float[npoint];
                (*interfaces)[ni].c13_grad = new float[npoint];
                (*interfaces)[ni].c13_pow  = new float[npoint];
                (*interfaces)[ni].dip      = new float[npoint];
                (*interfaces)[ni].dip_grad = new float[npoint];
                (*interfaces)[ni].dip_pow  = new float[npoint];

                for (int ip = 0; ip < npoint; ip++) {
                    int num_read = 0;
                    num_read = fscanf( tmp_file, 
                        "%f %f %f %f %f %f %f %f %f %f "\
                        "%f %f %f %f %f %f %f %f %f %f",
                        &(*interfaces)[ni].xloc[ip]     , 
                        &(*interfaces)[ni].elevation[ip], 
                        &(*interfaces)[ni].rho[ip]      ,
                        &(*interfaces)[ni].rho_grad[ip] ,
                        &(*interfaces)[ni].rho_pow[ip]  ,
                        &(*interfaces)[ni].c11[ip]      ,
                        &(*interfaces)[ni].c11_grad[ip] ,
                        &(*interfaces)[ni].c11_pow[ip]  ,
                        &(*interfaces)[ni].c33[ip]      ,
                        &(*interfaces)[ni].c33_grad[ip] ,
                        &(*interfaces)[ni].c33_pow[ip]  ,
                        &(*interfaces)[ni].c55[ip]      ,
                        &(*interfaces)[ni].c55_grad[ip] ,
                        &(*interfaces)[ni].c55_pow[ip]  ,
                        &(*interfaces)[ni].c13[ip]      ,
                        &(*interfaces)[ni].c13_grad[ip] ,
                        &(*interfaces)[ni].c13_pow[ip]  ,
                        &(*interfaces)[ni].dip[ip]      ,
                        &(*interfaces)[ni].dip_grad[ip] ,
                        &(*interfaces)[ni].dip_pow[ip]  );

                    if (num_read < 20) {
                        fprintf(stderr, "Error: Insufficient data in %s.\n", interface_file);
                        fflush(stderr);
                        exit(1);
                    }
                    if (ip > 0 && (*interfaces)[ni].xloc[ip] <= (*interfaces)[ni].xloc[ip-1]) {
                        fprintf(stderr, "Error: interface point must be sorted in increasing X, \n" \
                                        "       please check #%d point of #%d interface in %s",
                                        ip, ni, interface_file);
                        fflush(stderr);
                        exit(1);
                    }                   
                
                }
            break;

            case ELASTIC_ANISO_CIJ:
                (*interfaces)[ni].rho      = new float[npoint];
                (*interfaces)[ni].rho_grad = new float[npoint];
                (*interfaces)[ni].rho_pow  = new float[npoint];
                (*interfaces)[ni].c11      = new float[npoint]; 
                (*interfaces)[ni].c11_grad = new float[npoint];      
                (*interfaces)[ni].c11_pow  = new float[npoint];  
                (*interfaces)[ni].c13      = new float[npoint]; 
                (*interfaces)[ni].c13_grad = new float[npoint];  
                (*interfaces)[ni].c13_pow  = new float[npoint];      
                (*interfaces)[ni].c15      = new float[npoint]; 
                (*interfaces)[ni].c15_grad = new float[npoint];     
                (*interfaces)[ni].c15_pow  = new float[npoint];     
                (*interfaces)[ni].c33      = new float[npoint];                 
                (*interfaces)[ni].c33_grad = new float[npoint]; 
                (*interfaces)[ni].c33_pow  = new float[npoint]; 
                (*interfaces)[ni].c35      = new float[npoint];                 
                (*interfaces)[ni].c35_grad = new float[npoint]; 
                (*interfaces)[ni].c35_pow  = new float[npoint]; 
                (*interfaces)[ni].c55      = new float[npoint];
                (*interfaces)[ni].c55_grad = new float[npoint]; 
                (*interfaces)[ni].c55_pow  = new float[npoint];

                for (int ip = 0; ip < npoint; ip++) {
                    int num_read = 0;
                    num_read = fscanf( tmp_file, "%f %f "\
                        "%f %f %f %f %f %f %f "\
                        "%f %f %f %f %f %f %f "\
                        "%f %f %f %f %f %f %f",
                        &(*interfaces)[ni].xloc[ip],
                        &(*interfaces)[ni].elevation[ip],  
                        &(*interfaces)[ni].rho[ip], 
                        &(*interfaces)[ni].rho_grad[ip], 
                        &(*interfaces)[ni].rho_pow[ip],
                        &(*interfaces)[ni].c11[ip], 
                        &(*interfaces)[ni].c11_grad[ip], 
                        &(*interfaces)[ni].c11_pow[ip],
                        &(*interfaces)[ni].c13[ip], 
                        &(*interfaces)[ni].c13_grad[ip], 
                        &(*interfaces)[ni].c13_pow[ip],
                        &(*interfaces)[ni].c15[ip], 
                        &(*interfaces)[ni].c15_grad[ip], 
                        &(*interfaces)[ni].c15_pow[ip],
                        &(*interfaces)[ni].c33[ip], 
                        &(*interfaces)[ni].c33_grad[ip], 
                        &(*interfaces)[ni].c33_pow[ip],
                        &(*interfaces)[ni].c35[ip], 
                        &(*interfaces)[ni].c35_grad[ip], 
                        &(*interfaces)[ni].c35_pow[ip],
                        &(*interfaces)[ni].c55[ip], 
                        &(*interfaces)[ni].c55_grad[ip], 
                        &(*interfaces)[ni].c55_pow[ip] );

                    if (num_read < 23) {
                        fprintf(stderr, "Error: Insufficient data in %s.\n", interface_file);
                        fflush(stderr);
                        exit(1);
                    }
                    if (ip > 0 && (*interfaces)[ni].xloc[ip] <= (*interfaces)[ni].xloc[ip-1]) {
                        fprintf(stderr, "Error: interface point must be sorted in increasing X, \n" \
                                        "       please check #%d point of #%d interface in %s",
                                        ip, ni, interface_file);
                        fflush(stderr);
                        exit(1);
                    }
                }
            break;
            default: // for self-check
                fprintf(stderr,"Error: Unknow media #%d\n", 
                       md_type);
                fflush(stderr);
                exit(1);

            }
        }
        // end of the info
        break;
    } 
    fclose(file);   
    fclose(tmp_file);   
}


/* 
 * Just read the grid data within the given
 *  [Xmin, Xmax]\times[Ymin, Ymax] domain.
 */
void read_grid_file(
    const char *grid_file,
    float Xmin, float Xmax,
    int &NL,
    std::vector<int> &NGz, // how many z-grid in each layer
    inter_t *interfaces)
{
    FILE *file = gfopen(grid_file, "r");
    FILE *tmp_file = tmpfile();

    char  line[MAX_BUF_LEN];
    char  media_type[MAX_BUF_LEN];
    size_t NI = 0;
    int    NX = 0; 
    float  DX = 0;
    float  MINX = 0;

    // remove the annotation
    while(fgets(line, MAX_BUF_LEN, file) != NULL)
    {
        if (line[0] == '#' || line[0] == '\n')
            continue;
        fputs(line,tmp_file);
    } 

    rewind(tmp_file);

    while(feof(tmp_file) != EOF)
    {
    //- header: media type
        if (fscanf(tmp_file, "%s", media_type) < 1) {
            fprintf(stderr,"Error: The GRID MEIDA FILE is wrong, please give a media_type!\n");
            fflush(stderr);
            exit(1);
        } 

        if (strcmp(media_type, "one_component") == 0) {
            interfaces->media_type = ONE_COMPONENT; 
        } else if (strcmp(media_type, "acoustic_isotropic") == 0) {
            interfaces->media_type = ACOUSTIC_ISOTROPIC; 
        } else if (strcmp(media_type, "elastic_isotropic") == 0) {
            interfaces->media_type = ELASTIC_ISOTROPIC; 
        } else if (strcmp(media_type, "elastic_vti_prem") == 0) {
            interfaces->media_type = ELASTIC_VTI_PREM; 
        } else if (strcmp(media_type, "elastic_vti_thomsen") == 0) {
            interfaces->media_type = ELASTIC_VTI_THOMSEN; 
        } else if (strcmp(media_type, "elastic_vti_cij") == 0) {
            interfaces->media_type = ELASTIC_VTI_CIJ; 
        } else if (strcmp(media_type, "elastic_tti_thomsen") == 0) {
            interfaces->media_type = ELASTIC_TTI_THOMSEN; 
        } else if (strcmp(media_type, "elastic_tti_bond") == 0) {
            interfaces->media_type = ELASTIC_TTI_BOND;
        } else if (strcmp(media_type, "elastic_aniso_cij") == 0) {
            interfaces->media_type = ELASTIC_ANISO_CIJ; 
        } else {
            fprintf(stderr,"Error: media_type = %s is not supported, \n"\
                           "       please check %s!\n", 
                           media_type, grid_file);
            fflush(stderr);
            exit(1);
        }

        if (fscanf(tmp_file, "%d", &NL) < 1) {
            fprintf(stderr,"Error: The GRID MEIDA FILE is wrong, " \
                "please give a number of layers! \n");
            fflush(stderr);
            exit(1);
        }

        if (NL < 1) {
            fprintf(stderr, "Error: No enough layers (minimum is 1)!\n");
            fflush(stderr);
            exit(1);
        }

        for (int i = 0; i < NL; i++) {
            int ng_i = 0;
            if (fscanf(tmp_file, "%d",&ng_i) < 1) {
                fprintf(stderr,"Error: please give the number of grids in the %d-th layer! \n", i);
                fflush(stderr);
                exit(1);
            }
            NI += ng_i;
            NGz.push_back(ng_i);
        }

        if (fscanf(tmp_file, "%d %f %f", &NX, &MINX, &DX) < 3) {
            fprintf(stderr,"Error: please check the given interfaces mesh in %s! \n", grid_file);
            fflush(stderr);
            exit(1);
        }

        if (NX < 2) {
            fprintf(stderr,"Error: No enough point (NX >= 2)! \n");
            fflush(stderr); 
            exit(1);           
        }

        /* the range need to read */
        size_t ix0 = (Xmin-MINX)/DX;
        size_t ix1 = ceil((Xmax-MINX)/DX);
        size_t nx  = ix1-ix0+1;
        float minx = ix0 * DX + MINX;

        /* The given media domain must bigger than the calculation domain */
        if (ix0 < 0 || ix1 > NX ) {
            fprintf(stderr,"Error: The given media range is smaller than "\
                "the calculation grid range in x-direction! \n");
            fflush(stderr); 
            exit(1); 
        }

        interfaces -> NI = NI;
        interfaces -> NX = nx;
        interfaces -> DX = DX; 
        interfaces -> MINX = minx;

        // for read interface info
        size_t inter_line = nx;
        size_t inter_slice = inter_line*NI;

        interfaces->elevation = new float[inter_slice];

        switch (interfaces -> media_type) 
        {
        case ONE_COMPONENT:
            interfaces->var = new float[inter_slice];
            for (size_t ni = 0; ni < NI; ni++) {
                for (size_t ix = 0; ix < NX; ix++) {
                    if (ix >= ix0 && ix <= ix1) {
                        size_t ip = (ix-ix0) + ni*inter_line;
                        int num_read = 0;
                        num_read = fscanf( tmp_file, "%f %f",
                                           &(interfaces->elevation[ip]), 
                                           &(interfaces->var[ip]      ));
                        if (num_read < 2) {
                            fprintf(stderr, "Error: Insufficient data in %s.\n", grid_file);
                            fflush(stderr);
                            exit(1);
                        } 
                    } else { // read the data no used.
                        float tmp[2];
                        int u = fscanf(tmp_file, "%f %f", tmp, tmp+1);
                        if (u < 2) {
                            fprintf(stderr,"Error: Insufficient data in %s.\n", grid_file);
                            fflush(stderr);
                            exit(1);
                        }
                    
                    }
                }
            }
        break;

        case ACOUSTIC_ISOTROPIC:
            interfaces->rho = new float[inter_slice];
            interfaces->vp  = new float[inter_slice];
            for (size_t ni = 0; ni < NI; ni++) {
                for (size_t ix = 0; ix < NX; ix++) {
                    if (ix >= ix0 && ix <= ix1) {
                        size_t ip = (ix-ix0) + ni*inter_line;
                        int num_read = 0;
                        num_read = fscanf( tmp_file, "%f %f %f",
                                          &(interfaces->elevation[ip]), 
                                          &(interfaces->rho[ip]      ),
                                          &(interfaces->vp[ip]       ));
                        if (num_read < 3) {
                            fprintf(stderr, "Error: Insufficient data in %s.\n", grid_file);
                            fflush(stderr);
                            exit(1);
                        } 
                    } else { // read the data no used.
                        float tmp[3];
                        int u = fscanf(tmp_file, "%f %f %f", tmp, tmp+1, tmp+2);
                        if (u < 3) {
                            fprintf(stderr,"Error: Insufficient data in %s.\n", grid_file);
                            fflush(stderr);
                            exit(1);
                        }
                    
                    }
                }
            }
        break;

        case ELASTIC_ISOTROPIC:
            interfaces->rho = new float[inter_slice];
            interfaces->vp  = new float[inter_slice];
            interfaces->vs  = new float[inter_slice];
            for (size_t ni = 0; ni < NI; ni++) {
                for (size_t ix = 0; ix < NX; ix++) {
                    if (ix >= ix0 && ix <= ix1) {
                        size_t ip = (ix-ix0) + ni * inter_line;
                        int num_read = 0;
                        num_read = fscanf( tmp_file, "%f %f %f %f",
                                          &(interfaces->elevation[ip]), 
                                          &(interfaces->rho[ip]      ),
                                          &(interfaces->vp[ip]       ),
                                          &(interfaces->vs[ip]       ));
                        if (num_read < 4) {
                            fprintf(stderr, "Error: Insufficient data in %s.\n", grid_file);
                            fflush(stderr);
                            exit(1);
                        } 
                    } else { // read the data no used.
                        float tmp[4];
                        int u = fscanf(tmp_file, "%f %f %f %f", tmp, tmp+1, tmp+2, tmp+3);
                        if (u < 4) {
                            fprintf(stderr,"Error: Insufficient data in %s.\n", grid_file);
                            fflush(stderr);
                            exit(1);
                        }
                            
                    }
                }
            }
        break;

        case ELASTIC_VTI_PREM:
            interfaces->rho = new float[inter_slice];
            interfaces->vph = new float[inter_slice];
            interfaces->vpv = new float[inter_slice];
            interfaces->vsv = new float[inter_slice];
            interfaces->eta = new float[inter_slice];

            for (size_t ni = 0; ni < NI; ni++) {
                for (size_t ix = 0; ix < NX; ix++) {
                    if (ix >= ix0 && ix <= ix1) {
                        size_t ip = (ix-ix0) + ni * inter_line;
                        int num_read = 0;
                        num_read = fscanf( tmp_file, "%f %f %f %f %f %f ",
                                          &(interfaces->elevation[ip]), 
                                          &(interfaces->rho[ip]      ),
                                          &(interfaces->vph[ip]      ),
                                          &(interfaces->vpv[ip]      ),
                                          &(interfaces->vsv[ip]      ),
                                          &(interfaces->eta[ip]      ) );
                        if (num_read < 6) {
                            fprintf(stderr, "Error: Insufficient data in %s.\n", grid_file);
                            fflush(stderr);
                            exit(1);
                        } 
                    
                    } else { // read the data no used.
                        float tmp[6];
                        int u = fscanf(tmp_file, "%f %f %f %f %f %f", 
                            tmp, tmp+1, tmp+2, tmp+3, tmp+4, tmp+5);
                        if (u < 6) {
                            fprintf(stderr, "Error: Insufficient data in %s.\n", grid_file);
                            fflush(stderr);
                            exit(1);
                        }
                    }

                }
            }
        break;

        case ELASTIC_VTI_THOMSEN:
            interfaces->rho = new float[inter_slice];
            interfaces->vp0 = new float[inter_slice];
            interfaces->vs0 = new float[inter_slice];
            interfaces->epsilon = new float[inter_slice];
            interfaces->delta = new float[inter_slice];

            for (size_t ni = 0; ni < NI; ni++) {
                for (size_t ix = 0; ix < NX; ix++) {
                    if (ix >= ix0 && ix <= ix1) {
                        size_t ip = (ix-ix0) + ni * inter_line;
                        int num_read = 0;
                        num_read = fscanf( tmp_file, "%f %f %f %f %f %f ",
                                          &(interfaces->elevation[ip]), 
                                          &(interfaces->rho[ip]      ),
                                          &(interfaces->vp0[ip]      ),
                                          &(interfaces->vs0[ip]      ),
                                          &(interfaces->epsilon[ip]  ),
                                          &(interfaces->delta[ip]    ) );

                        if (num_read < 6) {
                            fprintf(stderr, "Error: Insufficient data in %s.\n", grid_file);
                            fflush(stderr);
                            exit(1);
                        }
                    } else { // read the data no used.
                        float tmp[6];
                        int u = fscanf(tmp_file, "%f %f %f %f %f %f ", 
                            tmp, tmp+1, tmp+2, tmp+3, tmp+4, tmp+5);
                        if (u < 6) {
                            fprintf(stderr,"Error: Insufficient data in %s.\n", grid_file);
                            fflush(stderr);
                            exit(1);
                        }
                    }
                }
            }
        break;

        case ELASTIC_VTI_CIJ:
            interfaces->rho = new float[inter_slice];
            interfaces->c11 = new float[inter_slice];
            interfaces->c33 = new float[inter_slice];
            interfaces->c55 = new float[inter_slice];
            interfaces->c13 = new float[inter_slice];

            for (size_t ni = 0; ni < NI; ni++) {
                for (size_t ix = 0; ix < NX; ix++) {
                    if (ix >= ix0 && ix <= ix1) {
                        size_t ip = (ix-ix0) + ni * inter_line;
                        int num_read = 0;
                        num_read = fscanf( tmp_file, "%f %f %f %f %f %f ",
                                          &(interfaces->elevation[ip]), 
                                          &(interfaces->rho[ip]      ),
                                          &(interfaces->c11[ip]      ),
                                          &(interfaces->c33[ip]      ),
                                          &(interfaces->c55[ip]      ),
                                          &(interfaces->c13[ip]      ) );
                        if (num_read < 6) {
                            fprintf(stderr, "Error: Insufficient data in %s.\n", grid_file);
                            fflush(stderr);
                            exit(1);
                        }
                    
                    } else { // read the data no used.
                        float tmp[6];
                        int u = fscanf(tmp_file, "%f %f %f %f %f %f", 
                            tmp, tmp+1, tmp+2, tmp+3, tmp+4, tmp+5);
                        if (u < 6) {
                            fprintf(stderr, "Error: Insufficient data in %s.\n", grid_file);
                            fflush(stderr);
                            exit(1);
                        }
                    }
                }
            }
        break; 

        case ELASTIC_TTI_THOMSEN:
            interfaces->rho     = new float[inter_slice];
            interfaces->vp0     = new float[inter_slice];
            interfaces->vs0     = new float[inter_slice];
            interfaces->epsilon = new float[inter_slice];
            interfaces->delta   = new float[inter_slice];
            interfaces->dip     = new float[inter_slice];

            for (size_t ni = 0; ni < NI; ni++) {
                for (size_t ix = 0; ix < NX; ix++) {
                    if (ix >= ix0 && ix <= ix1) {
                        size_t ip = (ix-ix0) + ni * inter_line;
                        int num_read = 0;
                        num_read = fscanf( tmp_file, "%f %f %f %f %f %f %f ",
                                          &(interfaces->elevation[ip]), 
                                          &(interfaces->rho[ip]      ),
                                          &(interfaces->vp0[ip]      ),
                                          &(interfaces->vs0[ip]      ),
                                          &(interfaces->epsilon[ip]  ),
                                          &(interfaces->delta[ip]    ),
                                          &(interfaces->dip[ip]      ) );
                        if (num_read < 7) {
                            fprintf(stderr, "Error: Insufficient data in %s.\n", grid_file);
                            fflush(stderr);
                            exit(1);
                        } 
                    } else { // read the data no used.
                        float tmp[7];
                        int u = fscanf(tmp_file, "%f %f %f %f %f %f %f ", 
                            tmp, tmp+1, tmp+2, tmp+3, tmp+4, tmp+5, tmp+6);
                        if (u < 7) {
                            fprintf(stderr,"Error: Insufficient data in %s.\n", grid_file);
                            fflush(stderr);
                            exit(1);
                        }
                    }
                }
            }
        break;               

        case ELASTIC_TTI_BOND:
            interfaces->rho = new float[inter_slice];
            interfaces->c11 = new float[inter_slice];
            interfaces->c33 = new float[inter_slice];
            interfaces->c55 = new float[inter_slice];
            interfaces->c13 = new float[inter_slice];
            interfaces->dip = new float[inter_slice];

            for (size_t ni = 0; ni < NI; ni++) {
                for (size_t ix = 0; ix < NX; ix++) {
                    if (ix >= ix0 && ix <= ix1) {
                        size_t ip = (ix-ix0) + ni * inter_line;
                        int num_read = 0;
                        num_read = fscanf( tmp_file, "%f %f %f %f %f %f %f ",
                                          &(interfaces->elevation[ip]), 
                                          &(interfaces->rho[ip]      ),
                                          &(interfaces->c11[ip]      ),
                                          &(interfaces->c33[ip]      ),
                                          &(interfaces->c55[ip]      ),
                                          &(interfaces->c13[ip]      ),
                                          &(interfaces->dip[ip]      ) );
                        if (num_read < 7) {
                            fprintf(stderr, "Error: Insufficient data in %s.\n", grid_file);
                            fflush(stderr);
                            exit(1);
                        } 
    
                    } else { // read the data no used.
                        float tmp[7];
                        int u = fscanf(tmp_file, "%f %f %f %f %f %f %f ", 
                            tmp, tmp+1, tmp+2, tmp+3, tmp+4, tmp+5, tmp+6);
                        if (u < 7) {
                            fprintf(stderr,"Error: Insufficient data in %s.\n", grid_file);
                            fflush(stderr);
                            exit(1);
                        }
                    }
                }
            }
        break;

        case ELASTIC_ANISO_CIJ:
            interfaces->rho = new float[inter_slice];
            interfaces->c11 = new float[inter_slice];
            interfaces->c13 = new float[inter_slice];
            interfaces->c15 = new float[inter_slice];
            interfaces->c33 = new float[inter_slice];
            interfaces->c35 = new float[inter_slice];
            interfaces->c55 = new float[inter_slice];

            for (size_t ni = 0; ni < NI; ni++) {
                for (size_t ix = 0; ix < NX; ix++) {
                    if (ix >= ix0 && ix <= ix1) {
                        size_t ip = (ix-ix0) + ni*inter_line;
                        int num_read = 0;
                        num_read = fscanf( tmp_file, "%f %f %f %f %f %f %f %f ",
                                            &(interfaces->elevation[ip]),
                                            &(interfaces->rho[ip]),
                                            &(interfaces->c11[ip]),
                                            &(interfaces->c13[ip]),
                                            &(interfaces->c15[ip]),
                                            &(interfaces->c33[ip]),
                                            &(interfaces->c35[ip]),
                                            &(interfaces->c55[ip]) );

                        if (num_read < 8) {
                            fprintf(stderr, "Error: Insufficient data in %s.\n", grid_file);
                            fflush(stderr);
                            exit(1);
                        } 
                    
                    } else { // read the data no used.
                        float tmp[8];
                        int u = fscanf( tmp_file, "%f %f %f %f %f %f %f %f ",
                                        tmp, tmp+1, tmp+2, tmp+3, tmp+4, tmp+5, tmp+6, tmp+7);
                        if (u < 8) {
                            fprintf(stderr,"Error: Insufficient data in %s.\n", grid_file);
                            fflush(stderr);
                            exit(1);
                        }
                    }
                }
            }
        break;

        default:
            fprintf(stderr, "Error: Unknow media_type %s (for code check)\n", media_type);
            fflush(stderr);
            exit(1);
        }
        break;
    }

    checkGridData(NL, NGz, *interfaces, grid_file);

    fclose(file);   
    fclose(tmp_file);    
}

// check whether the elevation[ng[i]-1] == elevation[ng[i]] 
int checkGridData(int NL, 
    std::vector <int> &NGz,
    inter_t &interfaces, 
    const char *grid_file) 
{   
    if (NL < 2) return 0;

    int nx = interfaces.NX;
    int ipoint = 0;
    for (int i = 0; i < NL-1; i++) {
        ipoint += (NGz[i]); 
        for (int indx = 0; indx < nx; indx++) {
            int indx_top = indx + (ipoint-1)*nx;
            int indx_bot = indx + ipoint*nx;
            if (interfaces.elevation[indx_top] != interfaces.elevation[indx_bot]) {
                fprintf(stderr, "Error: The last elevations of #%d layer should equal to the first "\
                                "elevations of #%d layer!, please check %s!", i, i+1, grid_file);        
                fflush(stderr);
                exit(1);
            }
        }      
    }

    return 0;
}
