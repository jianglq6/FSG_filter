#ifndef _MEDIA_UTILITY_
#define _MEDIA_UTILITY_

#include "media_geometry2d.hpp"
#include <map>
// for equivalent medium parametrization
#define NG 8

#define PI acos(-1)

// for media type
#define ONE_COMPONENT       0  /* var 16*/
#define ACOUSTIC_ISOTROPIC  10 /* rho vp 2*/
#define ELASTIC_ISOTROPIC   21 /* rho vp vs 3*/
#define ELASTIC_VTI_PREM    22 /* rho vph vpv vsv eta 5*/
#define ELASTIC_VTI_THOMSEN 23 /* rho vp0 vs0 epsilon delta 5*/
#define ELASTIC_VTI_CIJ     24 /* rho c11 c33 c55 c13 5*/
#define ELASTIC_TTI_THOMSEN 25 /* rho vp0 vs0 epsilon delta dip 6 */
#define ELASTIC_TTI_BOND    26 /* rho c11 c33 c55 c13 dip 6*/
#define ELASTIC_ANISO_CIJ   27 /* rho c11 c13 c15 c33 c35 c55 7*/

#define GRID_CART 1
#define GRID_VMAP 2
#define GRID_CURV 3

// for report error: int -> string
std::map<int, std::string> create_md2str_map(); 

/* 
 * number used for divided the mesh, 
 * for equivalent medium parameterization
 */

/* Interface information from the interface file (different media type) */
struct inter_t{
    int media_type = -1;
    /* all the interfaces are given by the same interface_file mesh. */    
    int   NI = 0; // number of interfaces (layer2model), sum of ngz (grid2model)
    int   NX = 0; // = npoint in layer2model
    float DX = FLT_MAX;
    float MINX = FLT_MAX;

    // ni*slice
    float *elevation= nullptr;
    float *xloc = nullptr;

    // for acoustic + elastic
    float *vp       = nullptr;
    float *rho      = nullptr;
    float *vp_grad  = nullptr;
    float *rho_grad = nullptr;
    float *vp_pow   = nullptr;
    float *rho_pow  = nullptr;
    float *vs       = nullptr;
    float *vs_grad  = nullptr;
    float *vs_pow   = nullptr;

    // for anisotropy (thomsen)
    float *epsilon      = nullptr; 
    float *epsilon_grad = nullptr; 
    float *epsilon_pow  = nullptr; 
    float *delta        = nullptr;
    float *delta_grad   = nullptr;
    float *delta_pow    = nullptr;
    float *dip          = nullptr;
    float *dip_grad     = nullptr;
    float *dip_pow      = nullptr;

    // for thomsen
    float *vp0      = nullptr;
    float *vs0      = nullptr;
    float *vp0_grad = nullptr;
    float *vs0_grad = nullptr;
    float *vp0_pow  = nullptr;
    float *vs0_pow  = nullptr;

    // for vti_prem
    float *vph      = nullptr;
    float *vpv      = nullptr;
    float *vsv      = nullptr;
    float *eta      = nullptr;
    float *vph_grad = nullptr;
    float *vpv_grad = nullptr;
    float *vsv_grad = nullptr;
    float *eta_grad = nullptr;
    float *vph_pow  = nullptr;
    float *vpv_pow  = nullptr;
    float *vsv_pow  = nullptr;
    float *eta_pow  = nullptr;

    // for one component
    float *var      = nullptr;
    float *var_grad = nullptr;
    float *var_pow  = nullptr;

    // for anisotropy c_ij, call one component
    float *c11 = nullptr;
    float *c13 = nullptr;
    float *c15 = nullptr;
    float *c33 = nullptr;
    float *c35 = nullptr;
    float *c55 = nullptr;

    float *c11_grad = nullptr;
    float *c13_grad = nullptr;
    float *c15_grad = nullptr;
    float *c33_grad = nullptr;
    float *c35_grad = nullptr;
    float *c55_grad = nullptr;
    
    float *c11_pow = nullptr;
    float *c13_pow = nullptr;
    float *c15_pow = nullptr;
    float *c33_pow = nullptr;
    float *c35_pow = nullptr;
    float *c55_pow = nullptr;

    ~inter_t() {
        // ni*slice
        if (elevation != nullptr) delete [] elevation;
        if (xloc != nullptr) delete [] xloc;

        // for one component
        if (var      != nullptr) delete [] var     ; 
        if (var_grad != nullptr) delete [] var_grad;
        if (var_pow  != nullptr) delete [] var_pow ;

        // for acoustic + elastic
        if (vp       != nullptr) delete [] vp      ;
        if (rho      != nullptr) delete [] rho     ;
        if (vp_grad  != nullptr) delete [] vp_grad ;
        if (rho_grad != nullptr) delete [] rho_grad;
        if (vp_pow   != nullptr) delete [] vp_pow  ;
        if (rho_pow  != nullptr) delete [] rho_pow ;
        if (vs       != nullptr) delete [] vs      ;
        if (vs_grad  != nullptr) delete [] vs_grad ;
        if (vs_pow   != nullptr) delete [] vs_pow  ;
    
        // for anisotropy (thomsen)
        if (epsilon      != nullptr) delete [] epsilon     ;
        if (epsilon_grad != nullptr) delete [] epsilon_grad;
        if (epsilon_pow  != nullptr) delete [] epsilon_pow ;
        if (delta        != nullptr) delete [] delta       ;
        if (delta_grad   != nullptr) delete [] delta_grad  ;
        if (delta_pow    != nullptr) delete [] delta_pow   ;
        if (dip          != nullptr) delete [] dip         ;
        if (dip_grad     != nullptr) delete [] dip_grad    ;
        if (dip_pow      != nullptr) delete [] dip_pow     ;
    
        // for vti_prem
        if (vph      != nullptr) delete [] vph     ;
        if (vpv      != nullptr) delete [] vpv     ;
        if (vsv      != nullptr) delete [] vsv     ;
        if (eta      != nullptr) delete [] eta     ;
        if (vph_grad != nullptr) delete [] vph_grad;
        if (vpv_grad != nullptr) delete [] vpv_grad;
        if (vsv_grad != nullptr) delete [] vsv_grad;
        if (eta_grad != nullptr) delete [] eta_grad;
        if (vph_pow  != nullptr) delete [] vph_pow ;
        if (vpv_pow  != nullptr) delete [] vpv_pow ;
        if (vsv_pow  != nullptr) delete [] vsv_pow ;
        if (eta_pow  != nullptr) delete [] eta_pow ;
    
        // for thomsen
        if (vp0      != nullptr) delete [] vp0     ;
        if (vs0      != nullptr) delete [] vs0     ;
        if (vp0_grad != nullptr) delete [] vp0_grad;
        if (vs0_grad != nullptr) delete [] vs0_grad;
        if (vp0_pow  != nullptr) delete [] vp0_pow ;
        if (vs0_pow  != nullptr) delete [] vs0_pow ;

        // for anisotropy c_ij, call one component
        if (c11 != nullptr) delete [] c11;
        if (c13 != nullptr) delete [] c13;
        if (c15 != nullptr) delete [] c15;
        if (c33 != nullptr) delete [] c33;
        if (c35 != nullptr) delete [] c35;
        if (c55 != nullptr) delete [] c55;
        if (c11_grad != nullptr) delete [] c11_grad;
        if (c13_grad != nullptr) delete [] c13_grad;
        if (c15_grad != nullptr) delete [] c15_grad;
        if (c33_grad != nullptr) delete [] c33_grad;
        if (c35_grad != nullptr) delete [] c35_grad;
        if (c55_grad != nullptr) delete [] c55_grad;
        if (c11_pow != nullptr) delete [] c11_pow;
        if (c13_pow != nullptr) delete [] c13_pow;
        if (c15_pow != nullptr) delete [] c15_pow;
        if (c33_pow != nullptr) delete [] c33_pow;
        if (c35_pow != nullptr) delete [] c35_pow;
        if (c55_pow != nullptr) delete [] c55_pow;
    }

};


bool isEqual(float a, float b);

void printProgress(float slowk);

void PrintIsPointOutOfInterfaceRange(Point2 A, 
    int ix, int iz, 
    float MINX, float MAXX);

/*------ find point and interpolation -------*/
int findLastGreaterEqualIndex(
    float value, 
    std::vector<float> &x);

int findFirstGreaterEqualIndex(
    float value, 
    std::vector<float> &x);

int findNearestNeighborIndex(
    float value, std::vector<float> &x);

float LinearInterpolation(
    std::vector<float> &x, 
    float *v,
    float xq);

float LinearInterp(
    int n,
    const float *x, 
    const float *v,
    float xq);

/*---- matrix: just for Bond transform ----*/
template <typename T>
class Matrix;

template <typename T>
std::ostream &operator<<(std::ostream& out, const Matrix<T> &mat);

template <typename T>
class Matrix{
private:
    int row;
    int col;
    T* p;
public:
    Matrix(int r, int c);
    Matrix(int r, int c, T* initval);
    Matrix(const Matrix<T> &mat); // copy
    ~Matrix();
    Matrix<T> operator*(const Matrix &mat);
    Matrix<T> operator/(T f);
    Matrix<T> operator+(const Matrix &mat);
    Matrix<T> operator-(const Matrix &mat);
    Matrix<T> operator=(const Matrix &mat);
    Matrix<T> transpose();
    Matrix<T> inverse2x2();
    T &operator()(int i, int j)const;
    friend std::ostream &operator<< <T>(std::ostream& out, const Matrix<T> &mat);
};

/*-----------------------------------------*/

int NumOfValues(std::vector<int> v, int NI);

void GenerateHalfGrid(
    size_t nx, 
    size_t nz,
    int grid_type, 
    const float *Gx,   // gridx
    const float *Gz,   // gridz
    float **Hx,        // half gridx    
    float **Hz);

/*
 *    ↑ +z  
 *    |     
 *         3----2 
 *         |    | 
 *         |    |
 *         0----1
 *
 */
Mesh2 GenerateHalfMesh(int grid_type,
                size_t ix, size_t iz, size_t indx, 
                size_t nx, size_t siz_slice,
                float *Hx, float *Hz);

/*
 *    ↑ +z  
 *    |     
 *         3----2 
 *         |    | 
 *         |    |
 *         0----1
 *
 */
Point2 *MeshSubdivide(Mesh2 M);

//======================for vti and tti ====================================
void para2vti(
    std::vector<float> &var, // input var
    int media_type, // return cij
    float &c11_2d,
    float &c33_2d,
    float &c55_2d,
    float &c13_2d,
    float &rho_2d); 

/* Reference: Bond, 1943, The Mathematics of the Physical Properties of Crystals */
/* theta: dip, phi: azimuth*/
void BondTransform2d(float c11, float c13, float c15, 
                     float c33, float c35, float c55,
                     float theta,
                     float &c11_tti, float &c13_tti, float &c15_tti, 
                     float &c33_tti, float &c35_tti, float &c55_tti);

void para2tti(std::vector<float> const &var, // input var
             int media_type, // return cij
             float &c11_2d,
             float &c13_2d,
             float &c15_2d,
             float &c33_2d,
             float &c35_2d,
             float &c55_2d,
             float &rho_2d);

#endif
