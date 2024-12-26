#include <iostream>
#include <vector>
#include <cfloat>
#include <cmath>
#include <string.h>
//#include <Eigen>
#include "media_utility.hpp"


// for report error: int -> string
std::map<int, std::string> create_md2str_map() 
{
    std::map<int, std::string> m;
    m[ONE_COMPONENT] = "one_componet";
    m[ACOUSTIC_ISOTROPIC] = "acoustic_isotropic";
    m[ELASTIC_ISOTROPIC] = "elastic_isotropic"; 
    m[ELASTIC_VTI_PREM] = "elastic_vti_prem";
    m[ELASTIC_VTI_THOMSEN] = "elastic_vti_thomsen"; 
    m[ELASTIC_VTI_CIJ] = "elastic_vti_cij";
    m[ELASTIC_TTI_THOMSEN] = "elastic_tti_thomsen"; 
    m[ELASTIC_TTI_BOND] = "elastic_tti_bond";
    m[ELASTIC_ANISO_CIJ] = "elastic_aniso_cij";
    return m;
}

bool isEqual(float a, float b) {
    return abs(a-b) < FLT_EPSILON;
}

void printProgress(float slowk) {

    if (slowk > 1) slowk = 1;
    int p = slowk * 50;

    std::cout << "\33[1A";
    std::cout << "  [" + std::string(p, '=') + ">" + std::string(50-p, ' ') << "]" << std::endl;

    fflush(stdout);
}

void PrintIsPointOutOfInterfaceRange(Point2 A, 
    int ix, int iz, 
    float MINX, float MAXX) 
{
    if (A.x < MINX || A.x > MAXX) {
        fprintf(stderr,"Error: Grid(%d, %d) = (%f, %f) is out of the INTERFACES MESH (x in [%f %f]),\n"\
                       "       please check the interfaces file!\n", ix, iz, A.x, A.z, MINX, MAXX);
        fflush(stderr);
        exit(1);
    }
}


/* 
 * Find the last index of the x-vector greater than or equal to value, (x[i] >= value)
 *  Just for the ordered array x is from largest to smallest.
 *  (used to find the nearest elevation, from top to bottom)
 *  If value > x[all], return -1.
 */
int findLastGreaterEqualIndex(
    float value, 
    std::vector<float> &x)
{

    if (value > x[0]) return -1;

    float dist = FLT_MAX, newDist;
    int indx = -1, n = x.size();

    for (int i = 0; i < n; i++) {
        if (x[i] == -FLT_MAX) continue;

        newDist = x[i]-value;
        if (newDist < 0) break; 
        if (newDist >= 0 && newDist <= dist) {
            dist = newDist;
            indx = i;
        }
    }

    return indx;
}


 /* 
 * Find the first index of the x-vector greater than or equal to value, (x[i] >= value)
 *  Just for the ordered array x is from largest to smallest.
 *  Once there is x[i] = value, return.
 *  If value > x[all], return -1.
 */
int findFirstGreaterEqualIndex(
    float value, 
    std::vector<float> &x)
{
    float dist = FLT_MAX, newDist;
    int idx = -1, n = x.size();

    for(size_t i = 0; i < n; i++) {
        if (x[i] == -FLT_MAX) continue;
         
        newDist = x[i] - value;   
        if (newDist >= 0 && newDist <= dist) {
            // If the value is equal to the value in x, return the first index of the value,
            if (newDist <= 1e-6) 
                return i;
            dist = newDist;
            idx = i;
        }
    }

    return idx;
}


/* 
 * Find the index of the nearest value (x[i] <= value), 
 *  the x vector can be arbitrary
 */
int findNearestNeighborIndex(
    float value, 
    std::vector<float> &x)
{
    float dist, newDist;
    int idx = -1, nx = x.size();
    // expolation: INF(idx = -1) or

    dist = FLT_MAX;
    for(int i = 0; i < nx; i++) {
        if (x[i] == -FLT_MAX) continue;
        newDist = value - x[i];
        if (newDist >= 0 && newDist < dist) {
            dist = newDist;
            idx = i;
        }
    }

    return idx;
}

/* 
 * LinearInterpolate to the (xq, vq) point 
 * It's just for the problem: the x vector is increment 
 */
float LinearInterpolation(
    std::vector<float> &x, 
    float *v,
    float xq)
{    

    float vq;
    int n = x.size();

    int i = std::lower_bound(x.begin(), x.end(), xq) - x.begin();
    if (i == n) {
       return v[n-1]; 
    } else if(i == 0) {
        return v[i];
    } else {
        if ( isEqual(xq, x[i]) )
            return v[i];
        else {
            vq = v[i-1] + (v[i]-v[i-1])/(x[i]-x[i-1])*(xq-x[i-1]);
        }
    }

    return vq;
}


/* 
 * LinearInterpolate to the (xq, vq) point 
 * It's just for the problem: the x vector is increment 
 * for layer2model, the xloc is a dynamic array
 */
float LinearInterp(
    int n,
    const float *x, 
    const float *v,
    float xq)
{    

    float vq;

    int i = std::lower_bound(x, x+n, xq) - x;
    if (i == n) {
       return v[n-1]; 
    } else if(i == 0) {
        return v[i];
    } else {
        if ( isEqual(xq, x[i]) )
            return v[i];
        else {
            vq = v[i-1] + (v[i]-v[i-1])/(x[i]-x[i-1])*(xq-x[i-1]);
        }
    }

    return vq;
}

//============== for matrix: just used in bond transform ===============
// construct: init
template<typename T>
Matrix<T>::Matrix(int r, int c):row(r), col(c){   
    p = new T[row*col];
}

// construct: init
template<typename T>
Matrix<T>::Matrix(int r, int c, T* initval){
    row = r;
    col = c;
    int sz = row*col;
    p = new T[sz];
    for (int i = 0; i < sz; i++) 
        p[i] = initval[i];
}

// construct: copy
template<typename T>
Matrix<T>::Matrix(const Matrix<T> &mat){
    row = mat.row;
    col = mat.col;
    int sz = row*col;
    p = new T[sz];
    for (int i = 0; i < sz; i++) {
        p[i] = mat.p[i];
    }
}

// destruct
template<typename T>
Matrix<T>::~Matrix(){
    if (p != nullptr)
        delete[] p;
}

// get value
template<typename T>
T& Matrix<T>::operator()(int i, int j) const {
    if (i < 0 || i > row-1 || i < 0 || j > col-1) {
        fprintf(stderr,"i=%d, j=%d out of bound!\n", i, j);
        fflush(stderr);
        exit(1);
    }
    return p[i*col+j];
}

// operator: equal
template<typename T>
Matrix<T> Matrix<T>::operator=(const Matrix<T> &mat){
    if (row != mat.row || col != mat.col) {
        fprintf(stderr,"The row and column do not match!\n");
        fflush(stderr);
        exit(1);
    } 
    for (int i = 0; i < row*col; i++)
        p[i] = mat.p[i];
    return *this;
}

// multiply
template<typename T>
Matrix<T> Matrix<T>::operator*(const Matrix<T> &mat){
    Matrix<T> ans(row, mat.col);
    if (col != mat.row) {
        fprintf(stderr,"Can not multiply!");
        fflush(stderr);
        exit(1);
    } else {
        for (int i = 0; i < row; i++) {
            for (int j = 0; j < mat.col; j++) {
                ans.p[i*mat.col + j] = 0;
                for (int k = 0; k < col; k++) {
                    ans.p[i*mat.col + j] += p[i*col + k]*mat.p[k*mat.col + j];
                }
            }
        }
    }
    return ans;
}

// divided
template<typename T>
Matrix<T> Matrix<T>::operator/(T f){
  Matrix<T> ans(row, col);
  
  for (int i = 0; i < row; i++) {
    for (int j = 0; j < col; j++) {
      ans.p[i*col + j] = p[i*col+j]/f;
    }
  }
  
  return ans;
}

// addition
template<typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T> &mat){
    Matrix<T> ans(row, col);
    if (col != mat.col || row != mat.row) {
        fprintf(stderr,"ERROR: Matrix + Two matric should be the same size!");
        fflush(stderr);
        exit(1);
    } else {
      for (int i = 0; i < row; i++) {
        for (int j = 0; j < mat.col; j++) {
          ans.p[i*col + j] = p[i*col + j] + mat.p[i*col + j];
        }
      }
    }
    return ans;
}

// minus
template<typename T>
Matrix<T> Matrix<T>::operator-(const Matrix<T> &mat){
    Matrix<T> ans(row, col);
    if (col != mat.col || row != mat.row) {
        fprintf(stderr,"ERROR: Matrix - Two matric should be the same size!");
        fflush(stderr);
        exit(1);
    } else {
      for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
          ans.p[i*col + j] = p[i*col + j] - mat.p[i*col + j];
        }
      }
    }
    return ans;
}

// tmp here, just for S-M method
template<typename T>
Matrix<T> Matrix<T>::inverse2x2() {
  if (row != 2 || col != 2) {
      fprintf(stderr,"The row and column can just be 2!\n");
      fflush(stderr);
      exit(1);
  } 
  Matrix<T> inversed(2,2, p);
  for (int k=0; k<2; k++)
  {
     float con = p[k*col+k];
     inversed.p[k*col+k] = 1.0;
  
     for (int i=0; i<2; i++) {
       inversed.p[k*col+i] = inversed.p[k*col+i]/con;
     }
  
     for (int i=0; i<2; i++)
     {
        if (i!=k) {
          con = inversed.p[i*col+k];
          inversed.p[i*col+k] = 0.0;
          for (int j=0; j<2; j++) {
            inversed.p[i*col+j] = inversed.p[i*col+j] - 
                                  inversed.p[k*col+j] * con;
          }
        }
     }
  }

  return inversed;
}

template<typename T>
Matrix<T> Matrix<T>::transpose() {
    Matrix<T> transposed(col, row);
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            transposed(j,i) = (*this)(i,j);
        }
    }
    return transposed;
}

// print
template<typename T>
std::ostream& operator<<(std::ostream &out, const Matrix<T> &mat) {
    for (int i = 0; i < mat.row; i++) {
        for (int j = 0; j < mat.col; j++) {
            out << mat(i,j);
            if (j == mat.col-1) 
                out << std::endl;
            else
                out << " ";
        }
    }
    return out;
}

// re-statement for the compiler determines the type
template class Matrix<float>;
template class Matrix<int>;

//====================================================================


//========== for equivalent medium parameterization method ===========

/* 
 * How many different values in the vector, 
 *  NI is the upper limitation of the v[i].
 */
int NumOfValues(std::vector<int> v, int NI) 
{
    std::vector<int> tab(NI, 0);
    int num = 0;
    for (size_t i = 0; i < v.size(); i++) {
        // If higher than interface, do not apply equivalent medium
        if (v[i] == -1)
            return 0;

        if (tab[v[i]] == 0)
            num++;
        tab[v[i]]++;
    }
    return num;
}

/*
 *    ↑ +z  
 *    |     
 *         3----2 
 *         |    | 
 *         |    |
 *         0----1
 *
 */
// one mesh is subdivided into ng+1 points in one direction
Point2 *MeshSubdivide(Mesh2 M) {

    int siz_line = NG+1;
    int siz_slice = (NG+1) * siz_line;

    Point2 *gd = new Point2[siz_slice];

    for (int k = 0; k <= NG; k++) {
        for (int i = 0; i <= NG; i++) {
            int indx = i + k * siz_line;
            // 1st x
            Point2 Lx0 = M.v[0] + (M.v[1]-M.v[0])/NG*i;
            Point2 Lx1 = M.v[3] + (M.v[2]-M.v[3])/NG*i;
            // interpolation in the z-axis direction.
            gd[indx] = Lx0 + (Lx1-Lx0)/NG*k;           
        }
    }    
    return gd;
}

void GenerateHalfGrid(
    size_t nx, 
    size_t nz,
    int grid_type, 
    const float *Gx,   // gridx
    const float *Gz,   // gridz
    float **Hx,        // half gridx    
    float **Hz)        // half gridz
{

    size_t siz_slice = nx*nz;
    /* 
     * half grid point is inside the grid, 
     * so the number of grid point in one direction is np-1.
     */
    if (grid_type == GRID_CART) {
        
        *Hx = new float [nx]; 
        *Hz = new float [nz]; 

        for (size_t i = 0; i < nx-1; i++)
            (*Hx)[i] = (Gx[i+1] + Gx[i])/2.0;
        (*Hx)[nx-1] = Gx[nx-1];

        for (size_t i = 0; i < nz-1; i++)
            (*Hz)[i] = (Gz[i+1] + Gz[i])/2.0;
        (*Hz)[nz-1] = Gz[nz-1];

    } else if (grid_type == GRID_VMAP) {

        *Hx = new float [nx]; 
        *Hz = new float [siz_slice]; 

        for (size_t i = 0; i < nx-1; i++)
            (*Hx)[i] = (Gx[i+1] + Gx[i])/2.0;
        (*Hx)[nx-1] = Gx[nx-1];

        for (size_t i = 0; i < siz_slice; ++i)
            (*Hz)[i] = Gz[i];

        for (size_t k = 0; k < nz-1; k++) {
            for (size_t i = 0; i < nx-1; i++) {
                size_t indx =  i + k*nx;
                (*Hz)[indx] = (Gz[indx] + Gz[indx+1] + 
                          Gz[indx+nx+1] + Gz[indx+nx] )/4.0;
            }
        }
    } else {

        *Hx = new float [siz_slice]; 
        *Hz = new float [siz_slice]; 

        for (size_t i = 0; i < siz_slice; i++) {
            (*Hx)[i] = Gx[i];
            (*Hz)[i] = Gz[i];
        }

        for (size_t k = 0; k < nz-1; k++) {
            for (size_t i = 0; i < nx-1; i++) {
                size_t indx =  i + k*nx;
    
                (*Hx)[indx] = (Gx[indx] + Gx[indx+1] + 
                          Gx[indx+nx+1] + Gx[indx+nx] )/4.0;
    
                (*Hz)[indx] = (Gz[indx] + Gz[indx+1] + 
                          Gz[indx+nx+1] + Gz[indx+nx] )/4.0;
                
            }
            
        }
    }

}

Mesh2 GenerateHalfMesh(int grid_type,
                size_t ix, size_t iz, size_t indx, 
                size_t nx, size_t siz_slice,
                float *Hx, float *Hz) 
{
    /* 
     * The H(ix, iz) = Bilinear of G(ix, iz) -> G(ix+1, iz+1),
     *  so, get the averaging of the point (ix, iz), 
     *  Mesh H(ix-1, iz-1) -> H(ix, iz) is needed.
     */
    Mesh2 M( Point2(Hx[ix-1], Hz[iz-1]),
             Point2(Hx[ix  ], Hz[iz-1]),
             Point2(Hx[ix  ], Hz[iz  ]),
             Point2(Hx[ix-1], Hz[iz  ]) );
    
    if (grid_type == GRID_VMAP) {
        size_t indx0 = indx-1-nx;
        size_t indx1 = indx  -nx;
        size_t indx2 = indx     ;
        size_t indx3 = indx-1   ;

        Mesh2 M( Point2(Hx[ix-1], Hz[indx0]),
                 Point2(Hx[ix  ], Hz[indx1]),
                 Point2(Hx[ix  ], Hz[indx2]),
                 Point2(Hx[ix-1], Hz[indx3]) );
        return M;

    } else if (grid_type == GRID_CURV) {
        size_t indx0 = indx-1-nx;
        size_t indx1 = indx  -nx;
        size_t indx2 = indx     ;
        size_t indx3 = indx-1   ;
        Mesh2 M( Point2( Hx[indx0], Hz[indx0] ),
                 Point2( Hx[indx1], Hz[indx1] ),
                 Point2( Hx[indx2], Hz[indx2] ),
                 Point2( Hx[indx3], Hz[indx3] ) );
        return M;
    }

    return M;
}



//======================for vti and tti ====================================
void para2vti(
    std::vector<float> &var, // input var
    int media_type, // return cij
    float &c11_2d,
    float &c33_2d,
    float &c55_2d,
    float &c13_2d,
    float &rho_2d) 
{
    /* Calculate cij */
    if (media_type == ELASTIC_VTI_PREM) {
    //-r eference: Dziewonski and Anderson, 1981, Preliminary reference Earth model
    /* rho, vph, vpv, vsh, vsv, ets */
        float vph = var[1], vpv = var[2];
        float vsv = var[3];
        float rho = var[0], eta = var[4];
    
        c11_2d = vph*vph*rho;
        c33_2d = vpv*vpv*rho;
        c55_2d = vsv*vsv*rho;
        rho_2d = rho; 
        c13_2d = eta*(c11_2d - 2.0*c55_2d);

    } else if (media_type == ELASTIC_VTI_THOMSEN) {
    //- reference: Thomsen, 1986, Weak elastic anisotropy 
    /* rho, vp0 (alpha), vs0 (beta), epsilon, delta, gamma */
        float rho = var[0];
        float vp0 = var[1], vs0 = var[2];
        float epsil = var[3], delta = var[4];
        float c33 = vp0*vp0*rho; 
        float c55 = vs0*vs0*rho;
        rho_2d = rho; 
        c33_2d = c33;
        c55_2d = c55;                
        c11_2d = 2.0*epsil * c33 + c33;
        c13_2d = sqrt( 2.0*delta*c33*(c33-c55) + (c33-c55)*(c33-c55) ) - c55;

    } else if (media_type == ELASTIC_VTI_CIJ) {
    /* rho, c11 c33 c55 c13*/
        c11_2d = var[1];
        c33_2d = var[2];
        c55_2d = var[3];
        c13_2d = var[4];
        rho_2d = var[0];

    } else {
        fprintf(stderr,"Error: Unknow VTI type, for code check, please contact Luqian Jiang!\n");
        fflush(stderr);
        exit(1);
    }
}

/* 
 * Reference: Bond, 1943, The Mathematics of the Physical Properties of Crystals 
 * Zhu and Dorman, 2000, Two-dimensional, three-component wave propagation in a transversely 
 *  isotropic medium with arbitrary-orientation–finite-element modeling
 * 
 * theta: dip
 */
void BondTransform2d(float c11, float c13, float c15, 
                     float c33, float c35, float c55,
                     float theta,
                     float &c11_tti, float &c13_tti, float &c15_tti, 
                     float &c33_tti, float &c35_tti, float &c55_tti) 
{
    float a11 = cos(theta),  a13 = -sin(theta);
    float a31 = sin(theta),  a33 =  cos(theta);   

    float R_tmp[9] = {
        a11*a11, a13*a13, 2*a13*a11,       
        a31*a31, a33*a33, 2*a33*a31,       
        a31*a11, a33*a13, a13*a31+a11*a33};

    float C0_tmp[9] = {
        c11, c13, c15,
        c13, c33, c35,
        c15, c35, c55 };

    
    Matrix<float> R(3,3,R_tmp);
    Matrix<float> C0(3,3,C0_tmp);
    Matrix<float> C(3,3);
    // C = R * C_vti * R^T
    C = R * C0 *R.transpose();

    c11_tti = C(0,0); c13_tti = C(0,1); c15_tti = C(0,2); 
    c33_tti = C(1,1); c35_tti = C(1,2);
    c55_tti = C(2,2); 
}

void para2tti(std::vector<float> const &var, // input var
             int media_type, // return cij
             float &c11_2d,
             float &c13_2d,
             float &c15_2d,
             float &c33_2d,
             float &c35_2d,
             float &c55_2d,
             float &rho_2d)

{
    if (media_type == ELASTIC_TTI_THOMSEN) {
    /* reference: Thomsen, 1986, Weak elastic anisotropy */    
        /* rho, vp0, vs0, epsilon, delta, gamma, dip */
        /* Calculate cij */
        float rho = var[0];
        float vp0 = var[1], vs0 = var[2];
        float epsil = var[3], delta = var[4];
        float dip = var[5];

        float c33 = vp0*vp0*rho;
        float c55 = vs0*vs0*rho;
        float c11 = 2.0*epsil*c33 + c33;
        float c13 = sqrt( 2.0*delta*c33*(c33-c55) + (c33-c55)*(c33-c55) ) - c55;
        rho_2d = rho; 
    
        BondTransform2d(c11, c13, 0, 
                        c33, 0, c55,                  
                        dip*(PI/180.0),
                        c11_2d, c13_2d, c15_2d,
                        c33_2d, c35_2d, c55_2d);

    } else if (media_type == ELASTIC_TTI_BOND) {
        /* rho c11 c33 c55 c13 dip */
        float c11 = var[1];
        float c33 = var[2];
        float c55 = var[3];
        float c13 = var[4];
        float dip = var[5];
        rho_2d = var[0];
        
        BondTransform2d(c11, c13, 0,  
                        c33, 0, c55, 
                        dip*(PI/180.0),
                        c11_2d, c13_2d, c15_2d,
                        c33_2d, c35_2d, c55_2d);

    } else if (media_type == ELASTIC_ANISO_CIJ) {
    /* rho c11 c13 c15 c33 c35 c55 */
    
        /* Calculate cij */
        rho_2d = var[0];
        c11_2d = var[1];
        c13_2d = var[2];
        c15_2d = var[3];
        c33_2d = var[4];
        c35_2d = var[5];
        c55_2d = var[6];

    } else {
        fprintf(stderr,"Error: Unknow TTI type, for code check, please contact Luqian Jiang!\n");
        fflush(stderr);
        exit(1);
    }
}



