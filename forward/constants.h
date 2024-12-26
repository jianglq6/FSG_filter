#ifndef CONSTANTS_H
#define CONSTANTS_H

// consts
#define CONST_NDIM     2
#define CONST_NDIM_2   4 // 2 * ndim
#define CONST_TIJ_SIZE 3 // 
#define CONST_MAX_STRLEN 1024

#ifndef M_PI
#define PI 3.14159265358979323846264338327950288419716939937510
#else
#define PI M_PI
#endif

// medium type
#define CONST_MEDIUM_ACOUSTIC  1
#define CONST_MEDIUM_ELASTIC   2

#define CONST_ANISO_ISO  1
#define CONST_ANISO_VTI  2
#define CONST_ANISO_TTI  3
#define CONST_ANISO_ALL  4

// visco type
#define CONST_VISCO_NONE    1
#define CONST_VISCO_GRAVES  2

#endif
