#include <stdio.h>
#include <stdlib.h>

#include <filter_fsg.h>

#define SIGMA 0.2


// -----------------------------------------------------------------------------//
//            vel is on integer grid     
// -----------------------------------------------------------------------------//
void
filter_stress_intv(
float *restrict Txx1,  float *restrict Tzz1,  float *restrict Txz1,
float *restrict Txx2,  float *restrict Tzz2,  float *restrict Txz2,
int lni1, int lni2, int lnj1, int lnj2, size_t siz_line_half,
int half_len, float *coef)
{
    for (size_t j=lnj1; j<=lnj2; j++)
    {
      size_t iptr_j = j * siz_line_half;

      size_t iptr = iptr_j + lni1;

      for (size_t i=lni1; i<=lni2; i++)      
	    {

		    // NO.2

		    filter_lower_right135(Txx1, Txx2, iptr, half_len, coef, siz_line_half, 1);
		    filter_lower_right45(Txx1, Txx2, iptr, half_len, coef, siz_line_half, 1);

		    filter_lower_right135(Tzz1, Tzz2, iptr, half_len, coef, siz_line_half, 1);
		    filter_lower_right45(Tzz1, Tzz2, iptr, half_len, coef, siz_line_half, 1);

		    filter_lower_right135(Txz2, Txz1, iptr, half_len, coef, siz_line_half, 1);
		    filter_lower_right45(Txz2, Txz1, iptr, half_len, coef, siz_line_half, 1);

		    // NO.3
		    filter_upper_left135(Txx2, Txx1, iptr, half_len, coef, siz_line_half, 1);
		    filter_upper_left45(Txx2, Txx1, iptr, half_len, coef, siz_line_half, 1);

		    filter_upper_left135(Tzz2, Tzz1, iptr, half_len, coef, siz_line_half, 1);
		    filter_upper_left45(Tzz2, Tzz1, iptr, half_len, coef, siz_line_half, 1);

		    filter_upper_left135(Txz1, Txz2, iptr, half_len, coef, siz_line_half, 1);
		    filter_upper_left45(Txz1, Txz2, iptr, half_len, coef, siz_line_half, 1);


		    iptr += 1;

	   }
	}
}


void 
filter_velocity_intv(
float *restrict Vx1,  float *restrict Vz1,
float *restrict Vx2,  float *restrict Vz2,
int lni1, int lni2, int lnj1, int lnj2, size_t siz_line_half,
int half_len, float *coef)
{

    for (size_t j=lnj1; j<=lnj2; j++)
    {
      size_t iptr_j = j * siz_line_half;
      size_t iptr = iptr_j + lni1;

      for (size_t i=lni1; i<=lni2; i++)      
	    {

		    // NO.1
		    filter_lower_left135(Vx1, Vx2, iptr, half_len,  coef, siz_line_half, 1);
		    filter_lower_left45(Vx1, Vx2, iptr, half_len,  coef, siz_line_half, 1);

		    filter_lower_left135(Vz2, Vz1, iptr, half_len,  coef, siz_line_half, 1);
		    filter_lower_left45(Vz2, Vz1, iptr, half_len,  coef, siz_line_half, 1);
		    

		    // NO.4
		    filter_upper_right135(Vz1, Vz2, iptr, half_len, coef, siz_line_half, 1);
		    filter_upper_right45(Vz1, Vz2, iptr, half_len, coef, siz_line_half, 1);

		    filter_upper_right135(Vx2, Vx1, iptr, half_len, coef, siz_line_half, 1);
		    filter_upper_right45(Vx2, Vx1, iptr, half_len, coef, siz_line_half, 1);


		    iptr += 1;
	    }
	}
}


// left_lower 
int
filter_lower_left135( float *u0, float *u1, size_t iptr, int half_len, float *coef,
		size_t dim2, size_t dim1)
{
  float Du;
  int indx;
  Du = coef[0] * u0[iptr];
  for (int i=1; i<= half_len; i = i+2) {
  	indx = (i+1)/2;
  
  	Du += coef[i] * ( u1[iptr + (indx-1) * dim2 - indx     * dim1] 
			              + u1[iptr - indx     * dim2 + (indx-1) * dim1]);
                                           

  	if ( i < half_len ) {
  
  	    Du += coef[i+1] * ( u0[iptr - indx * dim2 + indx * dim1] 
				                  + u0[iptr + indx * dim2 - indx * dim1]);
  	}
  
  }
  
  u0[iptr] -= SIGMA * Du;

  return 0;
}
  
int
filter_lower_left45( float *u0, float *u1, size_t iptr, int half_len, float *coef,
		size_t dim2, size_t dim1)
{
  float Du;
  int indx;
  Du = coef[0] * u0[iptr];
  for (int i=1; i<= half_len; i = i+2) {
  	indx = (i+1)/2;
  
  	Du += coef[i] * ( u1[iptr - indx     * dim2 - indx     * dim1] 
					          + u1[iptr + (indx-1) * dim2 + (indx-1) * dim1] );
                                           

  	if ( i < half_len ) {
  
  	    Du += coef[i+1] * ( u0[iptr - indx * dim2 - indx * dim1]
						              + u0[iptr + indx * dim2 + indx * dim1] );
  	}
  }

  u0[iptr] -= SIGMA * Du;

  return 0;
}


// left_upper
int
filter_upper_left135( float *u0, float *u1, size_t iptr, int half_len, float *coef,
		size_t dim2, size_t dim1)
{
  float Du;
  int indx;
  
  Du = coef[0] * u0[iptr];
  for (int i=1; i<= half_len; i = i+2) {
  	indx = (i+1)/2;
	             
	// First two items is 135 degree; The last two items is 45 degree
  	Du += coef[i] * ( u1[iptr - (indx-1) * dim2 + (indx-1)* dim1] 
			              + u1[iptr + indx     * dim2 - indx    * dim1]);

  	if ( i < half_len ) {
  
	// First two items is 135 degree; The last two items is 45 degree
  	    Du += coef[i+1] * ( u0[iptr - indx * dim2 + indx * dim1] 
				                  + u0[iptr + indx * dim2 - indx * dim1]);
  	}
  }

  u0[iptr] -= SIGMA * Du;

  return 0;
}

int
filter_upper_left45( float *u0, float *u1, size_t iptr, int half_len, float *coef,
		size_t dim2, size_t dim1)
{
  float Du;
  int indx;
  
  Du = coef[0] * u0[iptr];
  for (int i=1; i<= half_len; i = i+2) {
  	indx = (i+1)/2;
	             
	// First two items is 135 degree; The last two items is 45 degree
  	Du += coef[i] * ( u1[iptr - (indx-1) * dim2 - indx    * dim1] 
				            + u1[iptr + indx     * dim2 + (indx-1)* dim1]);

  	if ( i < half_len ) {
  
	// First two items is 135 degree; The last two items is 45 degree
  	    Du += coef[i+1] * ( u0[iptr - indx * dim2 - indx * dim1] 
						              + u0[iptr + indx * dim2 + indx * dim1]);
  	}
  }

  u0[iptr] -= SIGMA * Du;

  return 0;
}



// right_upper
int
filter_upper_right135( float *u0, float *u1, size_t iptr, int half_len, float *coef,
		size_t dim2, size_t dim1)
{
  float Du;
  int indx;
  Du = coef[0] * u0[iptr];
  for (int i=1; i<= half_len; i = i+2) {
  	indx = (i+1)/2;
  
  	Du += coef[i] * ( u1[iptr - (indx-1) * dim2 + indx     * dim1] 
			              + u1[iptr + indx     * dim2 - (indx-1) * dim1]);
  
  	if ( i < half_len ) {
  
  	    Du += coef[i+1] * ( u0[iptr - indx * dim2 + indx * dim1] 
				                  + u0[iptr + indx * dim2 - indx * dim1]);
  	}
  
  }

  u0[iptr] -= SIGMA * Du;
  
  return 0;
}

int
filter_upper_right45( float *u0, float *u1, size_t iptr, int half_len, float *coef,
		size_t dim2, size_t dim1)
{
  float Du;
  int indx;
  Du = coef[0] * u0[iptr];
  for (int i=1; i<= half_len; i = i+2) {
  	indx = (i+1)/2;
  
  	Du += coef[i] * ( u1[iptr - (indx-1) * dim2 - (indx-1) * dim1] 
					+ u1[iptr + indx     * dim2 + indx     * dim1] );
  
  	if ( i < half_len ) {
  
  	    Du += coef[i+1] * ( u0[iptr - indx * dim2 - indx * dim1] 
						              + u0[iptr + indx * dim2 + indx * dim1] );
  	}
  
  }

  u0[iptr] -= SIGMA * Du;
  
  return 0;
}


// right_lower
int
filter_lower_right135( float *u0, float *u1, size_t iptr, int half_len, float *coef,
		size_t dim2, size_t dim1)
{
  float Du;
  int indx;
  Du = coef[0] * u0[iptr];
  for (int i=1; i<= half_len; i = i+2) {
  	indx = (i+1)/2;
  
  	Du += coef[i] * ( u1[iptr - indx     * dim2 + indx    * dim1] 
			              + u1[iptr + (indx-1) * dim2 - (indx-1)* dim1]);
	                
  
  	if ( i < half_len ) {
  
  	    Du += coef[i+1] * ( u0[iptr - indx * dim2 + indx * dim1] 
				                  + u0[iptr + indx * dim2 - indx * dim1]);
  	}
  
  }

  u0[iptr] -= SIGMA * Du;

  return 0;
}	


int
filter_lower_right45( float *u0, float *u1, size_t iptr, int half_len, float *coef,
		size_t dim2, size_t dim1)
{
  float Du;
  int indx;
  Du = coef[0] * u0[iptr];
  for (int i=1; i<= half_len; i = i+2) {
  	indx = (i+1)/2;
  
  	Du += coef[i] * ( u1[iptr - indx     * dim2 - (indx-1)* dim1] 
          					+ u1[iptr + (indx-1) * dim2 + indx    * dim1]);
	                
  
  	if ( i < half_len ) {
  
  	    Du += coef[i+1] * ( u0[iptr - indx * dim2 - indx * dim1] 
					            	  + u0[iptr + indx * dim2 + indx * dim1]);
  	}
  
  }

  u0[iptr] -= SIGMA * Du;

  return 0;
}	







