/*! \file dwt.c
    \brief Implementation of the discrete wavelet transform

    This file actually implements the wavelet transform. 
    It references rwt_platform.h and does not use any matlab-specific functions, 
    though it does use the mat() macro which describes the in-memory layout of 
    1D and 2D matrices in matlab
*/

#include "rwt_platform.h"

/*!
 * Perform convolution for dwt
 *
 * @param x_in input signal values
 * @param lx the length of x
 * @param h0 the low pass coefficients
 * @param h1 the high pass coefficients
 * @param lh_minus_one one less than the number of scaling coefficients
 * @param x_out_low low pass results
 * @param x_out_high high pass results
 * 
 * For the convolution we will calculate the output of the lowpass and highpass filters in parallel
 *
 * Normally we can describe the calculation of a convolution as
 * \f$ (\textbf{w} * \textbf{z})_k = \frac{1}{N} \sum\limits_{l=0}^{2N-1} w_{k-l} \cdot z_{l} \f$
 *
 */
void dwt_convolution(double *x_in, int lx, double *h0, double *h1, int lh_minus_one, double *x_out_low, double *x_out_high) {
  int i, j, ind;
  double x0, x1;
  for (i=lx; i<lx+lh_minus_one; i++) { 
    x_in[i] = *(x_in+(i-lx)); /*! extend x_in by creating a small mirror at the end of length lh_minus_one */
  }
  ind = 0;
  for (i=0; i<(lx); i+=2) {
    x0 = 0;
    x1 = 0;
    for (j=0; j<=lh_minus_one; j++) {
      x0 = x0 + x_in[i+j] * h0[lh_minus_one-j];
      x1 = x1 + x_in[i+j] * h1[lh_minus_one-j];
    }
    x_out_low[ind] = x0;
    x_out_high[ind++] = x1;
  }
}


/*!
 * Allocate memory for dwt
 *
 * @param m the number of rows of the input matrix
 * @param n the number of columns of the input matrix
 * @param lh the number of scaling coefficients
 * @param xdummy
 * @param y_dummy_low
 * @param y_dummy_high
 * @param h0
 * @param h1
 *
 * The low pass and high pass filter coefficients are the same size as the scaling coefficients
 * For the output storage area we will need as much space as the input: m*n
 * For the input storage area we will need the same plus one less than the length of the coeffiecients
 */
void dwt_allocate(int m, int n, int lh, double **xdummy, double **y_dummy_low, double **y_dummy_high, double **h0, double **h1) {
  *xdummy       = (double *) rwt_calloc(max(m,n)+lh-1, sizeof(double));
  *y_dummy_low  = (double *) rwt_calloc(max(m,n),      sizeof(double));
  *y_dummy_high = (double *) rwt_calloc(max(m,n),      sizeof(double));
  *h0           = (double *) rwt_calloc(lh,            sizeof(double));
  *h1           = (double *) rwt_calloc(lh,            sizeof(double));
}


/*!
 * Free memory that we allocated for dwt
 *
 * @param xdummy
 * @param y_dummy_low
 * @param y_dummy_high
 * @param h0
 * @param h1
 *
 */
void dwt_free(double **xdummy, double **y_dummy_low, double **y_dummy_high, double **h0, double **h1) {
  rwt_free(*xdummy);
  rwt_free(*y_dummy_low);
  rwt_free(*y_dummy_high);
  rwt_free(*h0);
  rwt_free(*h1);
}


/*!
 * Put the scaling coeffients into a form ready for use in the convolution function
 *
 * @param lh length of h / the number of scaling coefficients
 * @param h  the wavelet scaling coefficients
 * @param h0 the low pass coefficients - reversed h
 * @param h1 the high pass coefficients - forward h, even values are sign flipped
 *
 * The coefficients of our Quadrature Mirror Filter are described by
 * \f$ g\left[lh - 1 - n \right] = (-1)^n * h\left[n\right] \f$
 *
 */
void dwt_coefficients(int lh, double *h, double **h0, double **h1) {
  int i;
  for (i=0; i<lh; i++) {
    (*h0)[i] = h[(lh-i)-1];
    (*h1)[i] = h[i];
  }
  for (i=0; i<lh; i+=2)
    (*h1)[i] = -((*h1)[i]);
}


/*!
 * Perform the discrete wavelet transform
 *
 * @param x  the input signal
 * @param m  number of rows in the input
 * @param n  number of columns in the input
 * @param h  wavelet scaling coefficients
 * @param lh length of h / the number of scaling coefficients
 * @param L  the number of levels
 * @param y  the output signal with the wavelet transform applied
 *
 * The discrete wavelet transform begins with a set of samples of a signal whose length
 * is a power of 2. This exponent we shall call 'L' as it corresponds to the number of
 * levels in the filter bank for the calculation of the wavelet transform. 
 *
 * We shall use the name 'a' to refer to the approximation coefficients.
 * However, the actual implementation will not store this information separately, 
 * but rather will place it in the available space in the output array.
 *
 */
void dwt(double *x, int m, int n, double *h, int lh, int L, double *y) {
  double  *h0, *h1, *y_dummy_low, *y_dummy_high, *xdummy;
  long i;
  int current_level, lh_minus_one;
  int current_rows, current_cols, row_of_a, column_of_a, idx_rows, idx_columns;

  if (n==1) { /*! If passed a 1d column vector, just treat it as a row vector */
    n = m;
    m = 1;
  }
  
  dwt_allocate(m, n, lh, &xdummy, &y_dummy_low, &y_dummy_high, &h0, &h1);
  dwt_coefficients(lh, h, &h0, &h1); /*! For performance, calculate what we can outside the loops */
  lh_minus_one = lh - 1;
  current_rows = 2*m; /*! current_rows and current_cols start at 2x since we divide by 2 at the start of the loop */
  current_cols = 2*n;
 
  for (current_level=1; current_level<=L; current_level++) {
    if (m==1)
      current_rows = 1;
    else{
      current_rows = current_rows/2;
      row_of_a = current_rows/2;     
    }
    current_cols = current_cols/2;
    column_of_a = current_cols/2;

    for (idx_rows=0; idx_rows<current_rows; idx_rows++) {
      for (i=0; i<current_cols; i++)
	if (current_level==1)  
	  xdummy[i] = mat(x, idx_rows, i, m, n);  
	else 
	  xdummy[i] = mat(y, idx_rows, i, m, n);  
      /*! Perform filtering lowpass and highpass*/
      dwt_convolution(xdummy, current_cols, h0, h1, lh_minus_one, y_dummy_low, y_dummy_high); 
      /*! Restore dummy variables in matrices */
      idx_columns = column_of_a;
      for (i=0; i<column_of_a; i++) {    
	mat(y, idx_rows, i,             m, n) = y_dummy_low[i];  
	mat(y, idx_rows, idx_columns++, m, n) = y_dummy_high[i];  
      } 
    }  
    
    /*! For the 2d transform, we go through each of the columns after having gone through the rows */
    if (m>1) {
      for (idx_columns=0; idx_columns<current_cols; idx_columns++) { /* loop over columns */
	/*! Store in dummy variables */
	for (i=0; i<current_rows; i++)
	  xdummy[i] = mat(y, i, idx_columns, m, n);  
	/*! Perform filtering lowpass and highpass*/
	dwt_convolution(xdummy, current_rows, h0, h1, lh_minus_one, y_dummy_low, y_dummy_high); 
	/*! Restore dummy variables in matrix */
	idx_rows = row_of_a;
	for (i=0; i<row_of_a; i++) {
	  mat(y, i,          idx_columns, m, n) = y_dummy_low[i];  
	  mat(y, idx_rows++, idx_columns, m, n) = y_dummy_high[i];  
	}
      }
    }
  }
  dwt_free(&xdummy, &y_dummy_low, &y_dummy_high, &h0, &h1);
}

