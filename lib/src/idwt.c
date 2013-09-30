/*! \file idwt.c
    \brief Implementation of the inverse discrete wavelet transform

*/

#include "rwt_platform.h"

/*!
 * Perform convolution for idwt
 *
 * @param x_out
 * @param lx
 * @param g0
 * @param g1
 * @param lh_minus_one
 * @param lh_halved_minus_one
 * @param x_in_low
 * @param x_in_high
 * 
 */
void idwt_convolution(double *x_out, size_t lx, double *g0, double *g1, int lh_minus_one, int lh_halved_minus_one, double *x_in_low, double *x_in_high) {
  int k;
  size_t i, j, ind, tj;
  double x0, x1;

  for (k=lh_halved_minus_one-1; k > -1; k--) {
    x_in_low[k]  = x_in_low[lx+k];
    x_in_high[k] = x_in_high[lx+k];
  }

  ind = 0;
  for (i=0; i<(lx); i++) {
    x0 = 0;
    x1 = 0;
    tj = 0;
    for (j=0; j<=lh_halved_minus_one; j++) {
      x0 = x0 + x_in_low[i+j]*g0[lh_minus_one-1-tj] + x_in_high[i+j]*g1[lh_minus_one-1-tj];
      x1 = x1 + x_in_low[i+j]*g0[lh_minus_one-tj] + x_in_high[i+j]*g1[lh_minus_one-tj];
      tj += 2;
    }
    x_out[ind++] = x0;
    x_out[ind++] = x1;
  }
}


/*!
 * Allocate memory for idwt
 *
 * @param m the number of rows of the input matrix
 * @param n the number of columns of the input matrix
 * @param lh the number of scaling coefficients
 * @param xdummy
 * @param y_dummy_low
 * @param y_dummy_high
 * @param g0
 * @param g1
 *
 */
void idwt_allocate(size_t m, size_t n, int lh, double **xdummy, double **y_dummy_low, double **y_dummy_high, double **g0, double **g1) {
  *xdummy       = (double *) rwt_calloc(max(m,n),        sizeof(double));
  *y_dummy_low  = (double *) rwt_calloc(max(m,n)+lh/2-1, sizeof(double));
  *y_dummy_high = (double *) rwt_calloc(max(m,n)+lh/2-1, sizeof(double));
  *g0           = (double *) rwt_calloc(lh,              sizeof(double));
  *g1           = (double *) rwt_calloc(lh,              sizeof(double));
}


/*!
 * Free memory we allocated for idwt
 *
 * @param xdummy
 * @param y_dummy_low
 * @param y_dummy_high
 * @param g0
 * @param g1
 *
 */
void idwt_free(double **xdummy, double **y_dummy_low, double **y_dummy_high, double **g0, double **g1) {
  rwt_free(*xdummy);
  rwt_free(*y_dummy_low);
  rwt_free(*y_dummy_high);
  rwt_free(*g0);
  rwt_free(*g1);
}


/*!
 * Put the scaling coeffients into a form ready for use in the convolution function
 *
 * @param lh length of h / the number of scaling coefficients
 * @param h  the wavelet scaling coefficients
 * @param g0 same as h
 * @param g1 reversed h, even values are sign flipped
 *
 */
void idwt_coefficients(int lh, double *h, double **g0, double **g1) {
  int i;
  for (i=0; i<lh; i++){
    (*g0)[i] = h[i];
    (*g1)[i] = h[lh-i-1];
  }
  for (i=1; i<=lh; i+=2)
    (*g1)[i] = -((*g1)[i]);
}


/*!
 * Perform the inverse discrete wavelet transform
 *
 * @param x  the output signal with the inverse wavelet transform applied
 * @param m  number of rows in the input
 * @param n  number of columns in the input
 * @param h  wavelet scaling coefficients
 * @param lh length of h / the number of scaling coefficients
 * @param L  the number of levels
 * @param y  the input signal
 *
 */
void idwt(double *x, size_t m, size_t n, double *h, int lh, int L, double *y) {
  double  *g0, *g1, *y_dummy_low, *y_dummy_high, *xdummy;
  long i;
  int current_level, lh_minus_one, lh_halved_minus_one, sample_f;
  size_t current_rows, current_cols, row_cursor, column_cursor, idx_rows, idx_cols;

  idwt_allocate(m, n, lh, &xdummy, &y_dummy_low, &y_dummy_high, &g0, &g1);
  idwt_coefficients(lh, h, &g0, &g1);

  if (n==1) {
    n = m;
    m = 1;
  }
  
  lh_minus_one = lh - 1;
  lh_halved_minus_one = lh/2 - 1;
  /* 2^L */
  sample_f = 1;
  for (i=1; i<L; i++)
    sample_f = sample_f*2;
  
  if (m>1)
    current_rows = m/sample_f;
  else 
    current_rows = 1;
  current_cols = n/sample_f;

  for (i=0; i<(m*n); i++)
    x[i] = y[i];
  
  /* main loop */
  for (current_level=L; current_level >= 1; current_level--){
    row_cursor = current_rows/2;
    column_cursor = current_cols/2;
    
    /* go by columns in case of a 2D signal*/
    if (m>1){
      for (idx_cols=0; idx_cols<current_cols; idx_cols++){            /* loop over column */
	/* store in dummy variables */
	idx_rows = row_cursor;
	for (i=0; i<row_cursor; i++){    
	  y_dummy_low[i+lh_halved_minus_one]  = mat(x, i,          idx_cols, m, n);  
	  y_dummy_high[i+lh_halved_minus_one] = mat(x, idx_rows++, idx_cols, m, n);  
	}
	/* perform filtering lowpass and highpass*/
	idwt_convolution(xdummy, row_cursor, g0, g1, lh_minus_one, lh_halved_minus_one, y_dummy_low, y_dummy_high); 
	/* restore dummy variables in matrix */
	for (i=0; i<current_rows; i++)
	  mat(x, i, idx_cols, m, n) = xdummy[i];  
      }
    }
    /* go by rows */
    for (idx_rows=0; idx_rows<current_rows; idx_rows++){            /* loop over rows */
      /* store in dummy variable */
      idx_cols = column_cursor;
      for  (i=0; i<column_cursor; i++){    
	y_dummy_low[i+lh_halved_minus_one]  = mat(x, idx_rows, i,          m, n);  
	y_dummy_high[i+lh_halved_minus_one] = mat(x, idx_rows, idx_cols++, m, n);  
      } 
      /* perform filtering lowpass and highpass*/
      idwt_convolution(xdummy, column_cursor, g0, g1, lh_minus_one, lh_halved_minus_one, y_dummy_low, y_dummy_high); 
      /* restore dummy variables in matrices */
      for (i=0; i<current_cols; i++)
        mat(x, idx_rows, i, m, n) = xdummy[i];  
    }  
    if (m==1)
      current_rows = 1;
    else
      current_rows = current_rows*2;
    current_cols = current_cols*2;
  }
  idwt_free(&xdummy, &y_dummy_low, &y_dummy_high, &g0, &g1);
}

