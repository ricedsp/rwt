/*! \file rdwt.c
    \brief Implementation of the redundant discrete wavelet transform

*/

#include "rwt_platform.h"

/*!
 * Perform convolution for rdwt
 *
 * @param x_in input signal values
 * @param lx the length of x
 * @param h0 the low pass coefficients
 * @param h1 the high pass coefficients
 * @param lh the number of scaling coefficients
 * @param x_out_low low pass results
 * @param x_out_high high pass results
 * 
 */
void rdwt_convolution(double *x_in, int lx, double *h0, double *h1, int lh, double *x_out_low, double *x_out_high) {
  int i, j;
  double x0, x1;

  for (i=lx; i < lx+lh-1; i++)
    x_in[i] = x_in[i-lx];
  for (i=0; i<lx; i++){
    x0 = 0;
    x1 = 0;
    for (j=0; j<lh; j++){
      x0 = x0 + x_in[j+i]*h0[lh-1-j];
      x1 = x1 + x_in[j+i]*h1[lh-1-j];
    }
    x_out_low[i] = x0;
    x_out_high[i] = x1;
  }
}


/*!
 * Allocate memory for rdwt
 *
 * @param m the number of rows of the input matrix
 * @param n the number of columns of the input matrix
 * @param lh the number of scaling coefficients
 * @param x_dummy_low
 * @param x_dummy_high
 * @param y_dummy_low_low
 * @param y_dummy_low_high
 * @param y_dummy_high_low
 * @param y_dummy_high_high
 * @param h0
 * @param h1
 *
 */
void rdwt_allocate(int m, int n, int lh, double **x_dummy_low, double **x_dummy_high, double **y_dummy_low_low, 
  double **y_dummy_low_high, double **y_dummy_high_low, double **y_dummy_high_high, double **h0, double **h1) {
  *x_dummy_low       = (double *) rwt_calloc(max(m,n)+lh-1, sizeof(double));
  *x_dummy_high      = (double *) rwt_calloc(max(m,n)+lh-1, sizeof(double));
  *y_dummy_low_low   = (double *) rwt_calloc(max(m,n),      sizeof(double));
  *y_dummy_low_high  = (double *) rwt_calloc(max(m,n),      sizeof(double));
  *y_dummy_high_low  = (double *) rwt_calloc(max(m,n),      sizeof(double));
  *y_dummy_high_high = (double *) rwt_calloc(max(m,n),      sizeof(double));
  *h0                = (double *) rwt_calloc(lh,            sizeof(double));
  *h1                = (double *) rwt_calloc(lh,            sizeof(double));
}


/*!
 * Free memory that we allocated for dwt
 *
 * @param x_dummy_low
 * @param x_dummy_high
 * @param y_dummy_low_low
 * @param y_dummy_low_high
 * @param y_dummy_high_low
 * @param y_dummy_high_high
 * @param h0
 * @param h1
 *
 */
void rdwt_free(double **x_dummy_low, double **x_dummy_high, double **y_dummy_low_low, double **y_dummy_low_high, 
  double **y_dummy_high_low, double **y_dummy_high_high, double **h0, double **h1) {
  rwt_free(*x_dummy_low);
  rwt_free(*x_dummy_high);
  rwt_free(*y_dummy_low_low);
  rwt_free(*y_dummy_low_high);
  rwt_free(*y_dummy_high_low);
  rwt_free(*y_dummy_high_high);
  rwt_free(*h0);
  rwt_free(*h1);
}


/*!
 * Put the scaling coeffients into a form ready for use in the convolution function
 *
 * @param lh length of h / the number of scaling coefficients
 * @param h  the wavelet scaling coefficients
 * @param h0 the high pass coefficients - reversed h
 * @param h1 the high pass coefficients - forward h, even values are sign flipped
 *
 * The coefficients of our Quadrature Mirror Filter are described by
 * \f$ g\left[lh - 1 - n \right] = (-1)^n * h\left[n\right] \f$
 *
 * This is identical to dwt_coefficients() 
 *
 */
void rdwt_coefficients(int lh, double *h, double **h0, double **h1) {
  int i;
  for (i=0; i<lh; i++) {
    (*h0)[i] = h[(lh-i)-1];
    (*h1)[i] = h[i];
  }
  for (i=0; i<lh; i+=2)
    (*h1)[i] = -((*h1)[i]);
}


/*!
 * Perform the redundant discrete wavelet transform
 *
 * @param x  the input signal
 * @param m  number of rows in the input
 * @param n  number of columns in the input
 * @param h  wavelet scaling coefficients
 * @param lh length of h / the number of scaling coefficients
 * @param L  the number of levels
 * @param yl
 * @param yh
 *
 */
void rdwt(double *x, int m, int n, double *h, int lh, int L, double *yl, double *yh) {
  double  *h0, *h1, *y_dummy_low_low, *y_dummy_low_high, *y_dummy_high_low;
  double *y_dummy_high_high, *x_dummy_low , *x_dummy_high;
  long i;
  int current_level, current_rows, current_cols, column_cursor;
  int idx_cols, n_c, n_cb, idx_rows, n_r, n_rb;
  int column_cursor_plus_n, column_cursor_plus_double_n, sample_f;
  int three_n_L;

  rdwt_allocate(m, n, lh, &x_dummy_low, &x_dummy_high, &y_dummy_low_low, &y_dummy_low_high, 
    &y_dummy_high_low, &y_dummy_high_high, &h0, &h1);

  rdwt_coefficients(lh, h, &h0, &h1);

  if (n==1) {
    n = m;
    m = 1;
  }  

  /* analysis lowpass and highpass */
 
  three_n_L = 3*n*L;
  current_rows = 2*m;
  current_cols = 2*n;
  for (i=0; i<m*n; i++)
    yl[i] = x[i];
  
  /* main loop */
  sample_f = 1;
  for (current_level=1; current_level <= L; current_level++) {
    current_rows = current_rows/2;
    current_cols = current_cols/2;
    /* actual (level dependent) column offset */
    if (m==1)
      column_cursor = n*(current_level-1);
    else
      column_cursor = 3*n*(current_level-1);
    column_cursor_plus_n = column_cursor + n;
    column_cursor_plus_double_n = column_cursor_plus_n + n;
    
    /* go by rows */
    n_cb = n/current_cols;                 /* # of column blocks per row */
    for (idx_rows=0; idx_rows<m; idx_rows++) {           /* loop over rows */
      for (n_c=0; n_c<n_cb; n_c++) {   /* loop within one row */      
	/* store in dummy variable */
	idx_cols = -sample_f + n_c;
	for (i=0; i<current_cols; i++) {
	  idx_cols = idx_cols + sample_f;
	  x_dummy_low[i] = mat(yl, idx_rows, idx_cols, m, n);  
	}
	/* perform filtering lowpass/highpass */
	rdwt_convolution(x_dummy_low, current_cols, h0, h1, lh, y_dummy_low_low, y_dummy_high_high); 
	/* restore dummy variables in matridx_colses */
	idx_cols = -sample_f + n_c;
	for  (i=0; i<current_cols; i++) {
          idx_cols = idx_cols + sample_f;
          mat(yl, idx_rows, idx_cols,                 m, n)         = y_dummy_low_low[i];
          mat(yh, idx_rows, idx_cols + column_cursor, m, three_n_L) = y_dummy_high_high[i];  
	} 
      }
    }
      
    /* go by columns in case of a 2D signal*/
    if (m>1) {
      n_rb = m/current_rows;                 /* # of row blocks per column */
      for (idx_cols=0; idx_cols<n; idx_cols++) {           /* loop over column */
	for (n_r=0; n_r<n_rb; n_r++) {   /* loop within one column */
	  /* store in dummy variables */
	  idx_rows = -sample_f + n_r;
	  for (i=0; i<current_rows; i++) {    
	    idx_rows = idx_rows + sample_f;
	    x_dummy_low[i]  = mat(yl, idx_rows, idx_cols,               m, n);
	    x_dummy_high[i] = mat(yh, idx_rows, idx_cols + column_cursor, m, three_n_L);
	  }
	  /* perform filtering: first LL/LH, then HL/HH */
	  rdwt_convolution(x_dummy_low,  current_rows, h0, h1, lh, y_dummy_low_low,  y_dummy_low_high);
	  rdwt_convolution(x_dummy_high, current_rows, h0, h1, lh, y_dummy_high_low, y_dummy_high_high);
	  /* restore dummy variables in matrices */
	  idx_rows = -sample_f + n_r;
	  for (i=0; i<current_rows; i++) {
	    idx_rows = idx_rows + sample_f;
	    mat(yl, idx_rows, idx_cols,                               m, n)         = y_dummy_low_low[i];
	    mat(yh, idx_rows, idx_cols + column_cursor,               m, three_n_L) = y_dummy_low_high[i];
	    mat(yh, idx_rows, idx_cols + column_cursor_plus_n,        m, three_n_L) = y_dummy_high_low[i];
	    mat(yh, idx_rows, idx_cols + column_cursor_plus_double_n, m, three_n_L) = y_dummy_high_high[i];
	  }
	}
      }
    }
    sample_f = sample_f*2;
  }
  rdwt_free(&x_dummy_low, &x_dummy_high, &y_dummy_low_low, &y_dummy_low_high, &y_dummy_high_low, &y_dummy_high_high, &h0, &h1);
}
