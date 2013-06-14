/*! \file idwt.c
    \brief Implementation of the inverse discrete wavelet transform

*/

#include "rwt_platform.h"

void bpsconv(double *x_out, int lx, double *g0, double *g1, int lh_minus_one, int lh_halved_minus_one, double *x_inl, double *x_inh) {
  int i, j, ind, tj;
  double x0, x1;

  for (i=lh_halved_minus_one-1; i > -1; i--){
    x_inl[i] = x_inl[lx+i];
    x_inh[i] = x_inh[lx+i];
  }
  ind = 0;
  for (i=0; i<(lx); i++){
    x0 = 0;
    x1 = 0;
    tj = -2;
    for (j=0; j<=lh_halved_minus_one; j++){
      tj+=2;
      x0 = x0 + x_inl[i+j]*g0[lh_minus_one-1-tj] + x_inh[i+j]*g1[lh_minus_one-1-tj] ;
      x1 = x1 + x_inl[i+j]*g0[lh_minus_one-tj] + x_inh[i+j]*g1[lh_minus_one-tj] ;
    }
    x_out[ind++] = x0;
    x_out[ind++] = x1;
  }
}


void idwt_allocate(int m, int n, int lh, double **xdummy, double **y_dummy_low, double **y_dummy_high, double **g0, double **g1) {
  *xdummy       = (double *) rwt_calloc(max(m,n),        sizeof(double));
  *y_dummy_low  = (double *) rwt_calloc(max(m,n)+lh/2-1, sizeof(double));
  *y_dummy_high = (double *) rwt_calloc(max(m,n)+lh/2-1, sizeof(double));
  *g0           = (double *) rwt_calloc(lh,              sizeof(double));
  *g1           = (double *) rwt_calloc(lh,              sizeof(double));
}


void idwt_free(double **xdummy, double **y_dummy_low, double **y_dummy_high, double **g0, double **g1) {
  rwt_free(*xdummy);
  rwt_free(*y_dummy_low);
  rwt_free(*y_dummy_high);
  rwt_free(*g0);
  rwt_free(*g1);
}


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
void idwt(double *x, int m, int n, double *h, int lh, int L, double *y) {
  double  *g0, *g1, *y_dummy_low, *y_dummy_high, *xdummy;
  long i;
  int actual_L, actual_m, actual_n, row_of_a, column_of_a, ir, ic, lh_minus_one, lh_halved_minus_one, sample_f;

  idwt_allocate(m, n, lh, &xdummy, &y_dummy_low, &y_dummy_high, &g0, &g1);
  idwt_coefficients(lh, h, &g0, &g1);

  if (n==1){
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
    actual_m = m/sample_f;
  else 
    actual_m = 1;
  actual_n = n/sample_f;

  for (i=0; i<(m*n); i++)
    x[i] = y[i];
  
  /* main loop */
  for (actual_L=L; actual_L >= 1; actual_L--){
    row_of_a = actual_m/2;
    column_of_a = actual_n/2;
    
    /* go by columns in case of a 2D signal*/
    if (m>1){
      for (ic=0; ic<actual_n; ic++){            /* loop over column */
	/* store in dummy variables */
	ir = row_of_a;
	for (i=0; i<row_of_a; i++){    
	  y_dummy_low[i+lh_halved_minus_one]  = mat(x, i,    ic, m);  
	  y_dummy_high[i+lh_halved_minus_one] = mat(x, ir++, ic, m);  
	}
	/* perform filtering lowpass and highpass*/
	bpsconv(xdummy, row_of_a, g0, g1, lh_minus_one, lh_halved_minus_one, y_dummy_low, y_dummy_high); 
	/* restore dummy variables in matrix */
	for (i=0; i<actual_m; i++)
	  mat(x, i, ic, m) = xdummy[i];  
      }
    }
    /* go by rows */
    for (ir=0; ir<actual_m; ir++){            /* loop over rows */
      /* store in dummy variable */
      ic = column_of_a;
      for  (i=0; i<column_of_a; i++){    
	y_dummy_low[i+lh_halved_minus_one]  = mat(x, ir, i,    m);  
	y_dummy_high[i+lh_halved_minus_one] = mat(x, ir, ic++, m);  
      } 
      /* perform filtering lowpass and highpass*/
      bpsconv(xdummy, column_of_a, g0, g1, lh_minus_one, lh_halved_minus_one, y_dummy_low, y_dummy_high); 
      /* restore dummy variables in matrices */
      for (i=0; i<actual_n; i++)
        mat(x, ir, i, m) = xdummy[i];  
    }  
    if (m==1)
      actual_m = 1;
    else
      actual_m = actual_m*2;
    actual_n = actual_n*2;
  }
  idwt_free(&xdummy, &y_dummy_low, &y_dummy_high, &g0, &g1);
}

