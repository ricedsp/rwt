/*! \file rdwt.c
    \brief Implementation of the redundant discrete wavelet transform

*/

#include "rwt_platform.h"

void fpconv(double *x_in, int lx, double *h0, double *h1, int lh, double *x_outl, double *x_outh) {
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
    x_outl[i] = x0;
    x_outh[i] = x1;
  }
}


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

/* This is identical to dwt_coefficients */
void rdwt_coefficients(int lh, double *h, double **h0, double **h1) {
  int i;
  for (i=0; i<lh; i++) {
    (*h0)[i] = h[(lh-i)-1];
    (*h1)[i] = h[i];
  }
  for (i=0; i<lh; i+=2)
    (*h1)[i] = -((*h1)[i]);
}


void rdwt(double *x, int m, int n, double *h, int lh, int L, double *yl, double *yh) {
  double  *h0, *h1, *y_dummy_low_low, *y_dummy_low_high, *y_dummy_high_low;
  double *y_dummy_high_high, *x_dummy_low , *x_dummy_high;
  long i;
  int actual_L, actual_m, actual_n, column_of_a, ir, n_c, n_cb;
  int ic, n_r, n_rb, column_of_a_plus_double_n, sample_f;

  rdwt_allocate(m, n, lh, &x_dummy_low, &x_dummy_high, &y_dummy_low_low, &y_dummy_low_high, 
    &y_dummy_high_low, &y_dummy_high_high, &h0, &h1);

  rdwt_coefficients(lh, h, &h0, &h1);

  if (n==1){
    n = m;
    m = 1;
  }  
  /* analysis lowpass and highpass */
  
  actual_m = 2*m;
  actual_n = 2*n;
  for (i=0; i<m*n; i++)
    yl[i] = x[i];
  
  /* main loop */
  sample_f = 1;
  for (actual_L=1; actual_L <= L; actual_L++){
    actual_m = actual_m/2;
    actual_n = actual_n/2;
    /* actual (level dependent) column offset */
    if (m==1)
      column_of_a = n*(actual_L-1);
    else
      column_of_a = 3*n*(actual_L-1);
    column_of_a_plus_double_n = column_of_a + 2*n;
    
    /* go by rows */
    n_cb = n/actual_n;                 /* # of column blocks per row */
    for (ir=0; ir<m; ir++){            /* loop over rows */
      for (n_c=0; n_c<n_cb; n_c++){    /* loop within one row */      
	/* store in dummy variable */
	ic = -sample_f + n_c;
	for (i=0; i<actual_n; i++){    
	  ic = ic + sample_f;
	  x_dummy_low[i] = mat(yl, ir, ic, m);  
	}
	/* perform filtering lowpass/highpass */
	fpconv(x_dummy_low, actual_n, h0, h1, lh, y_dummy_low_low, y_dummy_high_high); 
	/* restore dummy variables in matrices */
	ic = -sample_f + n_c;
	for  (i=0; i<actual_n; i++){    
	  ic = ic + sample_f;
	  mat(yl, ir, ic, m) = y_dummy_low_low[i];  
	  mat(yh, ir, column_of_a+ic, m) = y_dummy_high_high[i];  
	} 
      }
    }
      
    /* go by columns in case of a 2D signal*/
    if (m>1){
      n_rb = m/actual_m;                 /* # of row blocks per column */
      for (ic=0; ic<n; ic++){            /* loop over column */
	for (n_r=0; n_r<n_rb; n_r++){    /* loop within one column */
	  /* store in dummy variables */
	  ir = -sample_f + n_r;
	  for (i=0; i<actual_m; i++){    
	    ir = ir + sample_f;
	    x_dummy_low[i] = mat(yl, ir, ic, m);  
	    x_dummy_high[i] = mat(yh, ir,column_of_a+ic, m);  
	  }
	  /* perform filtering: first LL/LH, then HL/HH */
	  fpconv(x_dummy_low, actual_m, h0, h1, lh, y_dummy_low_low, y_dummy_low_high); 
	  fpconv(x_dummy_high, actual_m, h0, h1, lh, y_dummy_high_low, y_dummy_high_high); 
	  /* restore dummy variables in matrices */
	  ir = -sample_f + n_r;
	  for (i=0; i<actual_m; i++){    
	    ir = ir + sample_f;
	    mat(yl, ir, ic,                           m) = y_dummy_low_low[i];  
	    mat(yh, ir, column_of_a+ic,               m) = y_dummy_low_high[i];  
	    mat(yh, ir, column_of_a+n+ic,             m) = y_dummy_high_low[i];  
	    mat(yh, ir, column_of_a_plus_double_n+ic, m) = y_dummy_high_high[i];  
	  }
	}
      }
    }
    sample_f = sample_f*2;
  }
  rdwt_free(&x_dummy_low, &x_dummy_high, &y_dummy_low_low, &y_dummy_low_high, &y_dummy_high_low, &y_dummy_high_high, &h0, &h1);
}


