/*! \file irdwt.c
    \brief Implementation of the inverse redundant discrete wavelet transform

*/

#include "rwt_platform.h"

void irdwt_convolution(double *x_out, int lx, double *g0, double *g1, int lh, double *x_inl, double *x_inh) {
  int i, j;
  double x0;

  for (i=lh-2; i > -1; i--){
    x_inl[i] = x_inl[lx+i];
    x_inh[i] = x_inh[lx+i];
  }
  for (i=0; i<lx; i++){
    x0 = 0;
    for (j=0; j<lh; j++)
      x0 = x0 + x_inl[j+i]*g0[lh-1-j] +
	x_inh[j+i]*g1[lh-1-j];
    x_out[i] = x0;
  }
}


void irdwt_allocate(int m, int n, int lh, double **x_high, double **x_dummy_low, double **x_dummy_high, double **y_dummy_low_low, 
  double **y_dummy_low_high, double **y_dummy_high_low, double **y_dummy_high_high, double **g0, double **g1) {
  *x_high            = (double *) rwt_calloc(m*n,           sizeof(double));
  *x_dummy_low       = (double *) rwt_calloc(max(m,n),      sizeof(double));
  *x_dummy_high      = (double *) rwt_calloc(max(m,n),      sizeof(double));
  *y_dummy_low_low   = (double *) rwt_calloc(max(m,n)+lh-1, sizeof(double));
  *y_dummy_low_high  = (double *) rwt_calloc(max(m,n)+lh-1, sizeof(double));
  *y_dummy_high_low  = (double *) rwt_calloc(max(m,n)+lh-1, sizeof(double));
  *y_dummy_high_high = (double *) rwt_calloc(max(m,n)+lh-1, sizeof(double));
  *g0                = (double *) rwt_calloc(lh,            sizeof(double));
  *g1                = (double *) rwt_calloc(lh,            sizeof(double));
}


void irdwt_free(double **x_dummy_low, double **x_dummy_high, double **y_dummy_low_low, double **y_dummy_low_high, 
  double **y_dummy_high_low, double **y_dummy_high_high, double **g0, double **g1) {
  rwt_free(*x_dummy_low);
  rwt_free(*x_dummy_high);
  rwt_free(*y_dummy_low_low);
  rwt_free(*y_dummy_low_high);
  rwt_free(*y_dummy_high_low);
  rwt_free(*y_dummy_high_high);
  rwt_free(*g0);
  rwt_free(*g1);
}

/* not the same as idwt_coefficients */
void irdwt_coefficients(int lh, double *h, double **g0, double **g1) {
  int i;
  for (i=0; i<lh; i++){
    (*g0)[i] = h[i]/2;
    (*g1)[i] = h[lh-i-1]/2;
  }
  for (i=1; i<=lh; i+=2)
    (*g1)[i] = -((*g1)[i]);
}


void irdwt(double *x, int m, int n, double *h, int lh, int L, double *y_low, double *y_high) {
  double  *g0, *g1, *y_dummy_low_low, *y_dummy_low_high, *y_dummy_high_low;
  double *y_dummy_high_high, *x_dummy_low , *x_dummy_high, *x_high;
  long i;
  int actual_L, actual_m, actual_n, column_of_a, ir, n_c, n_cb, lh_minus_one;
  int ic, n_r, n_rb, column_of_a_plus_double_n, sample_f;

  irdwt_allocate(m, n, lh, &x_high, &x_dummy_low, &x_dummy_high, &y_dummy_low_low, 
    &y_dummy_low_high, &y_dummy_high_low, &y_dummy_high_high, &g0, &g1);
  irdwt_coefficients(lh, h, &g0, &g1);
  
  if (n==1){
    n = m;
    m = 1;
  }
  /* analysis lowpass and highpass */
  
  lh_minus_one = lh - 1;
  /* 2^L */
  sample_f = 1;
  for (i=1; i<L; i++)
    sample_f = sample_f*2;
  actual_m = m/sample_f;
  actual_n = n/sample_f;
  /* restore y_low in x */
  for (i=0;i<m*n;i++)
    x[i] = y_low[i];
  
  /* main loop */
  for (actual_L=L; actual_L >= 1; actual_L--){
    /* actual (level dependent) column offset */
    if (m==1)
      column_of_a = n*(actual_L-1);
    else
      column_of_a = 3*n*(actual_L-1);
    column_of_a_plus_double_n = column_of_a + 2*n;
    
    /* go by columns in case of a 2D signal*/
    if (m>1){
      n_rb = m/actual_m;                 /* # of row blocks per column */
      for (ic=0; ic<n; ic++){            /* loop over column */
	for (n_r=0; n_r<n_rb; n_r++){    /* loop within one column */
	  /* store in dummy variables */
	  ir = -sample_f + n_r;
	  for (i=0; i<actual_m; i++){    
	    ir = ir + sample_f;
	    y_dummy_low_low[i+lh_minus_one]   = mat(x,                                  ir, ic, m, n);  
	    y_dummy_low_high[i+lh_minus_one]  = mat(y_high + m * (column_of_a        ), ir, ic, m, n);  
	    y_dummy_high_low[i+lh_minus_one]  = mat(y_high + m * (column_of_a + n    ), ir, ic, m, n);  
	    y_dummy_high_high[i+lh_minus_one] = mat(y_high + m * (column_of_a + 2 * n), ir, ic, m, n);   
	  }
	  /* perform filtering and adding: first LL/LH, then HL/HH */
	  irdwt_convolution(x_dummy_low,  actual_m, g0, g1, lh, y_dummy_low_low,  y_dummy_low_high); 
	  irdwt_convolution(x_dummy_high, actual_m, g0, g1, lh, y_dummy_high_low, y_dummy_high_high); 
	  /* store dummy variables in matrices */
	  ir = -sample_f + n_r;
	  for (i=0; i<actual_m; i++){    
	    ir = ir + sample_f;
	    mat(x,      ir, ic, m, n) = x_dummy_low[i];  
	    mat(x_high, ir, ic, m, n) = x_dummy_high[i];  
	  }
	}
      }
    }
    
    /* go by rows */
    n_cb = n/actual_n;                 /* # of column blocks per row */
    for (ir=0; ir<m; ir++){            /* loop over rows */
      for (n_c=0; n_c<n_cb; n_c++){    /* loop within one row */      
	/* store in dummy variable */
	ic = -sample_f + n_c;
	for  (i=0; i<actual_n; i++){    
	  ic = ic + sample_f;
	  y_dummy_low_low[i+lh_minus_one] = mat(x, ir, ic, m, n);  
	  if (m>1)
	    y_dummy_high_high[i+lh_minus_one] = mat(x_high, ir, ic, m, n);
	  else
            y_dummy_high_high[i+lh_minus_one] = mat(y_high + m * column_of_a, ir, ic, m, n);
	} 
	/* perform filtering lowpass/highpass */
	irdwt_convolution(x_dummy_low, actual_n, g0, g1, lh, y_dummy_low_low, y_dummy_high_high); 
	/* restore dummy variables in matrices */
	ic = -sample_f + n_c;
	for (i=0; i<actual_n; i++){    
	  ic = ic + sample_f;
	  mat(x, ir, ic, m, n) = x_dummy_low[i];  
	}
      }
    }
    sample_f = sample_f/2;
    actual_m = actual_m*2;
    actual_n = actual_n*2;
  }
  irdwt_free(&x_dummy_low, &x_dummy_high, &y_dummy_low_low, &y_dummy_low_high, &y_dummy_high_low, &y_dummy_high_high, &g0, &g1);
}

