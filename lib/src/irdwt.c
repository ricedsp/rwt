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
  int current_level, current_rows, current_cols, column_cursor, column_blocks_per_row, lh_minus_one;
  int idx_rows, idx_cols, n_r, n_c;
  int row_blocks_per_column, column_cursor_plus_n, column_cursor_plus_double_n, sample_f;
  int three_n_L;

  irdwt_allocate(m, n, lh, &x_high, &x_dummy_low, &x_dummy_high, &y_dummy_low_low, 
    &y_dummy_low_high, &y_dummy_high_low, &y_dummy_high_high, &g0, &g1);
  irdwt_coefficients(lh, h, &g0, &g1);
 
  if (n==1) {
    n = m;
    m = 1;
  }
  /* analysis lowpass and highpass */
  
  three_n_L = 3*n*L;
  lh_minus_one = lh - 1;
  /*! we calculate sample_f = 2^(L - 1) with a loop since that is actually the recommended method for whole numbers */
  sample_f = 1;
  for (i=1; i<L; i++)
    sample_f = sample_f*2;

  current_rows = m/sample_f;
  current_cols = n/sample_f;
  /* restore y_low in x */
  for (i=0; i<m*n; i++)
    x[i] = y_low[i];
  
  /* main loop */
  for (current_level=L; current_level >= 1; current_level--) {
    /* actual (level dependent) column offset */
    if (m==1)
      column_cursor = n*(current_level-1);
    else
      column_cursor = 3*n*(current_level-1);
    column_cursor_plus_n = column_cursor + n;
    column_cursor_plus_double_n = column_cursor_plus_n + n;
    
    /* go by columns in case of a 2D signal*/
    if (m>1) {
      row_blocks_per_column = m/current_rows;                /* # of row blocks per column */
      for (idx_cols=0; idx_cols<n; idx_cols++) {                       /* loop over column */
	for (n_r=0; n_r<row_blocks_per_column; n_r++) {          /* loop within one column */
	  /* store in dummy variables */
	  idx_rows = -sample_f + n_r;
	  for (i=0; i<current_rows; i++) {    
	    idx_rows = idx_rows + sample_f;
	    y_dummy_low_low[i+lh_minus_one]   = mat(x,      idx_rows, idx_cols,                               m, n);
	    y_dummy_low_high[i+lh_minus_one]  = mat(y_high, idx_rows, idx_cols + column_cursor,               m, three_n_L);
	    y_dummy_high_low[i+lh_minus_one]  = mat(y_high, idx_rows, idx_cols + column_cursor_plus_n,        m, three_n_L);
	    y_dummy_high_high[i+lh_minus_one] = mat(y_high, idx_rows, idx_cols + column_cursor_plus_double_n, m, three_n_L);
	  }
	  /* perform filtering and adding: first LL/LH, then HL/HH */
	  irdwt_convolution(x_dummy_low,  current_rows, g0, g1, lh, y_dummy_low_low,  y_dummy_low_high); 
	  irdwt_convolution(x_dummy_high, current_rows, g0, g1, lh, y_dummy_high_low, y_dummy_high_high); 
	  /* store dummy variables in matrices */
	  idx_rows = -sample_f + n_r;
	  for (i=0; i<current_rows; i++) {
	    idx_rows = idx_rows + sample_f;
	    mat(x,      idx_rows, idx_cols, m, n) = x_dummy_low[i];
	    mat(x_high, idx_rows, idx_cols, m, n) = x_dummy_high[i];
	  }
	}
      }
    }
    
    /* go by rows */
    column_blocks_per_row = n/current_cols;                /* # of column blocks per row */
    for (idx_rows=0; idx_rows<m; idx_rows++) {                           /* loop over rows */
      for (n_c=0; n_c<column_blocks_per_row; n_c++) {  /* loop within one row */      
	/* store in dummy variable */
	idx_cols = -sample_f + n_c;
	for  (i=0; i<current_cols; i++) {    
	  idx_cols = idx_cols + sample_f;
	  y_dummy_low_low[i+lh_minus_one] = mat(x, idx_rows, idx_cols, m, n);  
	  if (m>1)
	    y_dummy_high_high[i+lh_minus_one] = mat(x_high, idx_rows, idx_cols,                 m, n);
	  else
            y_dummy_high_high[i+lh_minus_one] = mat(y_high, idx_rows, idx_cols + column_cursor, m, three_n_L);
	} 
	/* perform filtering lowpass/highpass */
	irdwt_convolution(x_dummy_low, current_cols, g0, g1, lh, y_dummy_low_low, y_dummy_high_high); 
	/* restore dummy variables in matrices */
	idx_cols = -sample_f + n_c;
	for (i=0; i<current_cols; i++) {    
	  idx_cols = idx_cols + sample_f;
	  mat(x, idx_rows, idx_cols, m, n) = x_dummy_low[i];  
	}
      }
    }
    sample_f = sample_f/2;
    current_rows = current_rows*2;
    current_cols = current_cols*2;
  }
  irdwt_free(&x_dummy_low, &x_dummy_high, &y_dummy_low_low, &y_dummy_low_high, &y_dummy_high_low, &y_dummy_high_high, &g0, &g1);
}

