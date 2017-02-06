/*! \file idwt.cc
    \brief Implementation of the inverse discrete wavelet transform

*/
#include "rwt_platform.h"
#include "rwt_transforms.h"
#include <vector>
#include <algorithm>

namespace {
template <typename data>
class IDWT
{
  private:
    size_t nrows;
    size_t ncols;
    int ncoeff;
    int levels;
    std::vector<data> x_dummy;
    std::vector<data> y_dummy_low; // low pass results of convolution
    std::vector<data> y_dummy_high;// high pass results of convolution
    std::vector<data> coeff_low;///< the low pass coefficients - reversed h
    std::vector<data> coeff_high;///< the high pass coefficients - forward h, alternate values are sign flipped

  public:
    IDWT(size_t nrows_,size_t ncols_,const data * h,int ncoeff_,int levels_)
      :nrows(nrows_),ncols(ncols_), ncoeff(ncoeff_),levels(levels_),
      x_dummy(std::max(nrows,ncols)),
      y_dummy_low(std::max(nrows,ncols)+ncoeff/2-1),
      y_dummy_high(std::max(nrows,ncols)+ncoeff/2-1),
      coeff_low(ncoeff),
      coeff_high(ncoeff)
    {
      if (ncols==1)
        std::swap(nrows,ncols);
      idwt_coefficients(h);
    }

    /*!
     * Put the scaling coeffients into a form ready for use in the convolution function
     * @param h  the wavelet scaling coefficients
     */
    void idwt_coefficients(const data *h)
    {
      for (int i=0; i<ncoeff; i++) {
        coeff_low[i] = h[i];
        coeff_high[i] = h[ncoeff-i-1];
      }
      for (int i=1; i<=ncoeff; i+=2)
        coeff_high[i] = -coeff_high[i];
    }

    void idwt_convolution(size_t lx)
    {
      const int ncoeff_minus_one = ncoeff - 1;
      const int ncoeff_halved_minus_one = ncoeff/2 - 1;
      int k;
      size_t i, j, ind, tj;
      data x0, x1;

      for (k=ncoeff_halved_minus_one-1; k > -1; k--) {
        y_dummy_low[k]  = y_dummy_low[lx+k];
        y_dummy_high[k] = y_dummy_high[lx+k];
      }

      ind = 0;
      for (i=0; i<(lx); i++) {
        x0 = 0;
        x1 = 0;
        tj = 0;
        for (j=0; j<=ncoeff_halved_minus_one; j++) {
          x0 = x0 + (y_dummy_low[i+j] * coeff_low[ncoeff_minus_one-1-tj]) + (y_dummy_high[i+j] * coeff_high[ncoeff_minus_one-1-tj]);
          x1 = x1 + (y_dummy_low[i+j] * coeff_low[ncoeff_minus_one-tj])   + (y_dummy_high[i+j] * coeff_high[ncoeff_minus_one-tj]);
          tj += 2;
        }
        x_dummy[ind++] = x0;
        x_dummy[ind++] = x1;
      }
    }

  /*!
   * Perform the inverse discrete wavelet transform
   *
   * @param x      the output signal with the inverse wavelet transform applied
   * @param y      the input signal
   */
  void process(data *x, const data *y)
  {
    long i;
    int current_level, sample_f;
    size_t current_rows, current_cols, row_cursor, column_cursor, idx_rows, idx_cols;

    int ncoeff_halved_minus_one = ncoeff/2 - 1;
    /* 2^levels */
    sample_f = 1;
    for (i=1; i<levels; i++)
      sample_f = sample_f*2;

    if (nrows>1)
      current_rows = nrows/sample_f;
    else
      current_rows = 1;
    current_cols = ncols/sample_f;

    for (i=0; i<(nrows*ncols); i++)
      x[i] = y[i];

    /* main loop */
    for (current_level=levels; current_level >= 1; current_level--) {
      row_cursor = current_rows/2;
      column_cursor = current_cols/2;

      /* go by columns in case of a 2D signal*/
      if (nrows>1) {
        for (idx_cols=0; idx_cols<current_cols; idx_cols++) {         /* loop over columns */
          /* store in dummy variables */
          idx_rows = row_cursor;
          for (i=0; i<row_cursor; i++){
            y_dummy_low[i+ncoeff_halved_minus_one]  = mat(x, i,          idx_cols, nrows, ncols);
            y_dummy_high[i+ncoeff_halved_minus_one] = mat(x, idx_rows++, idx_cols, nrows, ncols);
          }
          /* perform filtering lowpass and highpass*/
          idwt_convolution(row_cursor);
          /* restore dummy variables in matrix */
          for (i=0; i<current_rows; i++)
            mat(x, i, idx_cols, nrows, ncols) = x_dummy[i];
        }
      }
      /* go by rows */
      for (idx_rows=0; idx_rows<current_rows; idx_rows++) {           /* loop over rows */
        /* store in dummy variable */
        idx_cols = column_cursor;
        for  (i=0; i<column_cursor; i++){
          y_dummy_low[i+ncoeff_halved_minus_one]  = mat(x, idx_rows, i,          nrows, ncols);
          y_dummy_high[i+ncoeff_halved_minus_one] = mat(x, idx_rows, idx_cols++, nrows, ncols);
        }
        /* perform filtering lowpass and highpass*/
        idwt_convolution(column_cursor);
        /* restore dummy variables in matrices */
        for (i=0; i<current_cols; i++)
          mat(x, idx_rows, i, nrows, ncols) = x_dummy[i];
      }
      if (nrows==1)
        current_rows = 1;
      else
        current_rows = current_rows*2;
      current_cols = current_cols*2;
    }
  }
};
}
void idwt_double(double *x, size_t nrows, size_t ncols, double *h, int ncoeff, int levels, double *y)
{
    IDWT<double>(nrows,ncols,h,ncoeff,levels).process(x,y);
}

void idwt_float(float *x, size_t nrows, size_t ncols, float *h, int ncoeff, int levels, float *y)
{
    IDWT<float>(nrows,ncols,h,ncoeff,levels).process(x,y);
}
