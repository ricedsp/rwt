/*! \file dwt.cc
    \brief Implementation of the discrete wavelet transform

*/

#include "rwt_platform.h"
#include "rwt_transforms.h"
#include <vector>
#include <algorithm>

namespace {
template <typename data>
class DWT
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
    DWT(size_t nrows_,size_t ncols_,const data * h,int ncoeff_,int levels_)
      :nrows(nrows_),ncols(ncols_), ncoeff(ncoeff_),levels(levels_),
      x_dummy(std::max(nrows,ncols)+ncoeff-1),
      y_dummy_low(std::max(nrows,ncols)),
      y_dummy_high(std::max(nrows,ncols)),
      coeff_low(ncoeff),
      coeff_high(ncoeff)
  {
    if (ncols==1)
      std::swap(nrows,ncols);
    dwt_coefficients(h);
  }

    /*!
     * Perform convolution for dwt
     *
     * @param x_in input signal values
     * @param lx the length of x_in
     *
     * For the convolution we will calculate the output of the lowpass and highpass filters in parallel
     *
     * Normally we can describe the calculation of a convolution as
     * \f$ (\textbf{w} * \textbf{z})_k = \frac{1}{N} \sum\limits_{l=0}^{2N-1} w_{k-l} \cdot z_{l} \f$
     */
    void dwt_convolution(data *x_in, size_t lx)
    {
      int ncoeff_minus_one = ncoeff-1;
      size_t i, j, ind;
      data x0, x1;
      for (i=lx; i<lx+ncoeff_minus_one; i++) {
        x_in[i] = *(x_in+(i-lx)); /*! extend x_in by creating a small mirror at the end of length ncoeff_minus_one */
      }
      ind = 0;
      for (i=0; i<(lx); i+=2) {   /*! Step through the input values, moving right 2 values each loop */
        x0 = 0;
        x1 = 0;
        for (j=0; j<=ncoeff_minus_one; j++) {                   /*! Take the high and low filters in reverse order */
          x0 = x0 + x_in[i+j] * coeff_low[ncoeff_minus_one-j];  /*! Sum the product of the next ncoeff values of x_in with the filter coefficients */
          x1 = x1 + x_in[i+j] * coeff_high[ncoeff_minus_one-j];
        }
        y_dummy_low[ind] = x0; /*! Place these calculated sums in the next position of the output */
        y_dummy_high[ind++] = x1;
      }
    }

    /*!
     * Put the scaling coeffients into a form ready for use in the convolution function
     *
     * @param h  the wavelet scaling coefficients
     *
     * The coefficients of our Quadrature Mirror Filter are described by
     * \f$ g\left[lh - 1 - n \right] = (-1)^n * h\left[n\right] \f$
     */
    void dwt_coefficients(const data *h)
    {
      int i;
      for (i=0; i<ncoeff; i++) {
        coeff_low[i] = h[(ncoeff-i)-1];
        coeff_high[i] = h[i];
      }
      for (i=0; i<ncoeff; i+=2)
        coeff_high[i] = -coeff_high[i];
    }

    /*!
     * Perform the discrete wavelet transform
     *
     * @param x      the input signal
     * @param y      the output signal with the wavelet transform applied
     *
     * The discrete wavelet transform begins with a set of samples of a signal whose length
     * is a power of 2. This exponent will be the maximum number of levels of the transform
     * that we can perform.
     *
     */
    void process(const data *x, data *y)
    {
      const int ncoeff_minus_one = ncoeff-1;

      size_t current_rows = 2*nrows; /*! current_rows and current_cols start at 2x since we divide by 2 at the start of the loop */
      size_t current_cols = 2*ncols;

      for (int current_level=1; current_level<=levels; current_level++ )
      {
        if (nrows==1)
          current_rows = 1;
        else
          current_rows = current_rows/2;
        current_cols = current_cols/2;

        for (size_t idx_rows=0; idx_rows<current_rows; idx_rows++) {
          if (current_level==1)
            for (size_t i=0; i<current_cols; i++)
              x_dummy[i] = mat(x, idx_rows, i, nrows, ncols);
          else
            for (size_t i=0; i<current_cols; i++)
              x_dummy[i] = mat(y, idx_rows, i, nrows, ncols);
          /*! Perform filtering lowpass and highpass*/
          dwt_convolution(&x_dummy[0], current_cols);
          /*! Restore dummy variables in matrices */
          const size_t column_cursor = current_cols/2;
          size_t idx_columns = column_cursor;
          for (size_t i=0; i<column_cursor; i++) {
            mat(y, idx_rows, i,             nrows, ncols) = y_dummy_low[i];
            mat(y, idx_rows, idx_columns++, nrows, ncols) = y_dummy_high[i];
          }
        }

        /*! For the 2d transform, we go through each of the columns after having gone through the rows */
        if (nrows>1) {
          for (size_t idx_columns=0; idx_columns<current_cols; idx_columns++) { /* loop over columns */
            /*! Store in dummy variables */
            for (size_t i=0; i<current_rows; i++)
              x_dummy[i] = mat(y, i, idx_columns, nrows, ncols);
            /*! Perform filtering lowpass and highpass*/
            dwt_convolution(&x_dummy[0], current_rows);
            /*! Restore dummy variables in matrix */
            size_t idx_rows = current_rows/2;
            const size_t row_cursor = current_rows/2;
            for (size_t i=0; i<row_cursor; i++) {
              mat(y, i,          idx_columns, nrows, ncols) = y_dummy_low[i];
              mat(y, idx_rows++, idx_columns, nrows, ncols) = y_dummy_high[i];
            }
          }
        }
      }
    }
};
}
void dwt_double(const double *x, size_t nrows, size_t ncols, double *h, int ncoeff, int levels, double *y)
{
    DWT<double>(nrows,ncols,h,ncoeff,levels).process(x,y);
}

void dwt_float(const float *x, size_t nrows, size_t ncols, float *h, int ncoeff, int levels, float *y)
{
    DWT<float>(nrows,ncols,h,ncoeff,levels).process(x,y);
}
