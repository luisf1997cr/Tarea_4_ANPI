/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @Author: 
 * @Date  : 03.03.2018
 */

#include <cmath>
#include <limits>
#include <functional>

#include "Exception.hpp"
#include "Matrix.hpp"

#ifndef ANPI_LU_CROUT_HPP
#define ANPI_LU_CROUT_HPP

namespace anpi
{

/**
   * Auxiliary method used to debug LU decomposition.
   *
   * It separates a packed LU matrix into the lower triangular matrix
   * L and the upper triangular matrix U, such that the diagonal of U
   * is composed by 1's.
   */
template <typename T>
void unpackCrout(const Matrix<T> &LU,
                 Matrix<T> &L,
                 Matrix<T> &U)
{
  if (LU.rows() != LU.cols())
  {
    throw anpi::Exception("Matrix for Crout unpacking must be square");
  }

  L = LU;
  U = LU;
  int size = LU.cols();
  //we make all upper values of L equal 0
  for (int i = 0; i < size; ++i)
  {
    L[i][i] = 1;
    for (int j = i + 1; j < size; ++j)
    {
      L[i][j] = 0;
    }
  }
  //we make all lower values of U equal 0, and the diagonal equal 1
  for (int i = 0; i < size; ++i)
  {
    // U[i][i] = 1;
    for (int j = 0; j < i; ++j)
    {
      U[i][j] = 0;
    }
  }

  // throw anpi::Exception("To be implemented yet");
}

/**
   * Decompose the matrix A into a lower triangular matrix L and an
   * upper triangular matrix U.  The matrices L and U are packed into
   * a single matrix LU.  
   *
   * Crout's way of packing assumes a diagonal of
   * 1's in the U matrix.
   *
   * @param[in] A a square matrix 
   * @param[out] LU matrix encoding the L and U matrices
   * @param[out] permut permutation vector, holding the indices of the
   *             original matrix falling into the corresponding element.
   *             For example if permut[5]==3 holds, then the fifth row
   *             of the LU decomposition in fact is dealing with the third
   *             row of the original matrix.
   *
   * @throws anpi::Exception if matrix cannot be decomposed, or input
   *         matrix is not square.
   */
template <typename T>
void luCrout(const Matrix<T> &A,
             Matrix<T> &LU,
             std::vector<size_t> &permut)
{

  if (A.rows() != A.cols())
  {
    throw anpi::Exception("Matrix for Crout LU decomposition must be square");
  }

  // LUdcmp::LUdcmp(MatDoub_I &a) : n(a.nrows()), lu(a), aref(a), indx(n) {
  // Given a matrix a[0..n-1][0..n-1], this routine replaces it by the LU decomposition of a
  // rowwise permutation of itself. a is input. On output, it is arranged as in equation (2.3.14)
  // above; indx[0..n-1] is an output vector that records the row permutation effected by the
  // partial pivoting; d is output as  ̇1 depending on whether the number of row interchanges
  // was even or odd, respectively. This routine is used in combination with solve to solve linear
  // equations or invert a matrix.

  //intialize the result matrix, permutation vector
  LU = A;
  int n = A.rows();
  permut.resize(n);

  for (int i = 0; i < n; ++i)
    permut[i] = i;

  const T TINY = T(1.0e-40); //A small number.
  int i, imax, j, k;
  T big, temp;
  std::vector<T> vv(n); //vv stores the implicit scaling of each row.
  // d = 1.0;              //No row interchanges yet.
  for (i = 0; i < n; i++)
  { //Loop over rows to get the implicit scaling info
    big = 0.0;
    for (j = 0; j < n; j++)
      if ((temp = abs(LU[i][j])) > big)
        big = temp;
    if (big == T(0.0))
      throw anpi::Exception("A is a singular matrix, unable to decompose into LU");

    //No nonzero largest element.
    vv[i] = T(1.0) / big; //Save the scaling.
  }
  for (k = 0; k < n; k++)
  { //This is the outermost kij loop.
    big = 0.0;
    //Initialize for the search for largest pivot element.
    for (i = k; i < n; i++)
    {
      temp = vv[i] * abs(LU[i][k]);
      if (temp > big) //Is the figure of merit for the pivot better than the best so far?
      {
        big = temp;

        imax = i;
      }
    }
    if (k != imax)
    { //Do we need to interchange rows?
      for (j = 0; j < n; j++)
      { //Yes, do so...
        temp = LU[imax][j];
        LU[imax][j] = LU[k][j];
        LU[k][j] = temp;
      }
      // d = -d;
      temp = vv[imax];  //...and change the parity of d.
      vv[imax] = vv[k]; //Also interchange the scale factor.
      vv[k] = temp;

      //do the change
      temp = permut[k];
      permut[k] = imax;
      permut[imax] = temp;
    }

    if (LU[k][k] == T(0.0))
      LU[k][k] = TINY;
    // If the pivot element is zero, the matrix is singular (at least to the precision of the
    // algorithm). For some applications on singular matrices, it is desirable to substitute2.3 LU Decomposition and Its Applications
    // 53
    // TINY for zero.
    for (i = k + 1; i < n; i++)
    {
      temp = LU[i][k] /= LU[k][k]; // Divide by the pivot element.
      for (j = k + 1; j < n; j++)  //Innermost loop: reduce remaining submatrix.
        LU[i][j] -= temp * LU[k][j];
    }
  }
} //end LUCrout

// int mxSize = A.rows();
// //initialize vector holding the pivot, which shows the corresponding place where a row
// //of the original matrix was transposed
// std::vector<size_t> pivot;
// pivot.resize(mxSize);
// for (int i = 0; i < mxSize; ++i)
// {
//   pivot[i] = i;
// }
} // namespace anpi
// namespace anpi

#endif
