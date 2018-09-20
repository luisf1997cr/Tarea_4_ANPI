/**
 * Copyright (C) 2017
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @author Pablo Alvarado
 * @date   15.12.2017
 */

#ifndef ANPI_SOLVER_HPP
#define ANPI_SOLVER_HPP

#include "LUDoolittle.hpp"

namespace anpi {

  /** faster method used for LU decomposition
   */
  template<typename T>
  inline void lu(const anpi::Matrix<T>& A,
                 anpi::Matrix<T> LU,
                 std::vector<size_t>& p) {
    anpi::luDoolittle(A,LU,p);
  }

  /** method used to create  the permutation matrix given a
   * permutation vector
  **/
  template<typename T>
  void permutationMatrix(const std::vector<size_t>&       p,
                         anpi::Matrix<T>&           pMatrix) {

    int n = p.size();
    anpi::Matrix<T> matrix(p.size(),p.size());

    for(int i = 0; i < n; i++) {
      matrix[i][p[i]] = T(1);
    }

    pMatrix = matrix;

  }

  /// method used to solve lower triangular matrices
  template<typename T>
  void forwardSubstitution(const anpi::Matrix<T>& L,
                           const std::vector<T>&  b,
                           std::vector<T>&        y) {

    int n = L.rows();
    std::vector<T> x(n);
    for(int i = 0; i < n; i++) {
      x[i] =  T(1);
    }

    T sum;

    for(int m = 0; m < n ; m++) {
      sum = T(0);
      for(int i = 0; i < m; i++) {
        sum += L[m][i] * x[i];
      }
      x[m] =  (b[m] - sum)/L[m][m];
    }

    y = x;

  }

  /// method used to solve upper triangular matrices
  template<typename T>
  void backwardSubstitution(const anpi::Matrix<T>& U,
                        const std::vector<T>&  y,
                        std::vector<T>&        x) {
    int n = U.rows();

    std::vector<T> w(n);
    for(int i = 0; i < n; i++) {
      w[i] =  T(1);
    }

    T sum;

    w[n-1] = y[n-1] / U[n-1][n-1];

    for(int i = (n-2); i >= 0; i--) {
      sum = T(0);
      for(int j = (n-1); j >= (i+1); j--) {
        sum += U[i][j] * w[j];
      }
      w[i] = (y[i] - sum) / U[i][i];
    }

    x = w;

  }

  template<typename T>
  bool solveLU(const anpi::Matrix<T>& A,
               std::vector<T>&        x,
               const std::vector<T>&  b) {

    anpi::Matrix<T> LU;
    std::vector<size_t> p;
    anpi::lu(A,LU,p);

    anpi::Matrix<T> L;
    anpi::Matrix<T> U;
    anpi::unpackDoolittle(LU,L,U);

    anpi::Matrix<T> P;
    anpi::permutationMatrix(p,P);

    anpi::Matrix<T>PB = P * b; //ERRROR

    std::vector<T> Pb(PB.rows());

    for(int i = 0; i < PB.rows(); i++) {
      Pb[i] = PB[i][0];
    }

    std::vector<T>y;
    anpi::forwardSubstitution(L,Pb,y);

    anpi::backwardSubstitution(U,y,x);

    return 1;
  }

}//anpi

#endif
