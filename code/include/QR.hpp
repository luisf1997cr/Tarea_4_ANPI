#include <cmath>
#include <limits>
#include <functional>

#include <iostream>

#include "Exception.hpp"
#include "Matrix.hpp"

#ifndef ANPI_LU_CROUT_HPP
#define ANPI_LU_CROUT_HPP

namespace anpi
{

template <typename T>
void qr(const anpi::Matrix<T> &A, anpi::Matrix<T> &Q, anpi::Matrix<T> &R)
{
    //start implementation

    // QRdcmp::QRdcmp(MatDoub_I &a)
    // : n(a.nrows()), qt(n,n), r(a), sing(false) {
    // Construct the QR decomposition of a[0..n-1][0..n-1]. The upper triangular matrix R and
    // the transpose of the orthogonal matrix Q are stored. sing is set to true if a singularity is
    // encountered during the decomposition, but the decomposition is still completed in this case;
    // otherwise it is set to false.

    if (A.rows() != A.cols())
    {
        throw anpi::Exception("Matrix for QR decomposition must be square");
    }
    int n = A.cols(); //size of matrix

    // initialize matrices
    R = A;
    Q.clear();
    Q.allocate(n, n);

    int i, j, k;
    std::vector<T> c(n), d(n);
    T scale, sigma, sum, tau, tempT;

    for (k = 0; k < n - 1; ++k)
    {
        scale = 0.0;
        for (i = k; i < n; ++i)
            scale = scale > abs(R[i][k]) ? scale : abs(R[i][k]); //returns the maximum value between the 2 arguments
        if (scale == 0.0)
        { //Singular case.
            // sing=true;
            c[k] = d[k] = 0.0;
        }
        else
        { //Form Q k and Q k  A.
            for (i = k; i < n; i++)
            {
                R[i][k] /= scale;
            }
            for (sum = 0.0, i = k; i < n; ++i)
            {
                sum += (R[i][k] * R[i][k]);
            }

            tempT = std::sqrt(sum);
            sigma = (tempT > 0 && R[k][k] > 0) || (tempT < 0 && R[k][k] < 0) ? tempT : T(-1) * tempT; //gives tempT the same sign as the value in r[k][k]

            R[k][k] += sigma;
            c[k] = sigma * R[k][k];
            d[k] = -scale * sigma;
            for (j = k + 1; j < n; ++j)
            {
                for (sum = 0.0, i = k; i < n; ++i)
                {
                    sum += R[i][k] * R[i][j];
                }
                tau = sum / c[k];
                for (i = k; i < n; ++i)
                {
                    R[i][j] -= tau * R[i][k];
                }
            }
        }
    }
    d[n - 1] = R[n - 1][n - 1];
    // if (d[n-1] == 0.0) sing=true;
    for (i = 0; i < n; i++)
    { //Form Q T explicitly.
        for (j = 0; j < n; j++)
            Q[i][j] = 0.0;
        Q[i][i] = 1.0;
    }
    for (k = 0; k < n - 1; k++)
    {
        if (c[k] != 0.0)
        {
            for (j = 0; j < n; j++)
            {
                sum = 0.0;
                for (i = k; i < n; i++)
                    sum += R[i][k] * Q[i][j];
                sum /= c[k];
                for (i = k; i < n; i++)
                    Q[i][j] -= sum * R[i][k];
            }
        }
    }
    for (i = 0; i < n; i++)
    { //Form R explicitly.
        R[i][i] = d[i];
        for (j = 0; j < i; j++)
            R[i][j] = 0.0;
    }

} //end QR decomposition

template <class T>
void printMatrix(anpi::Matrix<T> &A)
{
    int n = A.cols();
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            std::cout << A[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

template <class T>
void transposeMatrix(anpi::Matrix<T> &A)
{
    T swap;
    int n = A.cols();
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            swap = A[i][j];
            A[i][j] = A[j][i];
            A[j][i] = swap;
        }
        std::cout << std::endl;
    }
}
} // namespace anpi

#endif