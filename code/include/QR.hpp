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

//begin LUIS
/*
    size_t k = 0; //Variable while dicta cada ciclo del qr
    anpi::Matrix<T> Atemp = A;
    while (k < 2)
    {
        std::vector<T> fila;
        for (size_t i = 0 + k; i < A.rows(); i++)
        {
            fila.push_back(Atemp[i][k]);
        }

        T norma = T(0);
        for (size_t i = 0; i < fila.size(); i++)
        {
            norma = norma + fila[i] * fila[i];
        }
        norma = std::sqrt(norma);

        std::vector<T> uVector;
        uVector.push_back(fila[0] - norma);
        T uNorma = (fila[0] - norma) * (fila[0] - norma);
        for (size_t i = 1; i < fila.size(); i++)
        {
            uVector[i] = fila[i];
            uNorma = uNorma + fila[i] * fila[i];
        }
        uNorma = std::sqrt(uNorma);

        std::vector<T> vVector;
        for (size_t i = 0; i < fila.size(); i++)
        {
            vVector.push_back(7);
            vVector[i] = (1 / uNorma) * uVector[i];
        }

        ///Verifica la identidad del Qr
        //Verifica el lugar exacto en donde se debe sumar el valor de K
        anpi::Matrix<T> identidad(A.rows(), A.cols());
        for (size_t i = 0 + k; i < identidad.cols(); i++)
        {
            for (size_t j = 0; j < identidad.rows(); j++)
            {
                if (i == j)
                {
                    identidad[i][j] = 1;
                }
                else
                {
                    identidad[i][j] = 0;
                }
            }
        }

        ///METODO Q = I - 2VV^T
        anpi::Matrix<T> Qtemp(A.rows(), A.cols());
        //Verifica y resuelve el produto externo de la matriz
        anpi::Matrix<T> multiVVT(vVector.size(), vVector.size());
        for (size_t i = 0; i < multiVVT.rows(); i++)
        {
            for (size_t j = 0; j < multiVVT.cols(); j++)
            {
                multiVVT[i][j] = vVector[i] * vVector[j];
            }
        }

        //Escala la nueva matriz
        T escalar = T(2);
        for (size_t j = 0; j < multiVVT.rows(); j++)
        {
            for (size_t i = 0; i < multiVVT.cols(); i++)
            {
                multiVVT[i][j] = escalar * multiVVT[i][j];
            }
        }
        if (k > 0)
        {
            //Aunmenta el tama;o de la matriz
            anpi::Matrix<T> agrandada(multiVVT.rows() + 1, multiVVT.cols() + 1); ///VER BIEN DONDE SUMAR EL K
            for (size_t i = 0; i < agrandada.cols(); i++)
            {
                for (size_t j = 0; j < agrandada.rows(); j++)
                {
                    if (i < k or j < k)
                    {
                        if (i == j)
                        {
                            agrandada[i][j] = -1; /// PARA QUE Q2 SEA POSITIVA, ESTA MAL?
                        }
                        else
                        {
                            agrandada[i][j] = 0;
                        }
                    }
                    else
                    {
                        agrandada[i][j] = multiVVT[i - k][j - k];
                    }
                }
            }
            multiVVT = agrandada;
        }

        Qtemp = identidad - multiVVT;

        if (k == 0)
        {
            Q = Qtemp;
        }
        else
        {
            Q = Q * Qtemp;
        }
        Atemp = Qtemp * A;
        ///Achica matriz
        for (size_t i = 0; i < Atemp.rows(); i++)
        {
            if (i != k)
            {
                Atemp[i][k] = 0;
            }
            else
            {
                Atemp[i][k] = 1;
            }
        }
        for (size_t i = 1; i < Atemp.cols(); i++)
        {
            Atemp[k][i] = 0;
        } ///fin achica matriz

        k = k + 1;
    }

    ///METODO transpuesta
    anpi::Matrix<T> transpuestaQ(Q.cols(), Q.rows());
    for (size_t i = 0; i < Q.cols(); i++)
    {
        for (size_t j = 0; j < Q.rows(); j++)
        {
            transpuestaQ[i][j] = Q[j][i];
        }
    }

    R = transpuestaQ * A;

    */
// } // namespace anpi

template <typename T>
bool solveQR(const anpi::Matrix<T> &A, std::vector<T> &x, const std::vector<T> &b)
{

    Matrix<T> Q;
    Matrix<T> R;
    anpi::qr(A, Q, R);

    ///METODO transpuesta
    anpi::Matrix<T> transpuestaQ(Q.cols(), Q.rows());
    for (size_t i = 0; i < Q.cols(); i++)
    {
        for (size_t j = 0; j < Q.rows(); j++)
        {
            transpuestaQ[i][j] = Q[j][i];
        }
    }

    Matrix<T> Rx = transpuestaQ * b; ///vector

    ///METODO sustitucion hacia atras
    x.push_back(7);
    x[R.cols() - 1] = (Rx[Rx.rows() - 1][0]) / R[R.rows() - 1][R.cols() - 1];

    int i = R.cols() - 1;
    while (i > 0)
    {
        i = i - 1;
        T sumatoria = 0;
        for (size_t j = i + 1; j < R.cols(); j++)
        {
            sumatoria = sumatoria + R[i][j] * x[j];
        }
        x[i] = (1 / R[i][i]) * (Rx[i][0] - sumatoria);
    }

    return 1;
}
} // namespace anpi

#endif