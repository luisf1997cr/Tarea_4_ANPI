/**
 * Copyright (C) 2017 
 * Área Académica de Ingeniería en Computadoras, TEC, Costa Rica
 *
 * This file is part of the CE3102 Numerical Analysis lecture at TEC
 *
 * @Author: Pablo Alvarado
 * @Date  : 10.02.2018
 */

#include <boost/test/unit_test.hpp>

#include "QR.hpp"
// #include "LUDoolittle.hpp"

#include <iostream>
#include <exception>
#include <cstdlib>
#include <complex>

#include <functional>

#include <cmath>

namespace anpi
{
namespace test
{

/// Test the given closed root finder
template <typename T>
void qrTest()
{

    // The result
    Matrix<T> Q, R, QT, mult;

    // Test decomposition
    {
        anpi::Matrix<T> A = {{2, 0, 1, 2}, {-1, -2, 1, 2}, {1, 1, 1, 1}, {-1, -1, 0, 1}};
        anpi::Matrix<T> exQ = {{-0.63246, 0.62280, 0.38139, 0.25820},
                               {-0.63246, -0.54495, 0.19069, -0.5164},
                               {-0.31623, 0.31140, -0.85812, -0.25820},
                               {-0.31623, -0.46710, -0.28604, 0.77460}};

        anpi::Matrix<T> exR = {{-3.16228, -0.63246, 1.26491, -1.58114},
                               {0.00000, 2.56905, 1.86840, 0.38925},
                               {0.00000, 0.00000, -0.95346, -0.28604},
                               {0.00000, 0.00000, 0.00000, -0.51640}};
        // std::vector<size_t> p;

        qr(A, Q, R);
        std::cout << "the matrix A is:\n";
        anpi::printMatrix(A);
        std::cout << std::endl;
        std::cout << "the matrix Q is:\n";
        anpi::printMatrix(Q);
        std::cout << std::endl;
        std::cout << "the matrix R is:\n";
        anpi::printMatrix(R);

        QT = Q;
        anpi::transposeMatrix(QT);
        std::cout << std::endl;
        std::cout << "the matrix QT is:\n";
        anpi::printMatrix(QT);

        std::cout << std::endl;
        std::cout << "QT * Q is:\n";
        mult = Q * QT;
        anpi::printMatrix(mult);

        Matrix<T> Ar = Q * R;

        const T eps = 2.2e-5; //std::numeric_limits<T>::epsilon();

        BOOST_CHECK(Ar.rows() == A.rows());
        BOOST_CHECK(Ar.cols() == A.cols());

        for (size_t i = 0; i < Ar.rows(); ++i)
        {
            for (size_t j = 0; j < Ar.cols(); ++j)
            {
                BOOST_CHECK(std::abs(Ar(i, j) - A(i, j)) < eps);
            }
        }

        for (size_t i = 0; i < Ar.rows(); ++i)
        {
            for (size_t j = 0; j < Ar.cols(); ++j)
            {
                BOOST_CHECK(std::abs(Q(i, j) - exQ(i, j)) < eps);
            }
        }

        for (size_t i = 0; i < Ar.rows(); ++i)
        {
            for (size_t j = 0; j < Ar.cols(); ++j)
            {
                BOOST_CHECK(std::abs(R(i, j) - exR(i, j)) < eps);
            }
        }
    }
}

} // namespace test
} // namespace anpi

BOOST_AUTO_TEST_SUITE(QR)

BOOST_AUTO_TEST_CASE(simpleQR)
{
    anpi::test::qrTest<float>();
    anpi::test::qrTest<double>();
}

BOOST_AUTO_TEST_SUITE_END()
