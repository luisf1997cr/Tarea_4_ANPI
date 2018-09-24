#include "LUCrout.hpp"

namespace anpi
{

template <typename T>
inline void lu(const anpi::Matrix<T> &A, anpi::Matrix<T> &LU, std::vector<size_t> &p)
{
    //this was benchmarked to be a little faster than Doolittle
    anpi::luCrout(A, LU, p);
}

/*
// Function to solve an equation of the form A x = b
// 
// 
// 
// 
//
*/
template <typename T>
bool solveLU(const anpi::Matrix<T> &A, std::vector<T> &x, std::vector<T> &b)
{
    //size n of the Matrix
    int n = A.cols();

    //check for corect inputs
    if (int(A.rows()) != n) //check if it's square
        throw anpi::Exception("solveLU::Input Matrix A is not square");
    if (int(b.size()) != n) // check the solutions are the same amount as the equations
        throw anpi::Exception("solveLU::solve vector needs to be the same size as the input matrix");

    //the matrix that holds the LU decomposition of A
    anpi::Matrix<T> LU;
    //the permutations vector
    std::vector<size_t> p;
    // copy b into x and set size
    x = b;

    //decompose A into LU
    anpi::lu(A, LU, p);

    int i, ii = 0, ip, j;
    T sum;

    //we permute the b vector to match the permuted A matrix obtained from LU decomposition
    for (i = 0; i < n; ++i)
    {
        ip = p[i];    //index in the permutation vector
        sum = x[ip];  //permuted result from b
        x[ip] = x[i]; //move the current value to the permuted position
        // x[i] = sum;
        if (ii != 0)
            for (j = ii - 1; j < i; j++)
                sum -= LU[i][j] * x[j];
        else if (sum != 0.0)
            ii = i + 1;

        x[i] = sum;
    }

    // //forward substitution
    // for (i = 0; i < n; ++i)
    // {
    //     sum = x[i] / LU[i][i];
    //     for (j = 0; j <= i; ++j)
    //         sum -= LU[i][j] * x[j];
    //     x[i] = sum; //set the assigned
    // }

    // back substitution
    for (i = n - 1; i >= 0; --i)
    {
        sum = x[i];
        for (j = i + 1; j < n; ++j)
            sum -= LU[i][j] * x[j];
        x[i] = sum / LU[i][i]; //Store a component of the solution vector X .
    }
    return true;
} //end of solveLU

template <typename T>
void invert(const anpi::Matrix<T> &A, anpi::Matrix<T> &Ai)
{
    int i, j;

    Ai = A; //copy size of input matrix
    int n = A.cols();
    //create an I (identity) Matrix
    for (i = 0; i < n; ++i)
    {
        for (j = 0; j < n; ++j)
            Ai[i][j] = 0;
        Ai[i][i] = 1;
    }
    //vector to copy each row of the Identity matrix
    std::vector<T> xx(n);
    for (j = 0; j < n; j++)
    { //Copy and solve each column in turn.
        for (i = 0; i < n; i++)
            xx[i] = Ai[i][j];
        anpi::solveLU(A, xx, xx);
        for (i = 0; i < n; i++)
            Ai[i][j] = xx[i];
    }

    // solveLU(Ai, Ai);
}

} // namespace anpi