#include <iostream>
#include "Matrix.hpp"

namespace anpi
{

//---------------------------------------helping functions Tarea 4-----------------------------------//
template <class T>
void printMatrix(anpi::Matrix<T> &A)
{
    int n = A.cols();
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            std::cout << A[i][j] << "  ";
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
        for (int j = 0; j < i; ++j)
        {
            swap = A[i][j];
            A[i][j] = A[j][i];
            A[j][i] = swap;
        }
    }
}

} // namespace anpi