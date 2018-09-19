
#include "LUCrout.hpp"

namespace anpi
{

template <typename T>
inline void lu(const anpi::Matrix<T> &A, anpi::Matrix<T> &LU, std::vector<size_t> &p)
{

    try
    {
        anpi::luCrout(A, LU, p);
    }
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

    anpi::Matrix<T> LU;
    std::vector<size_t> p;
    anpi::lu(A, LU, p);

    auto n = A.cols();

    int i, ii = 0, ip, j;
    T sum;
    if (b.size() != n || x.size() != n)
        throw anpi::Exception("solveLU::solve vector bad sizes");
    for (i = 0; i < n; ++i)
        x[i] = b[i];
    for (i = 0; i < n; ++i)
    {              //When ii is set to a positive value, it will become the index of the first nonvanishing element of b. We now
        ip = p[i]; //do the forward substitution, equation (2.3.6). The only new wrinkle is to unscramble the permutation as we go.
        sum = x[ip];
        x[ip] = x[i];
        if (ii != 0)
            for (j = ii - 1; j < i; ++j)
                sum -= LU[i][j] * x[j];
        else if (sum != 0.0) //A nonzero element was encountered, so from now on we will have to do the sums in the loop above.
            ii = i + 1;

        x[i] = sum;
    }
    for (i = n - 1; i >= 0; i--)
    { // Now we do the backsubstitution, equation (2.3.7).
        sum = x[i];
        for (j = i + 1; j < n; ++j)
            sum -= LU[i][j] * x[j];
        x[i] = sum / LU[i][i]; //Store a component of the solution vector X .
    }
}

} // namespace anpi