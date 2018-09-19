#include <boost/test/unit_test.hpp>

// #include <iostream>
// #include <exception>
// #include <cstdlib>
// #include <complex>

/**
 * Unit tests for the matrix class
 */
#include "benchmarkFramework.hpp"
// #include <Matrix.hpp>
// #include "Allocator.hpp"
#include "LUCrout.hpp"
#include "LUDoolittle.hpp"

BOOST_AUTO_TEST_SUITE(LU)

/// Benchmark for LU decomposition
template <typename T>
class benchLU
{
  protected:
    /// Maximum allowed size for the square matrices
    const size_t _maxSize;

    /// A large matrix holding
    anpi::Matrix<T> _data;

    /// State of the benchmarked evaluation
    anpi::Matrix<T> _a;
    anpi::Matrix<T> _b;
    std::vector<size_t> _c;

  public:
    /// Construct
    benchLU(const size_t maxSize)
        : _maxSize(maxSize), _data(maxSize, maxSize, anpi::DoNotInitialize)
    {

        size_t idx = 0;
        for (size_t r = 0; r < _maxSize; ++r)
        {
            for (size_t c = 0; c < _maxSize; ++c)
            {
                _data(r, c) = idx++;
            }
        }
    }

    /// Prepare the evaluation of given size
    void prepare(const size_t size)
    {
        assert(size <= this->_maxSize);
        this->_a = std::move(anpi::Matrix<T>(size, size, _data.data()));
        this->_b = this->_a;
    }
};

/// Provide the evaluation method for Crout LU decomposition
template <typename T>
class benchCrout : public benchLU<T>
{
  public:
    /// Constructor
    benchCrout(const size_t n) : benchLU<T>(n) {}

    // Evaluate add in-place
    inline void eval()
    {
        anpi::luCrout(this->_a, this->_b, this->_c);
    }
};

/// Provide the evaluation method for Doolittle LU decomposition
template <typename T>
class benchDoolittle : public benchLU<T>
{
  public:
    /// Constructor
    benchDoolittle(const size_t n) : benchLU<T>(n) {}

    // Evaluate add in-place
    inline void eval()
    {
        anpi::luDoolittle(this->_a, this->_b, this->_c);
    }
};

/**
 * Instantiate and test the methods of the Matrix class
 */
BOOST_AUTO_TEST_CASE(LUbench)
{

    std::vector<size_t> sizes = {24, 32, 48, 64,
                                 96, 128, 192, 256,
                                 384, 512};

    const size_t n = sizes.back();
    const size_t repetitions = 100;
    std::vector<anpi::benchmark::measurement> times;

    { //benchmark for Crout
        benchCrout<float> bc(n);

        // Measure on-copy add
        ANPI_BENCHMARK(sizes, repetitions, times, bc);

        anpi::benchmark::write("LU_Crout_float_bench.txt", times);
        anpi::benchmark::plotRange(times, "Crout (float)", "r");
    }

    { //benchmark for Doolittle
        benchDoolittle<float> bdl(n);

        // Measure on-copy add
        ANPI_BENCHMARK(sizes, repetitions, times, bdl);

        anpi::benchmark::write("LU_Doolittle_float_bench.txt", times);
        anpi::benchmark::plotRange(times, "Doolittle (float)", "g");
    }

    anpi::benchmark::show();
}

BOOST_AUTO_TEST_SUITE_END()