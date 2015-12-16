// only _assembler_kernel() timed for now

#include <sys/time.h>
#include <iostream>
#include <vector>
#include <stdexcept>
#include "ch_vdif_assembler_kernels.hpp"

using namespace std;
using namespace ch_vdif_assembler;


inline double time_diff(const struct timeval &tv1, const struct timeval &tv2)
{
    return (tv2.tv_sec - tv1.tv_sec) + 1.0e-6 * (tv2.tv_usec - tv1.tv_usec);
}

inline struct timeval get_time()
{
    struct timeval ret;
    if (gettimeofday(&ret, NULL) < 0)
	throw std::runtime_error("gettimeofday() failed");
    return ret;
}


static void time_assembler_kernel()
{
    // _assembler_kernel() assumes unaligned input pointer
    vector<uint8_t> v(1024, 0);
    __m128i *src = reinterpret_cast<__m128i *> (&v[0]);
    
    __m128i x0, x1, x2, x3, x4, x5, x6, x7;
    int n = (1 << 28);

    struct timeval tv0 = get_time();

    for (int i = 0; i < n; i++)
	_assembler_kernel(x0, x1, x2, x3, x4, x5, x6, x7, src);

    double dt = time_diff(tv0, get_time());

    double kernels_per_sec = (double)n / dt;
    cout << "_assembler_kernel: " << kernels_per_sec << " kernels/sec\n";
}


int main(int argc, char **argv)
{
    time_assembler_kernel();
    return 0;
}
