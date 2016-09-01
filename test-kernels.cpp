// g++ -Wall -O3 -march=native -o test-kernels test-kernels.cpp

#include <vector>
#include <cstdlib>
#include <cassert>
#include <iostream>

#include "ch_vdif_assembler_kernels.hpp"

using namespace std;
using namespace ch_vdif_assembler;


static uint8_t rand8()
{
    double r = (rand() + 0.5) / (RAND_MAX + 1.0);
    return (uint8_t)(256*r);
}

static vector<uint8_t> randvec(int nelts)
{
    assert(nelts > 0);
    vector<uint8_t> ret(nelts);

    for (int i = 0; i < nelts; i++)
	ret[i] = rand8();

    return ret;
}


// -------------------------------------------------------------------------------------------------


static void test_assemble128(int stride, int n, int dst_offset, int src_offset)
{
    assert(stride >= n);

    int nsrc = src_offset + 8*n + 1;             // the "+1" is for testing the n==0 case
    int ndst = dst_offset + 7*stride + n + 64;   // the "+64" is to test for writing past end of buffer

    vector<uint8_t> dst1 = randvec(ndst);
    vector<uint8_t> src1 = randvec(nsrc);

    vector<uint8_t> dst0 = dst1;
    vector<uint8_t> dst2 = dst1;
    vector<uint8_t> src2 = src1;

    _assemble8(&dst1[dst_offset], stride, &src1[src_offset], n);
    _assemble128(&dst2[dst_offset], stride, &src2[src_offset], n);

    for (int i = 0; i < ndst; i++) {
	if (dst1[i] == dst2[i])
	    continue;

	cerr << "\ntest_assemble128() failed\n"
	     << "  stride=" << stride << ", n=" << n << ", dst_offset=" << dst_offset << ", src_offset=" << src_offset << "\n";
	    
	for (int ii = max(i-8,0); ii < min(i+9,ndst); ii++) {
	    cerr << "   i=" << ii
		 << ", dst0=" << hex << (unsigned int)(dst0[ii]) << dec
		 << ", dst1=" << hex << (unsigned int)(dst1[ii]) << dec
		 << ", dst2=" << hex << (unsigned int)(dst2[ii]) << dec;

	    if (ii == i)
		cerr << "   [ first failure here ]";

	    cerr << "\n";
	}

	exit(1);
    }
}

static void test_assemble128()
{
    const int stride = 512;          // OK to fix this to a large power of 2
    const int nmax = 200;            // Go to a fairly large value, in anticipation of 64-byte AVX512 kernels
    const int dst_offset_max = 64;   // also in anticipation of 64-byte kernels
    const int src_offset_max = 64;   // also in anticipation of 64-byte kernels    
    
    cerr << "test_assemble128";

    for (int dst_offset = 0; dst_offset < dst_offset_max; dst_offset++) {
	cerr << ".";
	for (int src_offset = 0; src_offset < src_offset_max; src_offset++)
	    for (int n = 0; n < nmax; n++)
		test_assemble128(stride, n, dst_offset, src_offset);
    }

    cerr << "pass\n";
}


// -------------------------------------------------------------------------------------------------


static void test_sum16_auto_correlations()
{
    for (int n = 0; n < 1000000; n++) {
	vector<uint8_t> data = randvec(16);
	vector<uint8_t> data2 = data;
	
	int sum1, count1;
	_sum16_auto_correlations_reference(sum1, count1, &data[0]);
	
	int sum2, count2;
	_sum16_auto_correlations(sum2, count2, &data2[0]);
	
	if ((sum1 == sum2) && (count1 == count2))
	    continue;

	cerr << "test_sum16_auto_correlations() failed\n";
	cerr << "    data = [";

	for (int i = 0; i < 16; i++)
	    cerr << " " << hex << (unsigned int)(data[i]) << dec;

	cerr << " ]\n"
	     << "    sum1=" << sum1 << " count1=" << count1 << "\n"
	     << "    sum2=" << sum2 << " count2=" << count2 << "\n";

	exit(1);
    }

    cerr << "test_sum16_auto_correlations: pass\n";
}


// -------------------------------------------------------------------------------------------------


static void test_sum16_liam_hack()
{
    for (int n = 0; n < 1000000; n++) {
	vector<uint8_t> data = randvec(16);
	vector<uint8_t> data2 = data;
	
	int sum1, count1;
	_sum16_liam_hack_reference(sum1, count1, &data[0]);
	
	int sum2, count2;
	_sum16_liam_hack(sum2, count2, &data2[0]);
	
	if ((sum1 == sum2) && (count1 == count2))
	    continue;

	cerr << "test_sum16_liam_hack() failed\n";
	cerr << "    data = [";

	for (int i = 0; i < 16; i++)
	    cerr << " " << hex << (unsigned int)(data[i]) << dec;

	cerr << " ]\n"
	     << "    sum1=" << sum1 << " count1=" << count1 << "\n"
	     << "    sum2=" << sum2 << " count2=" << count2 << "\n";

	exit(1);
    }

    cerr << "test_sum16_liam_hack: pass\n";
}


// -------------------------------------------------------------------------------------------------


int main(int argc, char **argv)
{
    test_sum16_auto_correlations();
    test_sum16_liam_hack();
    test_assemble128();

    return 0;
}
