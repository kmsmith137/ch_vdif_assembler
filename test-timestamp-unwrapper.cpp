#include <cassert>
#include "ch_vdif_assembler_internals.hpp"

using namespace std;
using namespace ch_vdif_assembler;


static void test_initialization(uint32_t fpga0)
{
    timestamp_unwrapper tsu;
    int64_t t0 = tsu.unwrap(fpga0);

    assert(t0 >= 0);
    assert(t0 < (int64_t(1) << 32));
    assert((uint32_t)t0 == fpga0);
}

int main(int argc, char **argv)
{
    srand(time(NULL));

    // If this fails, then the test needs to be rethought....
    assert(RAND_MAX == 0x7fffffff);

    test_initialization(0);
    test_initialization(1);
    test_initialization(1 << 10);
    test_initialization(1 << 20);
    test_initialization(1 << 31);
    test_initialization((uint32_t)(-1));
    test_initialization((uint32_t)(-2));

    for (int n = 0; n < 100; n++) {
	timestamp_unwrapper tsu;

	uint32_t fpga = (uint32_t)rand() + (uint32_t)rand();
	int64_t expected_ts = tsu.unwrap(fpga);

	for (int i = 0; i < 100000; i++) {
	    // this is the shift we'll apply to the FPGA count
	    int r = rand() - RAND_MAX/2;

	    //
	    // The following ugliness is morally equivalent to
	    //    fpga += r
	    // but makes no assumptions about how signs or overflows
	    // are handled in casts.  (To avoid a circular test, since
	    // timestamp_unwrapper::unwrap() does make assumptions!)
	    //
	    uint32_t ru = (r >= 0) ? (uint32_t)r : (~(uint32_t)(-r) + 1);
	    fpga = (uint32_t)(((uint64_t)fpga + ru) & 0xffffffffULL);

	    int64_t ts = tsu.unwrap(fpga);
	    expected_ts += (int64_t)r;

	    if (ts != expected_ts) {
		cerr << "timestamp_unwrapper test: fail\n";
		exit(1);
	    }
	}
    }

    cerr << "test_timestamp_unwrapper: pass\n";
    return 0;
}
