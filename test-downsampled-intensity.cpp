#include "ch_vdif_assembler_internals.hpp"

using namespace std;
using namespace ch_vdif_assembler;


int main(int argc, char **argv)
{
    cerr << "test-downsampled_intensity (this test takes a few minutes)";

    for (int iouter = 0; iouter < 100; iouter++) {
	cerr << ".";

	int nt_downsample = pow2(randint(4,7));
	int nt_assembler = pow2(randint(nt_downsample+3, 11));
	int64_t current_t0 = nt_assembler * randint(0, 10);
	int nchunks = randint(0, 100);

	shared_ptr<assembled_chunk_pool> pool = make_shared<assembled_chunk_pool> (nt_assembler);
	downsampled_intensity ds1(nt_downsample);
	downsampled_intensity ds2(nt_downsample);

	for (int ichunk = 0; ichunk < nchunks; ichunk++) {
	    shared_ptr<assembled_chunk> chunk = assembled_chunk::make_random(pool, current_t0);
	    current_t0 = chunk->t0 + chunk->nt;

	    ds1.process_chunk(chunk);
	    ds2.process_chunk_reference(chunk);

	    xassert(ds1.initial_t0 == ds2.initial_t0);
	    xassert(ds1.curr_chunk_t0 == ds2.curr_chunk_t0);
	    xassert(ds1.curr_chunk_nt == ds2.curr_chunk_nt);

	    int nbuf = constants::chime_nfreq * 2 * (nt_assembler/nt_downsample);

	    for (int i = 0; i < nbuf; i++) {
		xassert(fabs(ds1.intensity_buf[i] - ds2.intensity_buf[i]) < 1.0e-4);
		xassert(fabs(ds1.weights_buf[i] - ds2.weights_buf[i]) < 1.0e-4);
	    }
	}
    }

    cerr << "pass\n";
    return 0;
}
