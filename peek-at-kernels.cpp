// only _assembler_kernel() for now

#include <cassert>
#include <vector>
#include <iostream>
#include <stdexcept>
#include "ch_vdif_assembler_kernels.hpp"

using namespace std;
using namespace ch_vdif_assembler;


static void peek_at_assembler_kernel()
{
    // _assembler_kernel() assumes unaligned input pointer
    vector<uint8_t> v(256, 0);
    for (int i = 0; i < 256; i++)
	v[i] = i;
    
    __m128i *src = reinterpret_cast<__m128i *> (&v[0]);

    __m128i x0, x1, x2, x3, x4, x5, x6, x7;
    _assembler_kernel(x0, x1, x2, x3, x4, x5, x6, x7, src);

    cout << "_assembler_kernel\n";
    cout << "   input\n";

    for (int i = 0; i < 8; i++) {
	cout << "    ";
	for (int j = 16*i; j < 16*i+16; j++)
	    cout << " " << hex << (unsigned int)(v[j]) << dec;
	cout << "\n";
    }

    cout << "   output\n";
    cout << "       " << str8(x0) << endl;
    cout << "       " << str8(x1) << endl;
    cout << "       " << str8(x2) << endl;
    cout << "       " << str8(x3) << endl;
    cout << "       " << str8(x4) << endl;
    cout << "       " << str8(x5) << endl;
    cout << "       " << str8(x6) << endl;
    cout << "       " << str8(x7) << endl;
}


int main(int argc, char **argv)
{
    peek_at_assembler_kernel();
    return 0;
}
