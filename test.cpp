//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <iostream>

#include "palisade_header.h"

#include <getopt.h>

using namespace std;

int main()
{
	int64_t x1=1002345;
	int64_t x2=99878;

	CT c1=encrypt_plaintext_integer_to_ciphertext(x1);
	CT c2=encrypt_plaintext_integer_to_ciphertext(x2);

	CT c_res=cc->EvalAdd(c1,c2);

	cout<<"cres: "<<decrypt_ciphertext_to_plaintext_vector(c_res)<<endl;

	return 0;
}