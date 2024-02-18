/*************************************************************************************
						   The MIT License

   BWA-MEM2  (Sequence alignment using Burrows-Wheeler Transform),
   Copyright (C) 2019  Vasimuddin Md, Sanchit Misra, Intel Corporation, Heng Li.

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.

Authors: Vasimuddin Md <vasimuddin.md@intel.com>; Sanchit Misra <sanchit.misra@intel.com>.
*****************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <limits>
#include <unistd.h>
#include <omp.h>
#include <fstream>
#include "utils.h"
#include "bandedSWA.h"
#include "palisade_header.h"

#include <chrono>
using namespace std::chrono;

// #define VTUNE_ANALYSIS 1

#define CLMUL 8
#ifdef VTUNE_ANALYSIS
#include <ittnotify.h>
#endif

#define DEFAULT_MATCH 1
#define DEFAULT_MISMATCH 4
#define DEFAULT_OPEN 6
#define DEFAULT_EXTEND 1
#define DEFAULT_AMBIG -1

#undef MAX_SEQ_LEN_REF
#define MAX_SEQ_LEN_REF 2048
#undef MAX_SEQ_LEN_QER
#define MAX_SEQ_LEN_QER 256

// #define MAX_NUM_PAIRS 1000
// #define MATRIX_MIN_CUTOFF -100000000
// #define LOW_INIT_VALUE (INT32_MIN/2)
// #define AMBIG 52
double freq = 2.6 * 1e9;
int32_t w_match, w_mismatch, w_open, w_extend, w_ambig, numThreads = 1, batchSize = 0;
uint64_t SW_cells;
char *pairFileName;
FILE *pairFile;
int8_t h0 = 0;
double clock_freq;
uint64_t prof[10][112], data, SW_cells2;

void bwa_fill_scmat(int a, int b, int ambig, int8_t mat[25])
{
	int i, j, k;
	for (i = k = 0; i < 4; ++i)
	{
		for (j = 0; j < 4; ++j)
			mat[k++] = i == j ? a : -b;
		mat[k++] = ambig; // ambiguous base
	}
	for (j = 0; j < 5; ++j)
		mat[k++] = ambig;
}

void parseCmdLine(int argc, char *argv[])
{
	int i;
	w_match = DEFAULT_MATCH;
	w_mismatch = DEFAULT_MISMATCH;
	w_open = DEFAULT_OPEN;
	w_extend = DEFAULT_EXTEND;
	w_ambig = DEFAULT_AMBIG;

	int pairFlag = 0;
	for (i = 1; i < argc; i += 2)
	{
		if (strcmp(argv[i], "-match") == 0)
		{
			w_match = atoi(argv[i + 1]);
		}
		if (strcmp(argv[i], "-mismatch") == 0)
		{ // penalty, +ve number
			w_mismatch = atoi(argv[i + 1]);
		}
		if (strcmp(argv[i], "-ambig") == 0)
		{
			w_ambig = atoi(argv[i + 1]);
		}

		if (strcmp(argv[i], "-gapo") == 0)
		{
			w_open = atoi(argv[i + 1]);
		}
		if (strcmp(argv[i], "-gape") == 0)
		{
			w_extend = atoi(argv[i + 1]);
		}
		if (strcmp(argv[i], "-pairs") == 0)
		{
			pairFileName = argv[i + 1];
			pairFlag = 1;
		}
		if (strcmp(argv[i], "-h0") == 0)
		{
			h0 = atoi(argv[i + 1]);
		}
		if (strcmp(argv[i], "-t") == 0)
		{
			numThreads = atoi(argv[i + 1]);
		}
		if (strcmp(argv[i], "-b") == 0)
		{
			batchSize = atoi(argv[i + 1]);
		}
	}
	if (pairFlag == 0)
	{
		fprintf(stderr, "ERROR! pairFileName not specified.\n");
		exit(EXIT_FAILURE);
	}
}

// -------------------------------------------------------------------------
// INPUT FILE FORMAT
// -------------------------------------------------------------------------
// Line 1: Seed score
// Line 2: Reference string
// Line 3: Query string
//
// E.g:
// 19
// 01230123
// 0123

void loadPairs(SeqPair *seqPairArray, uint8_t *seqBufRef, uint8_t* seqBufQer, size_t numPairs)
{
	size_t numPairsRead = 0;
	while (numPairsRead < numPairs) {
		int32_t h0 = 0;
		char temp[10];
		fgets(temp, 10, pairFile);
		sscanf(temp, "%d", &h0);
		if (!fgets((char *)(seqBufRef + numPairsRead * (int64_t)(MAX_SEQ_LEN_REF)), MAX_SEQ_LEN_REF, pairFile)) {
			printf("WARNING! fgets returned NULL in %s. Num Pairs : %d\n", pairFileName, numPairsRead);
			break;
        }
		if (!fgets((char *)(seqBufQer + numPairsRead * (int64_t)(MAX_SEQ_LEN_QER)), MAX_SEQ_LEN_QER, pairFile)) {
			printf("WARNING! Odd number of sequences in %s\n", pairFileName);
			break;
        }
		SeqPair sp;
		sp.id = numPairsRead;
		sp.len1 = strnlen((char *)(seqBufRef + numPairsRead * MAX_SEQ_LEN_REF), MAX_SEQ_LEN_REF) - 1;
		sp.len2 = strnlen((char *)(seqBufQer + numPairsRead * MAX_SEQ_LEN_QER), MAX_SEQ_LEN_QER) - 1;
        if (sp.len1 <= 0 || sp.len2 <= 0) {
            fprintf(stderr, "%d\n", numPairsRead);
        }
        assert(sp.len1 > 0);
        assert(sp.len2 > 0);
		sp.h0 = h0;
		uint8_t *seq1 = seqBufRef + numPairsRead * MAX_SEQ_LEN_REF;
		uint8_t *seq2 = seqBufQer + numPairsRead * MAX_SEQ_LEN_QER;
		sp.idr =  numPairsRead * MAX_SEQ_LEN_REF;
		sp.idq =  numPairsRead * MAX_SEQ_LEN_QER;
		for (int l = 0; l < sp.len1; l++) {
			seq1[l] -= 48;
        }
		for (int l = 0; l < sp.len2; l++) {
			seq2[l] -= 48;
        }
		sp.seqid = sp.regid = sp.score = sp.tle = sp.gtle = sp.qle = -1;
		sp.gscore = sp.max_off = -1;
		seqPairArray[numPairsRead] = sp;
		numPairsRead++;
		// SW_cells += (sp.len1 * sp.len2);
	}
}

// new seqBufRef is passed for every call
void loadPairs(SeqPair *seqPairArray, uint8_t *seqBufRef, CT *seqBufRef_ct, int i1, uint8_t *seqBufQer, CT *seqBufQer_ct, int i2, size_t numPairs)
{
	// cout << "i1: " << i1 << endl;

	size_t numPairsRead = 0;
	while (numPairsRead < numPairs)
	{
		// cout << "numPairsRead: " << numPairsRead << endl;
		// goes from 0 to 511
		// Error at numPairsRead: 152 for i=99840

		int32_t h0 = 0;
		char temp[10];

		// 		char *fgets(char *str, int n, FILE *stream);
		// 		str is the buffer where the read string will be stored.
		// 		n is the maximum number of characters to be read, including the terminating null character.
		// 		stream is a pointer to the input stream.
		//		The function returns a pointer to the str buffer if successful, and NULL otherwise

		fgets(temp, 10, pairFile);
		// printf("temp: %s\n",temp);

		// int sscanf(const char *str, const char *format, ...);
		// str is the input string containing the data to be read.
		// format is a string that specifies the format of the data to be read, using format specifiers like %d, %f, %s, etc.
		// ... is a variable-length argument list that contains the addresses of the variables where the data will be stored.
		// The number of arguments after format must match the number of format specifiers in format.

		sscanf(temp, "%d", &h0);
		// printf("h0: %d\n",h0);

		// SeqPair *seqPairArray = (SeqPair *)_mm_malloc(roundNumPairs * sizeof(SeqPair), 64);				   // roundNumPairs: 100000
		// uint8_t *seqBufQer = (uint8_t*) _mm_malloc(MAX_SEQ_LEN_QER * roundNumPairs * sizeof(int8_t), 64);   // MAX_SEQ_LEN_QER=256
		// uint8_t *seqBufRef = (uint8_t*) _mm_malloc(MAX_SEQ_LEN_REF * roundNumPairs * sizeof(int8_t), 64);   // MAX_SEQ_LEN_REF=2048

		// fgets also reads 0x00 (or NULL) character
		if (!fgets((char *)(seqBufRef + numPairsRead * (int64_t)(MAX_SEQ_LEN_REF)), MAX_SEQ_LEN_REF, pairFile))
		{
			printf("WARNING! fgets returned NULL in %s. Num Pairs : %d\n", pairFileName, numPairsRead);
			break;
		}
		if (!fgets((char *)(seqBufQer + numPairsRead * (int64_t)(MAX_SEQ_LEN_QER)), MAX_SEQ_LEN_QER, pairFile))
		{
			printf("WARNING! Odd number of sequences in %s\n", pairFileName);
			break;
		}

		// for (int ii = 0; ii < MAX_SEQ_LEN_REF; ii++)
		// {
		// 	cout << (seqBufRef + numPairsRead * (int64_t)(MAX_SEQ_LEN_REF))[ii] << " ";
		// 	if ((seqBufRef + numPairsRead * (int64_t)(MAX_SEQ_LEN_REF))[ii] == '\n')
		// 		break;
		// }
		// // cout << endl;
		// cout << "------------------------------------------" << endl;

		// int seq_buf_ref_read_count = strlen((char *)(seqBufRef + numPairsRead * (int64_t)(MAX_SEQ_LEN_REF))); // strlen will include the final '\n' that was read by fgets (it doesnt read '\0')
		// int seq_buf_qer_read_count = strlen((char *)(seqBufQer + numPairsRead * (int64_t)(MAX_SEQ_LEN_QER)));
		// cout << "seq_buf_ref_read_count: " << seq_buf_ref_read_count << endl;
		// cout << "seq_buf_qer_read_count: " << seq_buf_qer_read_count << endl;

		// So, given any index in seqBufRef (i1 + numPairsRead * (int64_t)(MAX_SEQ_LEN_REF) in this case), we can find which
		// ciphertext ct_cnt_1 it belongs to and at what index ind1 in its vector, we need to start reading its string(till we hit '\n')
		int tot_1 = 0;
		int ct_cnt_1 = (i1 + numPairsRead * (int64_t)(MAX_SEQ_LEN_REF)) / 16384;
		int ct_cnt_1_copy = ct_cnt_1;
		int ind1 = (i1 + numPairsRead * (int64_t)(MAX_SEQ_LEN_REF)) % 16384; // As numPairsRead goes from 0 to numPairs, incrementing by 1, so here, ind1 increases by 2048 each time ind1 increases by 1
		int ind1_cpy = ind1;

		vecInt v_ref = decrypt_ciphertext_to_plaintext_vector(*(seqBufRef_ct + ct_cnt_1));

		// while (tot_1 < seq_buf_ref_read_count)
		while((seqBufRef + numPairsRead * (int64_t)(MAX_SEQ_LEN_REF))[tot_1])
		{
			if (ind1 == 16384)
			{
				*(seqBufRef_ct + ct_cnt_1) = encrypt_plaintext_vector_to_ciphertext(v_ref);
				ct_cnt_1++;
				v_ref = decrypt_ciphertext_to_plaintext_vector(*(seqBufRef_ct + ct_cnt_1));
				ind1 = 0;
			}

			v_ref[ind1] = (seqBufRef + numPairsRead * (int64_t)(MAX_SEQ_LEN_REF))[tot_1];
			ind1++;
			tot_1++;
		}
		*(seqBufRef_ct + ct_cnt_1) = encrypt_plaintext_vector_to_ciphertext(v_ref);
		ct_cnt_1++;
		ind1 = 0;

		int second_cnt = 0;
		// for (int ii = ct_cnt_1_copy; ii < ct_cnt_1 && second_cnt < seq_buf_ref_read_count; ii++)
		// {
		// 	vecInt v_ref = decrypt_ciphertext_to_plaintext_vector(*(seqBufRef_ct + ii));
		// 	// cout<<ii<<" ";
		// 	if (ii == ct_cnt_1_copy)
		// 	{
		// 		for (int jj = ind1_cpy; jj < v_ref.size(); jj++)
		// 		{
		// 			cout << (char)(v_ref[jj]) << " ";
		// 			second_cnt++;
		// 			if (second_cnt == seq_buf_ref_read_count)
		// 				break;
		// 		}
		// 		// second_cnt+=v_ref.size()-ind1_cpy;
		// 	}
		// 	else
		// 	{
		// 		for (int jj = 0; jj < v_ref.size(); jj++)
		// 		{
		// 			cout << (char)(v_ref[jj]) << " ";
		// 			second_cnt++;
		// 			if (second_cnt == seq_buf_ref_read_count)
		// 				break;
		// 		}
		// 		// second_cnt+=v_ref.size();
		// 	}
		// }
		// cout << endl;
		// cout << "###############################################" << endl;

		// for (int ii = 0; ii < MAX_SEQ_LEN_QER; ii++)
		// {
		// 	cout << (int)(seqBufQer + numPairsRead * (int64_t)(MAX_SEQ_LEN_QER))[ii] << " ";
		// 	if ((seqBufQer + numPairsRead * (int64_t)(MAX_SEQ_LEN_QER))[ii] == '\n')
		// 		break;
		// }
		// cout << endl;
		// cout << "------------------------------------------" << endl;

		// So, given any index in seqBufRef (i1 + numPairsRead * (int64_t)(MAX_SEQ_LEN_REF) in this case), we can find which
		// ciphertext ct_cnt_1 it belongs to and at what index ind1 in its vector, we need to start reading its string(till we hit '\n')
		int tot_2 = 0;
		int ct_cnt_2 = (i2 + numPairsRead * (int64_t)(MAX_SEQ_LEN_QER)) / 16384;
		int ct_cnt_2_copy = ct_cnt_2;
		int ind2 = (i2 + numPairsRead * (int64_t)(MAX_SEQ_LEN_QER)) % 16384; // As numPairsRead goes from 0 to numPairs, incrementing by 1, so here, ind1 increases by 2048 each time ind1 increases by 1
		int ind2_cpy = ind2;

		v_ref = decrypt_ciphertext_to_plaintext_vector(*(seqBufQer_ct + ct_cnt_2));

		// while (tot_2 < seq_buf_qer_read_count)
		while((seqBufQer + numPairsRead * (int64_t)(MAX_SEQ_LEN_QER))[tot_2])
		{
			if (ind2 == 16384)
			{
				*(seqBufQer_ct + ct_cnt_2) = encrypt_plaintext_vector_to_ciphertext(v_ref);
				ct_cnt_2++;
				v_ref = decrypt_ciphertext_to_plaintext_vector(*(seqBufRef_ct + ct_cnt_2));
				ind2 = 0;
			}

			v_ref[ind2] = (seqBufQer + numPairsRead * (int64_t)(MAX_SEQ_LEN_QER))[tot_2];
			ind2++;
			tot_2++;
		}
		*(seqBufQer_ct + ct_cnt_2) = encrypt_plaintext_vector_to_ciphertext(v_ref);
		ct_cnt_2++;
		ind2 = 0;

		int seq_buf_ref_read_count; int seq_buf_qer_read_count;

		int cntref= (i1 + numPairsRead * MAX_SEQ_LEN_REF) / 16384;
		int indref = (i1 + numPairsRead * MAX_SEQ_LEN_REF) % 16384;

		vecInt vv12(16384, 0);
		for (int ii = indref; ii < 16384 && ii < indref + MAX_SEQ_LEN_REF; ii++)
		{
			vv12[ii] = (int64_t)('\n');
		}
		CT c_subb = cc->EvalSub(*(seqBufRef_ct + cntref), encode_vector_to_plaintext(vv12));
		vecInt v_subb = decrypt_ciphertext_to_plaintext_vector(c_subb);
		// cout<<"Dec chars: "<<endl;
		for (int ii = indref; ii < 16384 && ii < indref + MAX_SEQ_LEN_REF; ii++)
		{
			if (v_subb[ii] == 0)
			{
				seq_buf_ref_read_count = ii - indref+1;
				break;
			}
		}

		cntref= (i2 + numPairsRead * MAX_SEQ_LEN_QER) / 16384;
		indref = (i2 + numPairsRead * MAX_SEQ_LEN_QER) % 16384;

		vv12.resize(16384,0);
		for (int ii = indref; ii < 16384 && ii < indref + MAX_SEQ_LEN_QER; ii++)
		{
			vv12[ii] = (int64_t)('\n');
		}
		c_subb = cc->EvalSub(*(seqBufQer_ct + cntref), encode_vector_to_plaintext(vv12));
		v_subb = decrypt_ciphertext_to_plaintext_vector(c_subb);
		// cout<<"Dec chars: "<<endl;
		for (int ii = indref; ii < 16384 && ii < indref + MAX_SEQ_LEN_QER; ii++)
		{
			if (v_subb[ii] == 0)
			{
				seq_buf_qer_read_count = ii - indref+1;
				break;
			}
		}

		// cout << "seq_buf_ref_read_count: " << seq_buf_ref_read_count << endl;
		// cout << "seq_buf_qer_read_count: " << seq_buf_qer_read_count << endl<<endl;

		second_cnt = 0;
		// // seq_buf_qer_read_count includes the final \n
		// for (int ii = ct_cnt_2_copy; ii < ct_cnt_2 && second_cnt < seq_buf_qer_read_count; ii++)
		// {
		// 	vecInt v_ref = decrypt_ciphertext_to_plaintext_vector(*(seqBufQer_ct + ii));
		// 	// cout<<ii<<" ";
		// 	if (ii == ct_cnt_2_copy)
		// 	{
		// 		for (int jj = ind2_cpy; jj < v_ref.size(); jj++)
		// 		{
		// 			cout <<(int)(v_ref[jj]);
		// 			second_cnt++;
		// 			if (second_cnt == seq_buf_qer_read_count)
		// 				break;
		// 			cout<<" ";
		// 		}
		// 		// second_cnt+=v_ref.size()-ind1_cpy;
		// 	}
		// 	else
		// 	{
		// 		for (int jj = 0; jj < v_ref.size(); jj++)
		// 		{
		// 			cout <<(int)(v_ref[jj]);
		// 			second_cnt++;
		// 			if (second_cnt == seq_buf_qer_read_count)
		// 				break;
		// 			cout<<" ";
		// 		}
		// 		// second_cnt+=v_ref.size();
		// 	}
		// }
		// cout<<endl;

		// typedef struct dnaSeqPair
		// {
		// 	int64_t idr, idq, id;
		// 	int32_t len1, len2;
		// 	int32_t h0;
		// 	int seqid, regid;
		// 	int32_t score, tle, gtle, qle;
		// 	int32_t gscore, max_off;
		// }SeqPair;

		SeqPair sp;
		sp.id = numPairsRead;
		// // strnlen gives the length of null-terminated string; It also stops when it hits \n
		// sp.len1 = strnlen((char *)(seqBufRef + numPairsRead * MAX_SEQ_LEN_REF), MAX_SEQ_LEN_REF) - 1;
		// sp.len2 = strnlen((char *)(seqBufQer + numPairsRead * MAX_SEQ_LEN_QER), MAX_SEQ_LEN_QER) - 1;

		// cout<<"sp.len1: "<<sp.len1<<"; sp.len2: "<<sp.len2<<endl;
		// cout<<"Original chars"<<endl;
		// for(int ii=0;ii<MAX_SEQ_LEN_REF && *(seqBufRef + numPairsRead * MAX_SEQ_LEN_REF+ii)!='\n';ii++)
		// {
		// 	printf("%c ",*(seqBufRef + numPairsRead * MAX_SEQ_LEN_REF+ii));
		// }
		// cout<<endl;

		ct_cnt_1 = (i1 + numPairsRead * MAX_SEQ_LEN_REF) / 16384;
		ind1 = (i1 + numPairsRead * MAX_SEQ_LEN_REF) % 16384;

		// CT c_orig=*(seqBufRef_ct + ct_cnt_1);
		// vecInt v_origg=decrypt_ciphertext_to_plaintext_vector(c_orig);
		// for (int ii = ind1; ii < 16384 && v_origg[ii]!='\n' ; ii++)
		// {
		// 	printf("%c ",v_origg[ii]);
		// }
		// cout<<endl;

		vecInt vv1(16384, 0);
		for (int ii = ind1; ii < 16384 && ii < ind1 + MAX_SEQ_LEN_REF; ii++)
		{
			vv1[ii] = (int64_t)('\n');
		}
		CT c_sub = cc->EvalSub(*(seqBufRef_ct + ct_cnt_1), encode_vector_to_plaintext(vv1));
		vecInt v_sub = decrypt_ciphertext_to_plaintext_vector(c_sub);
		// cout<<"Dec chars: "<<endl;
		for (int ii = ind1; ii < 16384 && ii < ind1 + MAX_SEQ_LEN_REF; ii++)
		{
			// cout<<(char)(v_sub[ii])<<" ";
			if (v_sub[ii] == 0)
			{
				sp.len1 = ii - ind1;
				break;
			}
		}

		// cout<<"Original chars qer"<<endl;
		// for(int ii=0;ii<MAX_SEQ_LEN_QER && *(seqBufQer + numPairsRead * MAX_SEQ_LEN_QER+ii)!='\n';ii++)
		// {
		// 	printf("%c ",*(seqBufQer + numPairsRead * MAX_SEQ_LEN_QER+ii));
		// }
		// cout<<"Qer_orig_done"<<endl;

		// ct_cnt_1 = (i2 + numPairsRead * MAX_SEQ_LEN_QER) / 16384;
		// ind1 = (i2 + numPairsRead * MAX_SEQ_LEN_QER) % 16384;

		// c_orig=*(seqBufQer_ct + ct_cnt_1);
		// v_origg=decrypt_ciphertext_to_plaintext_vector(c_orig);
		// for (int ii = ind1; ii < 16384 && v_origg[ii]!='\n' ; ii++)
		// {
		// 	printf("%c ",v_origg[ii]);
		// }
		// cout<<endl;

		ct_cnt_2 = (i2 + numPairsRead * MAX_SEQ_LEN_QER) / 16384;
		ind2 = (i2 + numPairsRead * MAX_SEQ_LEN_QER) % 16384;
		vecInt vv2(16384, 0);
		for (int ii = ind2; ii < 16384 && ii < ind2 + MAX_SEQ_LEN_QER; ii++)
		{
			vv2[ii] = (int64_t)('\n');
		}
		c_sub = cc->EvalSub(*(seqBufQer_ct + ct_cnt_2), encode_vector_to_plaintext(vv2));
		v_sub = decrypt_ciphertext_to_plaintext_vector(c_sub);
		for (int ii = ind2; ii < 16384 && ii < ind2 + MAX_SEQ_LEN_QER; ii++)
		{
			if (v_sub[ii] == 0)
			{
				sp.len2 = ii - ind2;
				break;
			}
		}
		// cout<<"new sp.len1: "<<sp.len1<<"; sp.len2: "<<sp.len2<<endl<<endl;

		if (sp.len1 <= 0 || sp.len2 <= 0)
		{
			fprintf(stderr, "%d\n", numPairsRead);
		}
		assert(sp.len1 > 0);
		assert(sp.len2 > 0);
		sp.h0 = h0;
		uint8_t *seq1 = seqBufRef + numPairsRead * MAX_SEQ_LEN_REF;
		uint8_t *seq2 = seqBufQer + numPairsRead * MAX_SEQ_LEN_QER;
		sp.idr = numPairsRead * MAX_SEQ_LEN_REF;
		sp.idq = numPairsRead * MAX_SEQ_LEN_QER;

		// max_seq1=max(max_seq1,sp.len1);
		// max_seq2=max(max_seq2,sp.len2);

		for (int l = 0; l < sp.len1; l++)
		{
			seq1[l] -= 48;
			// printf("%u ",seq1[l]);
		}
		// cout<<endl<<endl;

		int i_seq1 = i1 + numPairsRead * MAX_SEQ_LEN_REF;
		int ct_seq1 = i_seq1 / 16384;
		int ind_seq1 = i_seq1 % 16384;
		int ct_seq2 = -1;
		int ind_seq2 = -1;

		int num_of_ints_in_first_ct = 16384 - ind_seq1;
		int num_of_ints_in_second_ct = 0;
		if (num_of_ints_in_first_ct < sp.len1)
		{
			num_of_ints_in_second_ct = sp.len1 - num_of_ints_in_first_ct;
			ct_seq2 = ct_seq1 + 1;
			ind_seq2 = 0;
		}
		// CT seq1_ciphertext_1=*(seqBufRef_ct+ct_seq1);
		vecInt v_seq1(16384, 0);
		for (int kk = ind_seq1; kk < 16384 && kk < ind_seq1 + sp.len1; kk++)
			v_seq1[kk] = 48;
		*(seqBufRef_ct + ct_seq1) = cc->EvalSub(*(seqBufRef_ct + ct_seq1), encode_vector_to_plaintext(v_seq1));
		// vecInt v_remaining_1 = decrypt_ciphertext_to_plaintext_vector(*(seqBufRef_ct + ct_seq1));
		// for(int kk=ind_seq1;kk<16384 && kk<ind_seq1+sp.len1;kk++)
		// 	cout<<v_remaining_1[kk]<<" ";
		// cout<<endl;
		if (ct_seq2 != -1)
		{
			// cout<<"Second!"<<endl;
			vecInt v_seq2(16384, 0);
			for (int kk = 0; kk < num_of_ints_in_second_ct; kk++)
				v_seq2[kk] = 48;
			*(seqBufRef_ct + ct_seq2) = cc->EvalSub(*(seqBufRef_ct + ct_seq2), encode_vector_to_plaintext(v_seq2));
			// vecInt v_remaining_2 = decrypt_ciphertext_to_plaintext_vector(*(seqBufRef_ct + ct_seq2));
			// for(int kk=ind_seq2;kk<num_of_ints_in_second_ct;kk++)
			// 	cout<<v_remaining_2[kk]<<" ";
			// cout<<endl;
		}
		// cout<<endl;
		// cout<<"seq2: "<<endl;
		for (int l = 0; l < sp.len2; l++)
		{
			seq2[l] -= 48;
			// printf("%u ", seq2[l]);
		}
		// cout << endl
		// 	 << endl;

		int i_seq2 = i2 + numPairsRead * MAX_SEQ_LEN_QER;
		ct_seq1 = i_seq2 / 16384;
		ind_seq1 = i_seq2 % 16384;
		ct_seq2 = -1;
		ind_seq2 = -1;

		num_of_ints_in_first_ct = 16384 - ind_seq1;
		num_of_ints_in_second_ct = 0;
		if (num_of_ints_in_first_ct < sp.len2)
		{
			num_of_ints_in_second_ct = sp.len2 - num_of_ints_in_first_ct;
			ct_seq2 = ct_seq1 + 1;
			ind_seq2 = 0;
		}
		// CT seq2_ciphertext_1=*(seqBufQer_ct+ct_seq1);
		vecInt v_seq2(16384, 0);
		for (int kk = ind_seq1; kk < 16384 && kk < ind_seq1 + sp.len2; kk++)
			v_seq2[kk] = 48;
		// vecInt v_r_2 = decrypt_ciphertext_to_plaintext_vector(*(seqBufQer_ct + ct_seq1));
		// cout << "Earlier" << endl;
		// for (int kk = ind_seq1; kk < 16384 && kk < ind_seq1 + sp.len2; kk++)
		// 	cout <<(int)(v_r_2[kk])<<" ";
		// cout << endl;
		*(seqBufQer_ct + ct_seq1) = cc->EvalSub(*(seqBufQer_ct + ct_seq1), encode_vector_to_plaintext(v_seq2));
		// vecInt v_remaining_2 = decrypt_ciphertext_to_plaintext_vector(*(seqBufQer_ct + ct_seq1));
		// for (int kk = ind_seq1; kk < 16384 && kk < ind_seq1 + sp.len2; kk++)
		// 	cout << (int)(v_remaining_2[kk]) << " ";
		// cout << endl;

		if (ct_seq2 != -1)
		{
			// cout << "Second!" << endl;
			vecInt v_seq2(16384, 0);
			for (int kk = 0; kk < num_of_ints_in_second_ct; kk++)
				v_seq2[kk] = 48;
			*(seqBufQer_ct + ct_seq2) = cc->EvalSub(*(seqBufQer_ct + ct_seq2), encode_vector_to_plaintext(v_seq2));
			// vecInt v_remaining_2 = decrypt_ciphertext_to_plaintext_vector(*(seqBufQer_ct + ct_seq2));
			// for (int kk = ind_seq2; kk < num_of_ints_in_second_ct; kk++)
			// 	cout << v_remaining_2[kk] << " ";
			// cout << endl;
		}
		// cout << endl;

		sp.seqid = sp.regid = sp.score = sp.tle = sp.gtle = sp.qle = -1;
		sp.gscore = sp.max_off = -1;
		seqPairArray[numPairsRead] = sp;
		numPairsRead++;
		// SW_cells += (sp.len1 * sp.len2);
		// cout << "---------------------------------------------" << endl;
	}
}

// profiling stats
uint64_t find_stats(uint64_t *val, int nt, double &min, double &max, double &avg)
{
	min = 1e10;
	max = 0;
	avg = 0;
	for (int i = 0; i < nt; i++)
	{
		avg += val[i];
		if (max < val[i])
			max = val[i];
		if (min > val[i])
			min = val[i];
	}
	avg /= nt;

	return 1;
}

// ../benchmarks/bsw/bsw -pairs $INPUTS_DIR/bsw/small/bandedSWA_SRR7733443_100k_input.txt -t 1 -b 512

int main(int argc, char *argv[])
{
#ifdef VTUNE_ANALYSIS
	__itt_pause();
#endif
	if (argc < 3)
	{
		fprintf(stderr, "usage: bsw -pairs <InSeqFile> -t <threads> -b <batch_size>\n");
		exit(EXIT_FAILURE);
	}

	parseCmdLine(argc, argv);

	pairFile = fopen(pairFileName, "r");
	if (pairFile == NULL)
	{
		fprintf(stderr, "Could not open file: %s\n", pairFileName);
		exit(EXIT_FAILURE);
	}

	// (1024*1024)/16384=64
	const int bufSize = 1024 * 1024;
	char *buffer = (char *)malloc(bufSize * sizeof(char));
	vecCT enc_buf; // enc_buf.resize(64,0);

	size_t numLines = 0;
	size_t n;

	// size_t fread(void *ptr, size_t size, size_t count, FILE *stream);
	// ptr is a pointer to the buffer where the data will be stored
	// size is the size (in bytes) of each element to be read
	// count is the number of elements to be read
	// stream is a pointer to the input file stream
	// The function returns the total number of elements successfully read. If an error occurs, the return value will be less than count.

	int ct_count = 0;
	vecInt vv;
	int flg = 0;

	auto phase1_start=high_resolution_clock::now();

	while (n = fread(buffer, sizeof(char), bufSize, pairFile))
	{

		cout<<"n: "<<n<<endl;

		// enc_buf.resize(64,0);
		enc_buf.clear();

		// initializing this with enc('a') because if you initialize with
		ct_count = 0;

		for (int i = 0; i < n; i++) // so, inclduing all, enc_buf has n elements
		{
			if (vv.size() < 16384)
			{
				vv.push_back((int64_t)(buffer[i]));
			}
			else
			{
				// cout<<"i: "<<i<<"; ct_count: "<<ct_count<<endl;
				enc_buf.push_back(encrypt_plaintext_vector_to_ciphertext(vv));
				ct_count++;
				vv.clear();
				vv.push_back((int64_t)(buffer[i]));
			}
		}
		if (vv.size() != 0)
		{
			int sz = vv.size();
			if (sz != 16384)
				flg = 1;
			// cout<<"sz: "<<sz<<"; ct_count: "<<ct_count<<endl;
			vv.resize(16384); // Not needed as vv.size()==16384
			for (int kk = sz; kk < 16384; kk++)
				vv[kk] = 97;							// lowercase 'a'
			// for(int kk=0;kk<16384;kk++)
			// {
			// 	cout<<(char)(vv[kk])<<" ";
			// }
			// cout<<endl;
			enc_buf.push_back(encrypt_plaintext_vector_to_ciphertext(vv));
			ct_count++;
			vv.clear();
		}

		int cnt = 0;
		CT c1;
		vecInt v2;
		for (int i = 0; i < enc_buf.size(); i++)
		{ //	buf[i] corresponds to dec(enc_buf[cnt])[index]
			CT c1 = enc_buf[i];
			vecInt v1 = decrypt_ciphertext_to_plaintext_vector(c1);
			// cout<<"v1: "<<endl;
			// if(flg==1)
			// {
			// 	for(int j=0;j<16384;j++){
			// 		cout<<(char)(v1[j])<<" ";
			// 	}
			// }
			// cout<<endl;

			vecInt v2(16384, (int64_t)('\n'));
			// v2[index]=(int64_t)('\n');

			CT c2 = encrypt_plaintext_vector_to_ciphertext(v2);
			vecInt v3 = decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(c1, c2));
			// cout<<"v3: "<<endl;
			// for(int j=0;j<16384;j++){
			// 	cout<<(char)(v3[j])<<" ";
			// }
			// cout<<endl;

			for (int j = 0; j < 16384; j++)
			{
				if (v3[j] == 0)
					numLines++;
			}
		}

		// for (int i = 0; i < n; i++) {
		// 	// cout<<buffer[i]<<" ";
		//     if (buffer[i] == '\n') {
		//         numLines++;
		//     }
		// }
		// cout<<endl;

		// cout<<"\nnum: "<<numLines<<endl;
		// cout<<"------------------------------------------------\n";
	}

	auto phase1_end=high_resolution_clock::now();
	auto phase1_duration = duration_cast<microseconds>(phase1_end - phase1_start);
	printf("Time taken by phase1 is %lld microseconds\n",phase1_duration.count());

	enc_buf.resize(0);
	free(buffer);

	cout << "numLines: " << numLines << endl; 			// numLines: 300000

	// Reset file pointer back to the beginning
	fseek(pairFile, 0L, SEEK_SET);

	size_t numPairs = numLines / 3;
	numPairs-=50000;
	printf("numPairs orig: %d\n",numPairs);
	size_t roundNumPairs = ((numPairs + SIMD_WIDTH16 - 1) / SIMD_WIDTH16) * SIMD_WIDTH16;
	if (batchSize == 0)
	{
		batchSize = roundNumPairs;
	}

	printf("Number of input pairs: %ld\n", numPairs);
	// roundNumPairs-=50000;										// 100000
	printf("roundNumPairs: %zu\n", roundNumPairs);			

	size_t memAlloc = roundNumPairs * (sizeof(SeqPair) + MAX_SEQ_LEN_QER * sizeof(int8_t) + MAX_SEQ_LEN_REF * sizeof(int8_t));
	printf("Allocating %.3f GB memory for input buffers...\n", (memAlloc * 1.0) / (1024 * 1024 * 1024));
	SeqPair *seqPairArray = (SeqPair *)_mm_malloc(roundNumPairs * sizeof(SeqPair), 64); // roundNumPairs: 100000
	uint8_t *seqBufQer = (uint8_t *)_mm_malloc(MAX_SEQ_LEN_QER * roundNumPairs * sizeof(int8_t), 64);
	vecCT seqBufQer_enc;
	vecInt v1(16384, 65);
	seqBufQer_enc.resize(ceil(MAX_SEQ_LEN_QER * roundNumPairs * 1.0 / 16384), encrypt_plaintext_vector_to_ciphertext(v1));
	cout << "seqBufQer_enc.size(): " << seqBufQer_enc.size() << endl; // 1562
	uint8_t *seqBufRef = (uint8_t *)_mm_malloc(MAX_SEQ_LEN_REF * roundNumPairs * sizeof(int8_t), 64);
	vecCT seqBufRef_enc;
	vecInt v2(16384, 65);
	seqBufRef_enc.resize(ceil(MAX_SEQ_LEN_REF * roundNumPairs * 1.0 / 16384), encrypt_plaintext_vector_to_ciphertext(v2));
	cout << "seqBufRef_enc.size(): " << seqBufRef_enc.size() << endl; // 12500

	int8_t mat[25];
	bwa_fill_scmat(w_match, w_mismatch, w_ambig, mat);
	int zdrop = 100, w = 100, end_bonus = 5;

	// OutScore *outScoreArray = (OutScore *)_mm_malloc(MAX_NUM_PAIRS * sizeof(OutScore), 64);
	BandedPairWiseSW *bsw[numThreads];
	for (int i = 0; i < numThreads; i++)
	{
		bsw[i] = new BandedPairWiseSW(w_open, w_extend, w_open, w_extend,
									  zdrop, end_bonus, mat,
									  w_match, w_mismatch, 1);
	}

	int64_t startTick, totalTicks = 0, readTim = 0;

	uint64_t tim = __rdtsc();

	int max_seq1 = -10000;
	int max_seq2 = -10000;

	int64_t numPairsIndex = 0;
	auto st = high_resolution_clock::now();

	for (int64_t i = 0; i < roundNumPairs; i += batchSize) // batchSize=512; last i=99840
	{
		auto i_start = high_resolution_clock::now();

		cout << "i: " << i << "; numPairsIndex: " << numPairsIndex << endl;
		int nPairsBatch = (numPairs - i) >= batchSize ? batchSize : numPairs - i;
		loadPairs(seqPairArray + numPairsIndex, seqBufRef + numPairsIndex * MAX_SEQ_LEN_REF, &seqBufRef_enc[0], (numPairsIndex * MAX_SEQ_LEN_REF), seqBufQer + numPairsIndex * MAX_SEQ_LEN_QER, &seqBufQer_enc[0], (numPairsIndex * MAX_SEQ_LEN_QER), nPairsBatch);
		// loadPairs(seqPairArray + numPairsIndex, seqBufRef + numPairsIndex * MAX_SEQ_LEN_REF, seqBufQer + numPairsIndex * MAX_SEQ_LEN_QER, nPairsBatch);
		numPairsIndex += nPairsBatch;

		auto i_end = high_resolution_clock::now();

		auto duration = duration_cast<microseconds>(i_end - i_start);
		printf("Time taken by i: %d is %lld microseconds\n",i,duration.count());
	}

	auto en = high_resolution_clock::now();
	auto duration_loops = duration_cast<microseconds>(en - st);
	printf("Time taken by loops is %lld microseconds\n",duration_loops.count());

	readTim += __rdtsc() - tim;

	startTick = __rdtsc();

#ifdef VTUNE_ANALYSIS
	__itt_resume();
#endif
	int64_t workTicks[CLMUL * numThreads];
	memset(workTicks, 0, CLMUL * numThreads * sizeof(int64_t));
#pragma omp parallel num_threads(numThreads)
	{
		int tid = omp_get_thread_num();
		auto scores_start=high_resolution_clock::now();

#pragma omp for schedule(dynamic, 1)
		for (int64_t i = 0; i < roundNumPairs; i += batchSize)
		{
			auto scores_iter_start=high_resolution_clock::now();
			
			int nPairsBatch = (numPairs - i) >= batchSize ? batchSize : numPairs - i;
			int64_t st1 = __rdtsc();
			// bsw[tid]->print_ct(100);
			bsw[tid]->getScores16(seqPairArray + i, seqBufRef + i * MAX_SEQ_LEN_REF, &seqBufRef_enc[0], i * MAX_SEQ_LEN_REF, seqBufQer + i * MAX_SEQ_LEN_QER, &seqBufQer_enc[0], i * MAX_SEQ_LEN_QER, nPairsBatch, 1, w);
            // bsw[tid]->getScores16(seqPairArray + i, seqBufRef + i * MAX_SEQ_LEN_REF, seqBufQer + i * MAX_SEQ_LEN_QER, nPairsBatch, 1, w);
			int64_t et1 = __rdtsc();
			workTicks[CLMUL * tid] += (et1 - st1);

			auto scores_iter_end=high_resolution_clock::now();
			auto scores_iter_duration = duration_cast<microseconds>(scores_iter_end - scores_iter_start);
			printf("Time taken by scores phase iteration #%d is %lld microseconds\n",i,scores_iter_duration.count());
		}

		auto scores_end=high_resolution_clock::now();
		auto scores_duration = duration_cast<microseconds>(scores_end - scores_start);
		printf("Time taken by scores phase is %lld microseconds\n",scores_duration.count());

		printf("%d] workTicks = %ld\n", tid, workTicks[CLMUL * tid]);
	}

#ifdef VTUNE_ANALYSIS
	__itt_pause();
#endif
	totalTicks += __rdtsc() - startTick;
	printf("Executed AVX2 vector code...\n");

	tim = __rdtsc();
	sleep(1);
	freq = __rdtsc() - tim;

	printf("Processor freq: %0.2lf MHz\n", freq / 1e6);

	// int64_t myTicks = bsw->getTicks();
	printf("Read time = %0.2lf s\n", readTim / freq);
	printf("Overall SW cycles = %ld, %0.2lf s\n", totalTicks, totalTicks * 1.0 / freq);
	printf("Total Pairs processed: %d\n", numPairs);

	int64_t sumTicks = 0;
	int64_t maxTicks = 0;
	for (int i = 0; i < numThreads; i++)
	{
		sumTicks += workTicks[CLMUL * i];
		if (workTicks[CLMUL * i] > maxTicks)
			maxTicks = workTicks[CLMUL * i];
	}
	double avgTicks = (sumTicks * 1.0) / numThreads;
	printf("avgTicks = %lf, maxTicks = %ld, load imbalance = %lf\n", avgTicks, maxTicks, maxTicks / avgTicks);

	// printf("SW cells(T)  = %ld\n", SW_cells);
	// printf("SW cells(||)  = %ld\n", SW_cells2);
	// printf("SW GCUPS  = %lf\n", SW_cells * freq/1e9 / myTicks);

	// printf("More stats:\n");
	// double freq = 2.3*1e9;
	// double min, max, avg;
	// find_stats(prof[1], numThreads, min, max, avg);
	// printf("Time in pre-processing: %0.2lf (%0.2lf, %0.2lf)\n",
	// avg*1.0/freq, min*1.0/freq, max*1.0/freq);
	// find_stats(prof[0], numThreads, min, max, avg);
	// printf("Time spent in smithWaterman(): %0.2lf (%0.2lf, %0.2lf)\n",
	// 	avg*1.0/freq, min*1.0/freq, max*1.0/freq);

	// printf("\nDebugging info:\n");
	// printf("Time taken for DP loop: %0.2lf\n", prof[DP][0]*1.0/freq);
	// printf("Time taken for DP loop upper part: %0.2lf\n", prof[DP3][0]*1.0/freq);
	// printf("Time taken for DP inner loop: %0.2lf\n", prof[DP1][0]*1.0/freq);
	// printf("Time taken for DP loop lower part: %0.2lf\n", prof[DP2][0]*1.0/freq);

	/**** free memory *****/
	_mm_free(seqPairArray);
	_mm_free(seqBufRef);
	_mm_free(seqBufQer);
	// bsw->getTicks();
	for (int i = 0; i < numThreads; i++)
	{
		bsw[i]->getTicks();
		delete bsw[i];
	}

	fclose(pairFile);
	return 1;
}
