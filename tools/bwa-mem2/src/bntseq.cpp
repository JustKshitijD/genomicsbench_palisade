/* The MIT License

   Copyright (c) 2008 Genome Research Ltd (GRL).

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

   Modified Copyright (C) 2020 Intel Corporation, Heng Li.
   Contacts: Vasimuddin Md <vasimuddin.md@intel.com>; Sanchit Misra <sanchit.misra@intel.com>;
   Heng Li <hli@jimmy.harvard.edu>
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <unistd.h>
#include <errno.h>
#include <limits.h>
#include "bntseq.h"
#include "utils.h"
#include "macro.h"

#include <chrono>
using namespace std::chrono;

#include "kseq.h"
KSEQ_DECLARE(gzFile)

#include "khash.h"
KHASH_MAP_INIT_STR(str, int)

#ifdef USE_MALLOC_WRAPPERS
#include "malloc_wrap.h"
#endif

#include "../../../palisade_header.h"

#include <chrono>
using namespace std::chrono;

extern uint64_t tprof[LIM_R][LIM_C];

unsigned char nst_nt4_table[256] = {
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5 /*'-'*/, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};

void bns_dump(const bntseq_t *bns, const char *prefix)
{
	char str[PATH_MAX];
	FILE *fp;
	int i;
	// assert(strlen(prefix) + 4 < 1024);
	{ // dump .ann
		strcpy_s(str, PATH_MAX, prefix);
		strcat_s(str, PATH_MAX, ".ann");
		fp = xopen(str, "w");
		err_fprintf(fp, "%lld %d %u\n", (long long)bns->l_pac, bns->n_seqs, bns->seed);
		for (i = 0; i != bns->n_seqs; ++i)
		{
			bntann1_t *p = bns->anns + i;
			err_fprintf(fp, "%d %s", p->gi, p->name);
			if (p->anno[0])
				err_fprintf(fp, " %s\n", p->anno);
			else
				err_fprintf(fp, "\n");
			err_fprintf(fp, "%lld %d %d\n", (long long)p->offset, p->len, p->n_ambs);
		}
		err_fflush(fp);
		err_fclose(fp);
	}
	{ // dump .amb
		strcpy_s(str, PATH_MAX, prefix);
		strcat_s(str, PATH_MAX, ".amb");
		fp = xopen(str, "w");
		err_fprintf(fp, "%lld %d %u\n", (long long)bns->l_pac, bns->n_seqs, bns->n_holes);
		for (i = 0; i != bns->n_holes; ++i)
		{
			bntamb1_t *p = bns->ambs + i;
			err_fprintf(fp, "%lld %d %c\n", (long long)p->offset, p->len, p->amb);
		}
		err_fflush(fp);
		err_fclose(fp);
	}
}

// bntseq_t *bns_restore_core(const char *ann_filename, const char *amb_filename, const char *pac_filename)
// {
// 	char str[8193];
// 	FILE *fp;
// 	const char *fname;
// 	bntseq_t *bns;
// 	long long xx;
// 	int i;
// 	int scanres;
// 	bns = (bntseq_t *)calloc(1, sizeof(bntseq_t));
// 	assert(bns != 0);
// 	{ // read .ann
// 		fp = xopen(fname = ann_filename, "r");
// 		scanres = fscanf(fp, "%lld%d%u", &xx, &bns->n_seqs, &bns->seed);
// 		assert(bns->n_seqs >= 0 && bns->n_seqs <= INT_MAX);
// 		if (scanres != 3)
// 			goto badread;
// 		bns->l_pac = xx;
// 		bns->anns = (bntann1_t *)calloc(bns->n_seqs, sizeof(bntann1_t));
// 		assert(bns->anns != NULL);
// 		for (i = 0; i < bns->n_seqs; ++i)
// 		{
// 			printf("i: %d\n", i);

// 			bntann1_t *p = bns->anns + i;
// 			char *q = str;
// 			int c;
// 			// read gi and sequence name
// 			scanres = fscanf(fp, "%u%8192s", &p->gi, str);
// 			if (scanres != 2)
// 				goto badread;
// 			p->name = strdup(str);
// 			// read fasta comments
// 			printf("sizeof(str): %lu\n", sizeof(str));
// 			int ind1 = 0;
// 			int ind2 = 0;
// 			printf("strlen(str): %d\n", strlen(str));
// 			printf("str: %s\n", str);

// 			while (q - str < sizeof(str) - 1 && (c = fgetc(fp)) != '\n' && c != EOF)
// 			{
// 				printf("ind1: %d\n", ind1);
// 				ind1++;
// 				*q++ = c;
// 				printf("c: %d; %c\n", c, c);
// 				// cout<<"c: "<<c<<"; "<<char(c)<<endl;
// 			}

// 			while (c != '\n' && c != EOF)
// 			{
// 				printf("ind2: %d\n", ind2);
// 				c = fgetc(fp);
// 				printf("c: %d; %c\n", c, c);
// 				ind2++;
// 			}

// 			if (c == EOF)
// 			{
// 				scanres = EOF;
// 				goto badread;
// 			}
// 			*q = 0;
// 			assert(strlen(str) < 8192);
// 			if (q - str > 1 && strcmp(str, " (null)") != 0)
// 				p->anno = strdup(str + 1); // skip leading space
// 			else
// 				p->anno = strdup("");
// 			// read the rest
// 			scanres = fscanf(fp, "%lld%d%d", &xx, &p->len, &p->n_ambs);
// 			if (scanres != 3)
// 				goto badread;
// 			p->offset = xx;
// 		}
// 		err_fclose(fp);
// 	}
// 	{ // read .amb
// 		int64_t l_pac;
// 		int32_t n_seqs;
// 		fp = xopen(fname = amb_filename, "r");
// 		scanres = fscanf(fp, "%lld%d%d", &xx, &n_seqs, &bns->n_holes);
// 		assert(bns->n_holes >= 0 && bns->n_holes <= INT_MAX);
// 		if (scanres != 3)
// 			goto badread;
// 		l_pac = xx;
// 		xassert(l_pac == bns->l_pac && n_seqs == bns->n_seqs, "inconsistent .ann and .amb files.");
// 		if (bns->n_holes)
// 		{
// 			bns->ambs = (bntamb1_t *)calloc(bns->n_holes, sizeof(bntamb1_t));
// 			assert(bns->ambs != NULL);
// 		}
// 		else
// 		{
// 			bns->ambs = 0;
// 		}
// 		for (i = 0; i < bns->n_holes; ++i)
// 		{
// 			bntamb1_t *p = bns->ambs + i;
// 			scanres = fscanf(fp, "%lld%d%8192s", &xx, &p->len, str);
// 			if (scanres != 3)
// 				goto badread;
// 			p->offset = xx;
// 			p->amb = str[0];
// 		}
// 		err_fclose(fp);
// 	}
// 	{ // open .pac
// 		bns->fp_pac = xopen(pac_filename, "rb");
// 	}
// 	return bns;

// badread:
// 	if (EOF == scanres)
// 	{
// 		err_fatal(__func__, "Error reading %s : %s\n", fname, ferror(fp) ? strerror(errno) : "Unexpected end of file");
// 	}
// 	err_fatal(__func__, "Parse error reading %s\n", fname);
// }

// bntseq_t *bns_restore(const char *prefix)
// {
// 	printf("prefix: %s\n",prefix);
// 	char ann_filename[PATH_MAX], amb_filename[PATH_MAX], pac_filename[PATH_MAX], alt_filename[PATH_MAX];
// 	FILE *fp;
// 	bntseq_t *bns;
// 	//assert(strlen(prefix) + 4 < 1024);
// 	strcpy_s(ann_filename, PATH_MAX, prefix); strcat_s(ann_filename, PATH_MAX, ".ann");
// 	strcpy_s(amb_filename, PATH_MAX, prefix); strcat_s(amb_filename, PATH_MAX, ".amb");
// 	strcpy_s(pac_filename, PATH_MAX, prefix); strcat_s(pac_filename, PATH_MAX, ".pac");
// 	bns = bns_restore_core(ann_filename, amb_filename, pac_filename);
// 	if (bns == 0) return 0;
//     strcpy_s(alt_filename, PATH_MAX, prefix); strcat_s(alt_filename, PATH_MAX, ".alt");
// 	if ((fp = fopen(alt_filename, "r")) != 0) { // read .alt file if present
// 		char str[1024];
// 		khash_t(str) *h;									// #define khash_t(name) kh_##name##_t

// 		// typedef unsigned int khint32_t;
// 		// typedef khint32_t khint_t;
// 		// typedef khint_t khiter_t;

// 		int c, i, absent;
// 		khint_t k;
// 		h = kh_init(str);
//         assert(h != NULL);
// 		for (i = 0; i < bns->n_seqs; ++i) {
// 			k = kh_put(str, h, bns->anns[i].name, &absent);
// 			kh_val(h, k) = i;
// 		}
// 		i = 0;
// 		while ((c = fgetc(fp)) != EOF) {
// 			if (c == '\t' || c == '\n' || c == '\r') {
// 				str[i] = 0;
// 				if (str[0] != '@') {
// 					k = kh_get(str, h, str);
// 					if (k != kh_end(h))
// 						bns->anns[kh_val(h, k)].is_alt = 1;
// 				}
// 				while (c != '\n' && c != EOF) c = fgetc(fp);
// 				i = 0;
// 			} else str[i++] = c; // FIXME: potential segfault here
// 		}
// 		kh_destroy(str, h);
// 		fclose(fp);
// 	}
// 	return bns;
// }

// bntseq_t *bns_restore_core(const char *ann_filename, const char *amb_filename, const char *pac_filename)
// {
// 	printf("In bns_retore_core\n");
// 	char str[8193];
// 	vecCT str_enc; // str_enc.resize(8193);
// 	FILE *fp;
// 	const char *fname;
// 	bntseq_t *bns;
// 	long long xx;
// 	CT xx_enc;
// 	int i;
// 	int scanres;
// 	bns = (bntseq_t *)calloc(1, sizeof(bntseq_t));
// 	assert(bns != 0);
// 	{ // read .ann
// 		fp = xopen(fname = ann_filename, "r");
// 		// int32_t bns_nseqs;
// 		uint32_t bns_seed;
// 		// scanres = fscanf(fp, "%lld%d%u", &xx, &bns->n_seqs, &bns->seed);	// bns->n_seqs and bns->seed read from ann_file
// 		scanres = fscanf(fp, "%lld%d%u", &xx, &bns->n_seqs, &bns_seed);
// 		xx_enc = encrypt_plaintext_integer_to_ciphertext(xx);
// 		bns->n_seqs_enc = encrypt_plaintext_integer_to_ciphertext(bns->n_seqs);
// 		// bns->seed_enc = encrypt_plaintext_integer_to_ciphertext(bns_seed);
// 		// printf("decrypt_ciphertext_to_plaintext_vector(xx_enc)[0]: %d, decrypt_ciphertext_to_plaintext_vector(bns->n_seqs_enc)[0]: %d, decrypt_ciphertext_to_plaintext_vector(bns->seed_enc)[0]: %d\n", decrypt_ciphertext_to_plaintext_vector(xx_enc)[0], decrypt_ciphertext_to_plaintext_vector(bns->n_seqs_enc)[0], decrypt_ciphertext_to_plaintext_vector(bns->seed_enc)[0]);
// 		assert(operate_and_decrypt(bns->n_seqs_enc, "-", 0) >= 0 && operate_and_decrypt(bns->n_seqs_enc, "-", INT_MAX) <= 0);
// 		// assert(decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(bns->n_seqs_enc, encode_integer_to_plaintext(0)))[0]>=0 && decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(bns->n_seqs_enc,encode_integer_to_plaintext(INT_MAX)))[0]<=0);
// 		if (scanres != 3)
// 			goto badread;
// 		bns->l_pac = xx;
// 		bns->l_pac_enc = xx_enc;
// 		cout << "bns->n_seqs: " << bns->n_seqs << endl;
// 		bns->anns = (bntann1_t *)calloc(decrypt_ciphertext_to_plaintext_vector(bns->n_seqs_enc)[0], sizeof(bntann1_t));
// 		assert(bns->anns != NULL);

// 		ofstream outfile_p_gi("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/p_gi.txt");
// 		ofstream outfile_p_str("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/p_str.txt");
// 		ofstream outfile_p_anno("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/p_anno.txt");
// 		ofstream outfile_p_len("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/p_len.txt");
// 		ofstream outfile_p_n_ambs("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/p_n_ambs.txt");
// 		ofstream outfile_p_offset("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/p_offset.txt");

// 		for (i = 0; i < decrypt_ciphertext_to_plaintext_vector(bns->n_seqs_enc)[0]; ++i)
// 		{
// 			cout << "-----------------------------------------------------" << endl;
// 			cout << "i: " << i << endl;
// 			bntann1_t *p = &bns->anns[i];
// 			char *q = str;
// 			int c;
// 			uint32_t gi;
// 			scanres = fscanf(fp, "%u%8192s", &gi, str); // reads 8192 characters; the 8193th character is made '\0'
// 			// if (scanres != 2) goto badread;

// 			str_enc.resize(8193,0);
// 			// So, strlen(str)=8192
// 			p->gi = gi;
// 			// p->gi_enc = encrypt_plaintext_integer_to_ciphertext(gi);
// 			outfile_p_gi << gi << endl;

// 			// After below, str_enc gets all till 0s
// 			for (int h = 0; h < 8193; h++) // copy from str to str_enc till str[i]!=0
// 			{
// 				// if (str[h] == 0)
// 				// {
// 				// 	str_enc[h]=encrypt_plaintext_integer_to_ciphertext(0);
// 				// 	cout << "Zero encountered h: " << h << "; str[h]: " << int(str[h]) <<"; decrypt_ciphertext_to_plaintext_vector(str_enc[h])[0]: "<<(char)(decrypt_ciphertext_to_plaintext_vector(str_enc[h])[0])<< endl;

// 				// 	break;
// 				// }

// 				// if (h + 1 > str_enc.size())
// 				// 	str_enc.push_back(encrypt_plaintext_integer_to_ciphertext(str[h]));
// 				// else
// 				str_enc[h] = encrypt_plaintext_integer_to_ciphertext(str[h]);
// 				// if(str[h]==0)
// 				// {
// 				// 	str_enc[h]=encrypt_plaintext_integer_to_ciphertext('0');
// 				// }

// 				cout << "h: " << h << "; str[h]: " << int(str[h]) << "; decrypt_ciphertext_to_plaintext_vector(str_enc[h])[0]: " << decrypt_ciphertext_to_plaintext_vector(str_enc[h])[0] << endl;

// 				// if(h==sizeof(str)-1)
// 				// {
// 				// 	if (h + 2 > str_enc.size())
// 				// 		str_enc.push_back(0);
// 				// 	else
// 				// 		str_enc[h+1] = 0;

// 				// 	cout << "Last h with 0, decrypt_ciphertext_to_plaintext_vector(str_enc[h+1])[0]: "<<decrypt_ciphertext_to_plaintext_vector(str_enc[h+1])[0]<< endl;

// 				// }
// 			}

// 			// cout << "3 p->name_enc.size(): " << p->name_enc.size() << endl;

// 			if (scanres != 2)
// 				goto badread;
// 			p->name = strdup(str);

// 			printf("p->name: %s; strlen(p->name): %d\n", p->name, strlen(p->name));
// 			printf("convert_ciphertext_vector_to_plaintext_string(str_enc): %s\n", convert_ciphertext_vector_to_plaintext_string(str_enc));
// 			cout << "str_enc.size(): " << str_enc.size() << "; strlen_string_enc(str_enc): " << strlen_string_enc(str_enc) << endl;
// 			// p->name_enc=vecCT(1,0);
// 			// printf("&p->name_enc: %p\n",&p->name_enc);
// 			// cout<<"4 p->name_enc.size(): "<<p->name_enc.size()<<endl; //<<"; strlen_enc(p->name_enc): "<<strlen_enc(p->name_enc)<<endl;

// 			vecCT name_enc;
// 			name_enc.resize(0);
// 			strdup_enc(str_enc, name_enc);							// str_enc has 0/NULL at end
// 			cout << "name_enc.size(): " << name_enc.size() <<"; str_enc.size(): "<<str_enc.size()<<endl; //<<"; strlen_enc(p->name_enc): "<<strlen_enc(p->name_enc)<<endl;
// 			printf("convert_ciphertext_vector_to_plaintext_string(name_enc): %s\n", convert_ciphertext_vector_to_plaintext_string(name_enc));

// 			// for (int kk = 0; kk < name_enc.size(); kk++)
// 			// {
// 			// 	if (!Serial::SerializeToFile("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/serializer_folder/p_name/" + to_string(i) + "/" + to_string(kk) + ".txt", str_enc[kk], SerType::BINARY))
// 			// 	{
// 			// 		std::cerr
// 			// 			<< "Error writing serialization of ciphertext_bns_anns_" + to_string(i)
// 			// 			<< std::endl;
// 			// 		return;
// 			// 	}
// 			// }
// 			outfile_p_str << str << endl;

// 			// p->name_enc.resize(str_enc.size()); // p->name_enc.shrink_to_fit();

// 			// for(int i=0;i<str_enc.size();i++)
// 			// 	p->name_enc[i]=str_enc[i];

// 			// read fasta comments
// 			// while (q - str < sizeof(str) - 1 && (c = fgetc(fp)) != '\n' && c != EOF) *q++ = c;

// 			int i_comment = 0;
// 			CT c_enc;
// 			// extra -1 because str_enc.size()=strlen(str)+1, due to extra enc('\n')
// 			while (i_comment < (str_enc).size() - 1 && compare_enc((c_enc = encrypt_plaintext_integer_to_ciphertext(fgetc(fp))), '\n') == 0 && compare_enc(c_enc, EOF) == 0)
// 			{
// 				cout << "i_comment: " << i_comment << endl;
// 				printf("c: %d; %c\n", decrypt_ciphertext_to_plaintext_vector(c_enc)[0], decrypt_ciphertext_to_plaintext_vector(c_enc)[0]);
// 				str_enc[i_comment++] = c_enc;
// 			}

// 			cout << "Done1" << endl;
// 			// printf("out decrypt_ciphertext_to_plaintext_vector(c_enc)[0]: %d\n", decrypt_ciphertext_to_plaintext_vector(c_enc)[0]);

// 			// while (c != '\n' && c != EOF) c = fgetc(fp);
// 			int ind2 = 0;
// 			while (compare_enc(c_enc, '\n') == 0 && compare_enc(c_enc, EOF) == 0)
// 			{
// 				cout << "ind2: " << ind2 << endl;
// 				// cout << "In loop" << endl;
// 				c_enc = encrypt_plaintext_integer_to_ciphertext(fgetc(fp));
// 				printf("decrypt_ciphertext_to_plaintext_vector(c_enc)[0]: %d\n", decrypt_ciphertext_to_plaintext_vector(c_enc)[0]);
// 				ind2++;
// 			}
// 			if (compare_enc(c_enc, EOF) == 1)
// 			{ // returns 1 when c_enc and EOF are equal
// 				scanres = EOF;
// 				printf("Bad decrypt_ciphertext_to_plaintext_vector(c_enc)[0]: %d\n", decrypt_ciphertext_to_plaintext_vector(c_enc)[0]);
// 				goto badread;
// 			}
// 			// *q = 0;
// 			str_enc[i_comment] = encrypt_plaintext_integer_to_ciphertext('\0');
// 			for (int h = 0; h < i_comment; h++)
// 				str[h] = decrypt_ciphertext_to_plaintext_vector(str_enc[h])[0];
// 			str[i_comment] = 0;

// 			cout << "strlen(new str): " << strlen(str) << endl;
// 			printf("new str: %s\n", str);
// 			printf("new str_enc: %s\n", convert_ciphertext_vector_to_plaintext_string(str_enc));
// 			cout << "strlen_enc(new str_enc): " << strlen_enc(str_enc) << endl;
// 			// cout << "strlen_string_enc(new str_enc): " << strlen_string_enc(str_enc) << endl;

// 			// assert(strlen(str) < 8192);
// 			assert(strlen_string_enc(str_enc) < 8192);
// 			vecCT p_anno_enc;

// 			// if (q - str > 1 && strcmp(str, " (null)") != 0) p->anno = strdup(str + 1); // skip leading space
// 			// else p->anno = strdup("");
// 			// vecCT str_enc_2; str_enc_2.assign(str_enc.begin()+1,str_enc.end());
// 			// if(i_comment>1 &&  strcmp_enc(str_enc," (null)")!=0 ) strdup_enc(str_enc_2,p->anno_enc);
// 			if (i_comment > 1 && strcmp_enc(str_enc, " (null)") != 0)
// 				assign_string_to_vecCT(p_anno_enc, str + 1, -1);
// 			else
// 				assign_string_to_vecCT(p_anno_enc, "", -1);

// 			printf("p_anno_enc: %s\n", convert_ciphertext_vector_to_plaintext_string(p_anno_enc));

// 			// for (int kk = 0; kk < p_anno_enc.size(); kk++)
// 			// {
// 			// 	if (!Serial::SerializeToFile("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/serializer_folder/p_anno_enc/" + to_string(i) + "/" + to_string(kk) + ".txt", p_anno_enc[kk], SerType::BINARY))
// 			// 	{
// 			// 		std::cerr
// 			// 			<< "Error writing serialization of ciphertext_bns_anns_" + to_string(i)
// 			// 			<< std::endl;
// 			// 		return;
// 			// 	}
// 			// }
// 			outfile_p_anno<<convert_ciphertext_vector_to_plaintext_string(p_anno_enc)<<endl;

// 			// cout<<"Earlier p->anno_enc.size(): "<<p->anno_enc.size()<<endl;
// 			// p->anno_enc.resize(strlen_enc(p->anno_enc)+1);
// 			// cout<<"Later p->anno_enc.size(): "<<p->anno_enc.size()<<endl;

// 			// else assign_string_to_vecCT(p->anno_enc,"",-1);
// 			int32_t p_len, p_n_ambs;
// 			cout << "Reading rest of broad.ann" << endl;
// 			// read the rest
// 			scanres = fscanf(fp, "%lld%d%d", &xx, &p_len, &p_n_ambs);
// 			cout << "Read broad.ann" << endl;
// 			// xx_enc = encrypt_plaintext_integer_to_ciphertext(xx);

// 			// if (!Serial::SerializeToFile("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/serializer_folder/p_len/" + to_string(i) + ".txt", encrypt_plaintext_integer_to_ciphertext(p_len), SerType::BINARY))
// 			// {
// 			// 	std::cerr
// 			// 		<< "Error writing serialization of ciphertext_bns_anns_" + to_string(i)
// 			// 		<< std::endl;
// 			// 	return;
// 			// }
// 			// if (!Serial::SerializeToFile("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/serializer_folder/p_n_ambs/" + to_string(i) + ".txt", encrypt_plaintext_integer_to_ciphertext(p_n_ambs), SerType::BINARY))
// 			// {
// 			// 	std::cerr
// 			// 		<< "Error writing serialization of ciphertext_bns_anns_" + to_string(i)
// 			// 		<< std::endl;
// 			// 	return;
// 			// }

// 			outfile_p_len<<p_len<<endl;
// 			outfile_p_n_ambs<<p_n_ambs<<endl;

// 			// p->n_ambs_enc = encrypt_plaintext_integer_to_ciphertext(p_n_ambs);
// 			if (scanres != 3)
// 				goto badread;
// 			p->offset = xx;

// 			// if (!Serial::SerializeToFile("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/serializer_folder/p_offset/" + to_string(i) + ".txt", encrypt_plaintext_integer_to_ciphertext(xx), SerType::BINARY))
// 			// {
// 			// 	std::cerr
// 			// 		<< "Error writing serialization of ciphertext_bns_anns_" + to_string(i)
// 			// 		<< std::endl;
// 			// 	return;
// 			// }
// 			outfile_p_offset<<xx<<endl;

// 			// p->offset_enc = xx_enc;

// 			// if (!Serial::SerializeToFile("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/serializer_folder/ciphertext_bns_anns_" + to_string(i) + ".txt", c1, SerType::BINARY))
// 			// {
// 			// 	std::cerr
// 			// 		<< "Error writing serialization of ciphertext_bns_anns_" + to_string(i)
// 			// 		<< std::endl;
// 			// 	return;
// 			// }
// 		}
// 		err_fclose(fp);
// 	}
// 	{ // read .amb
// 		int64_t l_pac;
// 		CT l_pac_enc;
// 		int32_t n_seqs;
// 		CT n_seqs_enc;
// 		int32_t n_holes;
// 		fp = xopen(fname = amb_filename, "r");
// 		scanres = fscanf(fp, "%lld%d%d", &xx, &n_seqs, &n_holes);
// 		bns->n_holes = n_holes;
// 		xx_enc = encrypt_plaintext_integer_to_ciphertext(xx);
// 		n_seqs_enc = encrypt_plaintext_integer_to_ciphertext(n_seqs);
// 		bns->n_holes_enc = encrypt_plaintext_integer_to_ciphertext(n_holes);
// 		bns->n_seqs_enc = encrypt_plaintext_integer_to_ciphertext(n_seqs);

// 		// assert(bns->n_holes >= 0 && bns->n_holes <= INT_MAX);
// 		assert(decrypt_ciphertext_to_plaintext_vector(bns->n_holes_enc)[0] >= 0 && operate_and_decrypt(bns->n_holes_enc, "-", INT_MAX) <= 0);
// 		if (scanres != 3)
// 			goto badread;
// 		// l_pac = xx;
// 		l_pac_enc = xx_enc;
// 		// xassert(l_pac == bns->l_pac && n_seqs == bns->n_seqs, "inconsistent .ann and .amb files.");
// 		xassert(compare_enc(l_pac_enc, bns->l_pac_enc) && compare_enc(n_seqs_enc, bns->n_seqs_enc), "inconsistent .ann and .amb files.");
// 		if (compare_enc(bns->n_holes_enc, 1) == 1)
// 		{
// 			// bns->ambs = (bntamb1_t*)calloc(bns->n_holes, sizeof(bntamb1_t));
// 			bns->ambs = (bntamb1_t *)calloc(decrypt_ciphertext_to_plaintext_vector(bns->n_holes_enc)[0], sizeof(bntamb1_t));
// 			assert(bns->ambs != NULL);
// 		}
// 		else
// 		{
// 			bns->ambs = 0;
// 		}

// 		cout<<"n_holes: "<<n_holes<<endl;

// 		ofstream outfile_ambs_len("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/ambs_len.txt");
// 		ofstream outfile_ambs_offset("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/ambs_offset.txt");
// 		ofstream outfile_ambs_amb("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/ambs_amb.txt");

// 		for (i = 0; operate_and_decrypt(bns->n_holes_enc, "-", i) > 0; ++i)
// 		{
// 			bntamb1_t *p = bns->ambs + i;
// 			int32_t p_len;
// 			scanres = fscanf(fp, "%lld%d%8192s", &xx, p_len, str);
// 			xx_enc = encrypt_plaintext_integer_to_ciphertext(xx);

// 			// if (!Serial::SerializeToFile("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/serializer_folder/ambs_len/" + to_string(i) + ".txt", encrypt_plaintext_integer_to_ciphertext(p_len), SerType::BINARY))
// 			// {
// 			// 	std::cerr
// 			// 		<< "Error writing serialization of ciphertext_bns_anns_" + to_string(i)
// 			// 		<< std::endl;
// 			// 	return;
// 			// }

// 			outfile_ambs_len<<p_len<<endl;

// 			// p->len_enc = encrypt_plaintext_integer_to_ciphertext(p_len);
// 			for (int h = 0; h < 8193; h++)
// 				str_enc[h] = encrypt_plaintext_integer_to_ciphertext(str[h]);

// 			if (scanres != 3)
// 				goto badread;
// 			// p->offset = xx;
// 			// p->amb = str[0];

// 			// if (!Serial::SerializeToFile("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/serializer_folder/ambs_offset/" + to_string(i) + ".txt", xx_enc, SerType::BINARY))
// 			// {
// 			// 	std::cerr
// 			// 		<< "Error writing serialization of ciphertext_bns_anns_" + to_string(i)
// 			// 		<< std::endl;
// 			// 	return;
// 			// }
// 			// if (!Serial::SerializeToFile("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/serializer_folder/ambs_amb/" + to_string(i) + ".txt", str_enc[0], SerType::BINARY))
// 			// {
// 			// 	std::cerr
// 			// 		<< "Error writing serialization of ciphertext_bns_anns_" + to_string(i)
// 			// 		<< std::endl;
// 			// 	return;
// 			// }

// 			outfile_ambs_offset<<xx<<endl;
// 			outfile_ambs_amb<<str[0]<<endl;

// 			// p->offset_enc = xx_enc;
// 			// p->amb_enc = str_enc[0];
// 		}
// 		err_fclose(fp);
// 	}
// 	{ // open .pac
// 		bns->fp_pac = xopen(pac_filename, "rb");
// 	}
// 	return bns;

// badread:
// 	if (EOF == scanres)
// 	{
// 		err_fatal(__func__, "Error reading %s : %s\n", fname, ferror(fp) ? strerror(errno) : "Unexpected end of file");
// 	}
// 	err_fatal(__func__, "Parse error reading %s\n", fname);
// }

bntseq_t *bns_restore_core(const char *ann_filename, const char *amb_filename, const char *pac_filename)
{
	printf("In bns_retore_core\n");
	char str[8193];
	vecCT str_enc; // str_enc.resize(8193);
	FILE *fp;
	const char *fname;
	bntseq_t *bns;
	long long xx;
	CT xx_enc;
	int i;
	int scanres;
	bns = (bntseq_t *)calloc(1, sizeof(bntseq_t));

	assert(bns != 0);
	{ // read .ann
		fp = xopen(fname = ann_filename, "r");
		// int32_t bns_nseqs;
		uint32_t bns_seed;
		// scanres = fscanf(fp, "%lld%d%u", &xx, &bns->n_seqs, &bns->seed);	// bns->n_seqs and bns->seed read from ann_file
		scanres = fscanf(fp, "%lld%d%u", &xx, &bns->n_seqs, &bns_seed);
		xx_enc = encrypt_plaintext_integer_to_ciphertext(xx);
		bns->n_seqs_enc = encrypt_plaintext_integer_to_ciphertext(bns->n_seqs);
		// bns->seed_enc = encrypt_plaintext_integer_to_ciphertext(bns_seed);
		// printf("decrypt_ciphertext_to_plaintext_vector(xx_enc)[0]: %d, decrypt_ciphertext_to_plaintext_vector(bns->n_seqs_enc)[0]: %d, decrypt_ciphertext_to_plaintext_vector(bns->seed_enc)[0]: %d\n", decrypt_ciphertext_to_plaintext_vector(xx_enc)[0], decrypt_ciphertext_to_plaintext_vector(bns->n_seqs_enc)[0], decrypt_ciphertext_to_plaintext_vector(bns->seed_enc)[0]);
		assert(operate_and_decrypt(bns->n_seqs_enc, "-", 0) >= 0 && operate_and_decrypt(bns->n_seqs_enc, "-", INT_MAX) <= 0);
		// assert(decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(bns->n_seqs_enc, encode_integer_to_plaintext(0)))[0]>=0 && decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(bns->n_seqs_enc,encode_integer_to_plaintext(INT_MAX)))[0]<=0);
		if (scanres != 3)
			goto badread;
		bns->l_pac = xx;
		bns->l_pac_enc = xx_enc;
		// cout << "bns->n_seqs: " << bns->n_seqs << endl;
		bns->anns = (bntann1_t *)calloc(decrypt_ciphertext_to_plaintext_vector(bns->n_seqs_enc)[0], sizeof(bntann1_t));
		for(int kk=0;kk<decrypt_ciphertext_to_plaintext_vector(bns->n_seqs_enc)[0];kk++)
			bns->anns[kk].is_alt_enc=encrypt_plaintext_integer_to_ciphertext(0);

		assert(bns->anns != NULL);

		ofstream outfile_p_gi("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/p_gi.txt");
		ofstream outfile_p_str("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/p_str.txt");
		ofstream outfile_p_anno("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/p_anno.txt");
		ofstream outfile_p_len("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/p_len.txt");
		ofstream outfile_p_n_ambs("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/p_n_ambs.txt");
		ofstream outfile_p_offset("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/p_offset.txt");

		printf("decrypt_ciphertext_to_plaintext_vector(bns->n_seqs_enc)[0]: %d\n",decrypt_ciphertext_to_plaintext_vector(bns->n_seqs_enc)[0]);

		for (i = 0; i < decrypt_ciphertext_to_plaintext_vector(bns->n_seqs_enc)[0]; ++i)
		{
			// cout << "-----------------------------------------------------" << endl;
			cout << "i: " << i << endl;
			bntann1_t *p = &bns->anns[i];
			char *q = str;
			int c;
			uint32_t gi;
			scanres = fscanf(fp, "%u%8192s", &gi, str); // Ex- str[0]='c',str[1]='h',str[2]='r',str[3]='1',rest str[i's] are 0x00
			// if (scanres != 2) goto badread;

			str_enc.resize(8193,0);
			// So, strlen(str)=8192
			p->gi = gi;
			// p->gi_enc = encrypt_plaintext_integer_to_ciphertext(gi);
			
			// outfile_p_gi << gi << endl;

			// After below, str_enc gets all till 0s
			for (int h = 0; h < 8193; h++) // copy from str to str_enc till str[i]!=0
			{
				// if (str[h] == 0)
				// {
				// 	str_enc[h]=encrypt_plaintext_integer_to_ciphertext(0);
				// 	cout << "Zero encountered h: " << h << "; str[h]: " << int(str[h]) <<"; decrypt_ciphertext_to_plaintext_vector(str_enc[h])[0]: "<<(char)(decrypt_ciphertext_to_plaintext_vector(str_enc[h])[0])<< endl;

				// 	break;
				// }

				// if (h + 1 > str_enc.size())
				// 	str_enc.push_back(encrypt_plaintext_integer_to_ciphertext(str[h]));
				// else
				str_enc[h] = encrypt_plaintext_integer_to_ciphertext(str[h]);

				// cout << "h: " << h << "; str[h]: " << int(str[h]) << "; decrypt_ciphertext_to_plaintext_vector(str_enc[h])[0]: " << decrypt_ciphertext_to_plaintext_vector(str_enc[h])[0] << endl;

				if(str[h]==0)
				{
					break;
				}

				// if(h==sizeof(str)-1)
				// {
				// 	if (h + 2 > str_enc.size())
				// 		str_enc.push_back(0);
				// 	else
				// 		str_enc[h+1] = 0;

				// 	cout << "Last h with 0, decrypt_ciphertext_to_plaintext_vector(str_enc[h+1])[0]: "<<decrypt_ciphertext_to_plaintext_vector(str_enc[h+1])[0]<< endl;

				// }
			}

			// cout << "3 p->name_enc.size(): " << p->name_enc.size() << endl;

			if (scanres != 2)
				goto badread;
			p->name = strdup(str);

			// vector<int64_t> vec;
			// vec.clear();
			// for(int kk=0;kk<8193;kk++)
			// {
			// 	// cout<<"str["<<kk<<"]= "<<str[kk]<<endl;
			// 	vec.push_back((int64_t)(str[kk]));
			// }
			// Plaintext pt = cc->MakePackedPlaintext(vec);
			// auto ct = cc->Encrypt(kp.publicKey, pt);

			// if (!Serial::SerializeToFile("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/packed_cts/p_name/"+to_string(i), ct,
            //                         SerType::BINARY)) {
			// 	std::cerr << "Error writing serialization of the crypto context to "
			// 				"p_name/to_string(i)"
			// 			<< std::endl;
			// 	return NULL;
			// }

			// printf("p->name: %s; strlen(p->name): %d\n", p->name, strlen(p->name));

			// printf("convert_ciphertext_vector_to_plaintext_string(str_enc): %s\n", convert_ciphertext_vector_to_plaintext_string(str_enc));
			// cout << "str_enc.size(): " << str_enc.size() << "; strlen_string_enc(str_enc): " << strlen_string_enc(str_enc) <<"; strlen_enc(str_enc): "<<strlen_enc(str_enc) << endl;
			// p->name_enc=vecCT(1,0);
			// printf("&p->name_enc: %p\n",&p->name_enc);
			// cout<<"4 p->name_enc.size(): "<<p->name_enc.size()<<endl; //<<"; strlen_enc(p->name_enc): "<<strlen_enc(p->name_enc)<<endl;

			vecCT name_enc;
			name_enc.resize(0);
			strdup_enc(str_enc, name_enc);						
			// cout << "name_enc.size(): " << name_enc.size() <<"; str_enc.size(): "<<str_enc.size()<<endl; //<<"; strlen_enc(p->name_enc): "<<strlen_enc(p->name_enc)<<endl;
			// printf("convert_ciphertext_vector_to_plaintext_string(name_enc): %s\n", convert_ciphertext_vector_to_plaintext_string(name_enc));

			// for (int kk = 0; kk < name_enc.size(); kk++)
			// {
			// 	if (!Serial::SerializeToFile("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/serializer_folder/p_name/" + to_string(i) + "/" + to_string(kk) + ".txt", str_enc[kk], SerType::BINARY))
			// 	{
			// 		std::cerr
			// 			<< "Error writing serialization of ciphertext_bns_anns_" + to_string(i)
			// 			<< std::endl;
			// 		return;
			// 	}
			// }
			
			// outfile_p_str << str << endl;

			// p->name_enc.resize(str_enc.size()); // p->name_enc.shrink_to_fit();

			// for(int i=0;i<str_enc.size();i++)
			// 	p->name_enc[i]=str_enc[i];

			// read fasta comments
			// while (q - str < sizeof(str) - 1 && (c = fgetc(fp)) != '\n' && c != EOF) *q++ = c;

			int i_comment = 0;
			CT c_enc;
			// extra -1 because str_enc.size()=strlen(str)+1, due to extra enc('\n')
			while (i_comment < (str_enc).size() - 1 && compare_enc((c_enc = encrypt_plaintext_integer_to_ciphertext(fgetc(fp))), '\n') == 0 && compare_enc(c_enc, EOF) == 0)
			{
				// cout << "i_comment: " << i_comment << endl;
				// printf("c: %d; %c\n", decrypt_ciphertext_to_plaintext_vector(c_enc)[0], decrypt_ciphertext_to_plaintext_vector(c_enc)[0]);
				str_enc[i_comment++] = c_enc;
			}

			// cout << "Done1" << endl;
			// printf("out decrypt_ciphertext_to_plaintext_vector(c_enc)[0]: %d\n", decrypt_ciphertext_to_plaintext_vector(c_enc)[0]);

			// while (c != '\n' && c != EOF) c = fgetc(fp);
			int ind2 = 0;
			while (compare_enc(c_enc, '\n') == 0 && compare_enc(c_enc, EOF) == 0)
			{
				// cout << "ind2: " << ind2 << endl;
				// cout << "In loop" << endl;
				c_enc = encrypt_plaintext_integer_to_ciphertext(fgetc(fp));
				// printf("decrypt_ciphertext_to_plaintext_vector(c_enc)[0]: %d\n", decrypt_ciphertext_to_plaintext_vector(c_enc)[0]);
				ind2++;
			}
			if (compare_enc(c_enc, EOF) == 1)
			{ // returns 1 when c_enc and EOF are equal
				scanres = EOF;
				// printf("Bad decrypt_ciphertext_to_plaintext_vector(c_enc)[0]: %d\n", decrypt_ciphertext_to_plaintext_vector(c_enc)[0]);
				goto badread;
			}
			// *q = 0;
			str_enc[i_comment] = encrypt_plaintext_integer_to_ciphertext('\0');
			// for (int h = 0; h < i_comment; h++)
			// 	str[h] = decrypt_ciphertext_to_plaintext_vector(str_enc[h])[0];
			// str[i_comment] = 0;

			for (int h = 0; h < i_comment; h++)
				str[h] = decrypt_ciphertext_to_plaintext_vector(str_enc[h])[0];
			str[i_comment] = 0;

			// cout << "strlen(new str): " << strlen(str) << endl;
			// printf("new str: %s\n", str);

			// printf("new str_enc: %s\n", convert_ciphertext_vector_to_plaintext_string(str_enc));
			// cout << "strlen_enc(new str_enc): " << strlen_enc(str_enc) << endl;

			// cout << "strlen_string_enc(new str_enc): " << strlen_string_enc(str_enc) << endl;

			// assert(strlen(str) < 8192);
			assert(strlen_string_enc(str_enc) < 8192);
			vecCT p_anno_enc;

			// if (q - str > 1 && strcmp(str, " (null)") != 0) p->anno = strdup(str + 1); // skip leading space
			// else p->anno = strdup("");
			// vecCT str_enc_2; str_enc_2.assign(str_enc.begin()+1,str_enc.end());
			// if(i_comment>1 &&  strcmp_enc(str_enc," (null)")!=0 ) strdup_enc(str_enc_2,p->anno_enc);
			if (i_comment > 1 && strcmp_enc(str_enc, " (null)") != 0)
				assign_string_to_vecCT(p_anno_enc, str + 1, -1);
			else
				assign_string_to_vecCT(p_anno_enc, "", -1);

			// printf("p_anno_enc: %s\n", convert_ciphertext_vector_to_plaintext_string(p_anno_enc));

			// for (int kk = 0; kk < p_anno_enc.size(); kk++)
			// {
			// 	if (!Serial::SerializeToFile("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/serializer_folder/p_anno_enc/" + to_string(i) + "/" + to_string(kk) + ".txt", p_anno_enc[kk], SerType::BINARY))
			// 	{
			// 		std::cerr
			// 			<< "Error writing serialization of ciphertext_bns_anns_" + to_string(i)
			// 			<< std::endl;
			// 		return;
			// 	}
			// }
			
			// outfile_p_anno<<convert_ciphertext_vector_to_plaintext_string(p_anno_enc)<<endl;

			// cout<<"Earlier p->anno_enc.size(): "<<p->anno_enc.size()<<endl;
			// p->anno_enc.resize(strlen_enc(p->anno_enc)+1);
			// cout<<"Later p->anno_enc.size(): "<<p->anno_enc.size()<<endl;

			// else assign_string_to_vecCT(p->anno_enc,"",-1);
			int32_t p_len, p_n_ambs;
			// cout << "Reading rest of broad.ann" << endl;
			// read the rest
			scanres = fscanf(fp, "%lld%d%d", &xx, &p_len, &p_n_ambs);
			// cout << "Read broad.ann" << endl;
			// xx_enc = encrypt_plaintext_integer_to_ciphertext(xx);

			// if (!Serial::SerializeToFile("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/serializer_folder/p_len/" + to_string(i) + ".txt", encrypt_plaintext_integer_to_ciphertext(p_len), SerType::BINARY))
			// {
			// 	std::cerr
			// 		<< "Error writing serialization of ciphertext_bns_anns_" + to_string(i)
			// 		<< std::endl;
			// 	return;
			// }
			// if (!Serial::SerializeToFile("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/serializer_folder/p_n_ambs/" + to_string(i) + ".txt", encrypt_plaintext_integer_to_ciphertext(p_n_ambs), SerType::BINARY))
			// {
			// 	std::cerr
			// 		<< "Error writing serialization of ciphertext_bns_anns_" + to_string(i)
			// 		<< std::endl;
			// 	return;
			// }

			// outfile_p_len<<p_len<<endl;
			// outfile_p_n_ambs<<p_n_ambs<<endl;

			// p->n_ambs_enc = encrypt_plaintext_integer_to_ciphertext(p_n_ambs);
			if (scanres != 3)
				goto badread;
			p->offset = xx;

			// if (!Serial::SerializeToFile("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/serializer_folder/p_offset/" + to_string(i) + ".txt", encrypt_plaintext_integer_to_ciphertext(xx), SerType::BINARY))
			// {
			// 	std::cerr
			// 		<< "Error writing serialization of ciphertext_bns_anns_" + to_string(i)
			// 		<< std::endl;
			// 	return;
			// }
			
			// outfile_p_offset<<xx<<endl;

			// p->offset_enc = xx_enc;

			// if (!Serial::SerializeToFile("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/serializer_folder/ciphertext_bns_anns_" + to_string(i) + ".txt", c1, SerType::BINARY))
			// {
			// 	std::cerr
			// 		<< "Error writing serialization of ciphertext_bns_anns_" + to_string(i)
			// 		<< std::endl;
			// 	return;
			// }
		}
		err_fclose(fp);
	}

	// { // read .ann
	// 	fp = xopen(fname = ann_filename, "r");
	// 	scanres = fscanf(fp, "%lld%d%u", &xx, &bns->n_seqs, &bns->seed);
	// 	assert(bns->n_seqs >= 0 && bns->n_seqs <= INT_MAX);
	// 	if (scanres != 3)
	// 		goto badread;
	// 	bns->l_pac = xx;
	// 	bns->l_pac_enc = encrypt_plaintext_integer_to_ciphertext(xx);
	// 	bns->anns = (bntann1_t *)calloc(bns->n_seqs, sizeof(bntann1_t));
	// 	for (int kk = 0; kk < bns->n_seqs; kk++)
	// 	{
	// 		bns->anns[kk].is_alt_enc = encrypt_plaintext_integer_to_ciphertext(0);
	// 	}

	// 	assert(bns->anns != NULL);
	// 	for (i = 0; i < bns->n_seqs; ++i)
	// 	{
	// 		// printf("-----------------------------------------\n");
	// 		// printf("i: %d\n",i);

	// 		bntann1_t *p = bns->anns + i;
	// 		char *q = str;
	// 		int c;
	// 		// read gi and sequence name
	// 		scanres = fscanf(fp, "%u%8192s", &p->gi, str);
	// 		if (scanres != 2)
	// 			goto badread;
	// 		p->name = strdup(str);

	// 		// vector<int64_t> vec;
	// 		// vec.clear();
	// 		// for(int kk=0;kk<8193;kk++)
	// 		// {
	// 		// 	vec.push_back((int64_t)(str[kk]));
	// 		// }
	// 		// Plaintext pt = cc->MakePackedPlaintext(vec);
	// 		// auto ct = cc->Encrypt(kp.publicKey, pt);

	// 		// if (!Serial::SerializeToFile("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/serializer_folder/p_name/"+to_string(i), ct,
    //         //                         SerType::BINARY)) {
	// 		// 	std::cerr << "Error writing serialization of the crypto context to "
	// 		// 				"p_name/to_string(i)"
	// 		// 			<< std::endl;
	// 		// 	return;
	// 		// }

	// 		// read fasta comments

	// 		// printf("sizeof(str): %lu\n",sizeof(str)); int ind1=0; int ind2=0;
	// 		// printf("strlen(str): %d\n",strlen(str));
	// 		// printf("str: %s\n",str);

	// 		while (q - str < sizeof(str) - 1 && (c = fgetc(fp)) != '\n' && c != EOF)
	// 		{
	// 			// printf("ind1: %d\n",ind1);
	// 			// ind1++;
	// 			*q++ = c;
	// 			// printf("c: %d; %c\n",c,c);
	// 			// cout<<"c: "<<c<<"; "<<char(c)<<endl;
	// 		}
	// 		// printf("new str: %s\n",str);

	// 		while (c != '\n' && c != EOF)
	// 		{
	// 			// printf("ind2: %d\n",ind2);
	// 			c = fgetc(fp);
	// 			// printf("c: %d; %c\n",c,c);
	// 			// ind2++;
	// 		}

	// 		if (c == EOF)
	// 		{
	// 			scanres = EOF;
	// 			goto badread;
	// 		}
	// 		*q = 0;
	// 		assert(strlen(str) < 8192);
	// 		if (q - str > 1 && strcmp(str, " (null)") != 0)
	// 			p->anno = strdup(str + 1); // skip leading space
	// 		else
	// 			p->anno = strdup("");

	// 		// printf("p->anno: %s\n",p->anno);
	// 		// read the rest
	// 		scanres = fscanf(fp, "%lld%d%d", &xx, &p->len, &p->n_ambs);
	// 		if (scanres != 3)
	// 			goto badread;
	// 		p->offset = xx;
	// 		// printf("offset: %d\n",xx);
	// 	}
	// 	err_fclose(fp);
	// }
	{ // read .amb
		int64_t l_pac;
		CT l_pac_enc;
		int32_t n_seqs;
		CT n_seqs_enc;
		int32_t n_holes;
		fp = xopen(fname = amb_filename, "r");
		scanres = fscanf(fp, "%lld%d%d", &xx, &n_seqs, &n_holes);
		bns->n_holes = n_holes;
		xx_enc = encrypt_plaintext_integer_to_ciphertext(xx);
		n_seqs_enc = encrypt_plaintext_integer_to_ciphertext(n_seqs);
		bns->n_holes_enc = encrypt_plaintext_integer_to_ciphertext(n_holes);
		bns->n_seqs_enc = encrypt_plaintext_integer_to_ciphertext(n_seqs);
		// printf("decrypt_ciphertext_to_plaintext_vector(bns->n_holes_enc)[0]: %d; sizeof(bntamb1_t): %zu\n",decrypt_ciphertext_to_plaintext_vector(bns->n_holes_enc)[0], sizeof(bntamb1_t));

		// assert(bns->n_holes >= 0 && bns->n_holes <= INT_MAX);
		assert(decrypt_ciphertext_to_plaintext_vector(bns->n_holes_enc)[0] >= 0 && operate_and_decrypt(bns->n_holes_enc, "-", INT_MAX) <= 0);
		if (scanres != 3)
			goto badread;
		// l_pac = xx;
		l_pac_enc = encrypt_plaintext_integer_to_ciphertext(xx);
		// xassert(l_pac == bns->l_pac && n_seqs == bns->n_seqs, "inconsistent .ann and .amb files.");
		xassert(compare_enc(l_pac_enc, bns->l_pac_enc) && compare_enc(n_seqs_enc, bns->n_seqs_enc), "inconsistent .ann and .amb files.");
		if (compare_enc(bns->n_holes_enc, 0) == 0)
		{
			// printf("In\n");
			// bns->ambs = (bntamb1_t*)calloc(bns->n_holes, sizeof(bntamb1_t));
			bns->ambs = (bntamb1_t *)calloc(decrypt_ciphertext_to_plaintext_vector(bns->n_holes_enc)[0], sizeof(bntamb1_t));
			assert(bns->ambs != NULL);
		}
		else
		{
			bns->ambs = 0;
		}

		// cout << "n_holes: " << n_holes << endl;

		ofstream outfile_ambs_len("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/ambs_len.txt");
		ofstream outfile_ambs_offset("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/ambs_offset.txt");
		ofstream outfile_ambs_amb("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/ambs_amb.txt");
		// cout << "Before" << endl;

		for (i = 0; operate_and_decrypt(bns->n_holes_enc, "-", i) > 0; ++i)
		{
			bntamb1_t *p = &bns->ambs[i];
			
			scanres = fscanf(fp, "%lld%d%8192s", &xx, &p->len, str);
			// printf("\ni: %d\n",i);
			// printf("p->len: %d\n", p->len);
			// printf("xx: %lld\n", xx);
			xx_enc = encrypt_plaintext_integer_to_ciphertext((int64_t)(xx));
			// printf("decrypt_ciphertext_to_plaintext_vector(xx_enc)[0]: %lld\n", decrypt_ciphertext_to_plaintext_vector(xx_enc)[0]);
			// printf("str: %s\n",str);

			// str[8192]=0;
			// for(int q=0;q<8193;q++)
			// {
			// 	printf("str[%d]=%c\n",q,str[q]);
			// }

			// printf("i: %d, xx: %lld, p_len: %d, str: %s\n",i,xx,p_len,str);

			// if (!Serial::SerializeToFile("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/serializer_folder/ambs_len/" + to_string(i) + ".txt", encrypt_plaintext_integer_to_ciphertext(p_len), SerType::BINARY))
			// {
			// 	std::cerr
			// 		<< "Error writing serialization of ciphertext_bns_anns_" + to_string(i)
			// 		<< std::endl;
			// 	return;
			// }

			outfile_ambs_len << p->len << endl;
			// printf("Wrote\n");
			str_enc.resize(8193, 0);

			// p->len_enc = encrypt_plaintext_integer_to_ciphertext(p_len);
			for (int h = 0; h < 8193; h++)
			{
				str_enc[h] = encrypt_plaintext_integer_to_ciphertext(str[h]);
				// printf("decrypt_ciphertext_to_plaintext_vector(str_enc[%d])[0]: %d\n", h, decrypt_ciphertext_to_plaintext_vector(str_enc[h])[0]);
				if (str[h] == 0)
				{
					break;
				}
			}

			// printf("convert_ciphertext_vector_to_plaintext_string(str_enc): %s\n",convert_ciphertext_vector_to_plaintext_string(str_enc));

			if (scanres != 3)
				goto badread;

			// p->offset = xx;
			// p->amb = str[0];

			// if (!Serial::SerializeToFile("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/serializer_folder/ambs_offset/" + to_string(i) + ".txt", xx_enc, SerType::BINARY))
			// {
			// 	std::cerr
			// 		<< "Error writing serialization of ciphertext_bns_anns_" + to_string(i)
			// 		<< std::endl;
			// 	return;
			// }
			// if (!Serial::SerializeToFile("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/serializer_folder/ambs_amb/" + to_string(i) + ".txt", str_enc[0], SerType::BINARY))
			// {
			// 	std::cerr
			// 		<< "Error writing serialization of ciphertext_bns_anns_" + to_string(i)
			// 		<< std::endl;
			// 	return;
			// }

			// outfile_ambs_offset << xx << endl;
			// outfile_ambs_amb << str[0] << endl;

			// p->offset_enc = xx_enc;
			// p->amb_enc = str_enc[0];
		}
		err_fclose(fp);
	}
	{ // open .pac
		bns->fp_pac = xopen(pac_filename, "rb");
	}
	return bns;

badread:
	if (EOF == scanres)
	{
		err_fatal(__func__, "Error reading %s : %s\n", fname, ferror(fp) ? strerror(errno) : "Unexpected end of file");
	}
	err_fatal(__func__, "Parse error reading %s\n", fname);
}

// bntseq_t *bns_restore_core(const char *ann_filename, const char *amb_filename, const char *pac_filename)
// {
// 	char str[8193];
// 	FILE *fp;
// 	const char *fname;
// 	bntseq_t *bns;
// 	long long xx;
// 	int i;
// 	int scanres;
// 	bns = (bntseq_t *)calloc(1, sizeof(bntseq_t));
// 	assert(bns != 0);
// 	{ // read .ann
// 		fp = xopen(fname = ann_filename, "r");
// 		scanres = fscanf(fp, "%lld%d%u", &xx, &bns->n_seqs, &bns->seed);
// 		assert(bns->n_seqs >= 0 && bns->n_seqs <= INT_MAX);
// 		if (scanres != 3)
// 			goto badread;
// 		bns->l_pac = xx;
// 		bns->l_pac_enc = encrypt_plaintext_integer_to_ciphertext(xx);
// 		bns->anns = (bntann1_t *)calloc(bns->n_seqs, sizeof(bntann1_t));
// 		for (int kk = 0; kk < bns->n_seqs; kk++)
// 		{
// 			bns->anns[kk].is_alt_enc = encrypt_plaintext_integer_to_ciphertext(0);
// 		}

// 		assert(bns->anns != NULL);
// 		for (i = 0; i < bns->n_seqs; ++i)
// 		{
// 			// printf("-----------------------------------------\n");
// 			// printf("i: %d\n",i);

// 			bntann1_t *p = bns->anns + i;
// 			char *q = str;
// 			int c;
// 			// read gi and sequence name
// 			scanres = fscanf(fp, "%u%8192s", &p->gi, str);
// 			if (scanres != 2)
// 				goto badread;
// 			p->name = strdup(str);

// 			// read fasta comments

// 			// printf("sizeof(str): %lu\n",sizeof(str)); int ind1=0; int ind2=0;
// 			// printf("strlen(str): %d\n",strlen(str));
// 			// printf("str: %s\n",str);

// 			while (q - str < sizeof(str) - 1 && (c = fgetc(fp)) != '\n' && c != EOF)
// 			{
// 				// printf("ind1: %d\n",ind1);
// 				// ind1++;
// 				*q++ = c;
// 				// printf("c: %d; %c\n",c,c);
// 				// cout<<"c: "<<c<<"; "<<char(c)<<endl;
// 			}
// 			// printf("new str: %s\n",str);

// 			while (c != '\n' && c != EOF)
// 			{
// 				// printf("ind2: %d\n",ind2);
// 				c = fgetc(fp);
// 				// printf("c: %d; %c\n",c,c);
// 				// ind2++;
// 			}

// 			if (c == EOF)
// 			{
// 				scanres = EOF;
// 				goto badread;
// 			}
// 			*q = 0;
// 			assert(strlen(str) < 8192);
// 			if (q - str > 1 && strcmp(str, " (null)") != 0)
// 				p->anno = strdup(str + 1); // skip leading space
// 			else
// 				p->anno = strdup("");

// 			// printf("p->anno: %s\n",p->anno);
// 			// read the rest
// 			scanres = fscanf(fp, "%lld%d%d", &xx, &p->len, &p->n_ambs);
// 			if (scanres != 3)
// 				goto badread;
// 			p->offset = xx;
// 			// printf("offset: %d\n",xx);
// 		}
// 		err_fclose(fp);
// 	}
// 	{ // read .amb
// 		int64_t l_pac;
// 		int32_t n_seqs;
// 		fp = xopen(fname = amb_filename, "r");
// 		scanres = fscanf(fp, "%lld%d%d", &xx, &n_seqs, &bns->n_holes);
// 		assert(bns->n_holes >= 0 && bns->n_holes <= INT_MAX);
// 		if (scanres != 3)
// 			goto badread;
// 		l_pac = xx;
// 		xassert(l_pac == bns->l_pac && n_seqs == bns->n_seqs, "inconsistent .ann and .amb files.");
// 		if (bns->n_holes)
// 		{
// 			bns->ambs = (bntamb1_t *)calloc(bns->n_holes, sizeof(bntamb1_t));
// 			assert(bns->ambs != NULL);
// 		}
// 		else
// 		{
// 			bns->ambs = 0;
// 		}
// 		for (i = 0; i < bns->n_holes; ++i)
// 		{
// 			bntamb1_t *p = bns->ambs + i;
// 			scanres = fscanf(fp, "%lld%d%8192s", &xx, &p->len, str);
// 			if (scanres != 3)
// 				goto badread;
// 			p->offset = xx;
// 			p->amb = str[0];
// 		}
// 		err_fclose(fp);
// 	}
// 	{ // open .pac
// 		bns->fp_pac = xopen(pac_filename, "rb");
// 	}
// 	return bns;

// badread:
// 	if (EOF == scanres)
// 	{
// 		err_fatal(__func__, "Error reading %s : %s\n", fname, ferror(fp) ? strerror(errno) : "Unexpected end of file");
// 	}
// 	err_fatal(__func__, "Parse error reading %s\n", fname);
// }

// bntseq_t *bns_restore_core(const char *ann_filename, const char* amb_filename, const char* pac_filename)
// {
// 	char str[8193];
// 	FILE *fp;
// 	const char *fname;
// 	bntseq_t *bns;
// 	long long xx;
// 	int i;
// 	int scanres;
// 	bns = (bntseq_t*)calloc(1, sizeof(bntseq_t));
// 	assert(bns != 0);
// 	{ // read .ann
// 		fp = xopen(fname = ann_filename, "r");
// 		scanres = fscanf(fp, "%lld%d%u", &xx, &bns->n_seqs, &bns->seed);
// 		assert(bns->n_seqs >= 0 && bns->n_seqs <= INT_MAX);
// 		if (scanres != 3) goto badread;
// 		bns->l_pac = xx;
// 		bns->anns = (bntann1_t*)calloc(bns->n_seqs, sizeof(bntann1_t));
//         assert(bns->anns != NULL);
// 		for (i = 0; i < bns->n_seqs; ++i) {
// 			printf("-----------------------------------------\n");
// 			printf("i: %d\n",i);

// 			bntann1_t *p = bns->anns + i;
// 			char *q = str;
// 			int c;
// 			// read gi and sequence name
// 			scanres = fscanf(fp, "%u%8192s", &p->gi, str);
// 			if (scanres != 2) goto badread;
// 			p->name = strdup(str);
// 			// read fasta comments
// 			printf("sizeof(str): %lu\n",sizeof(str)); int ind1=0; int ind2=0;
// 			printf("strlen(str): %d\n",strlen(str));
// 			printf("str: %s\n",str);

// 			while (q - str < sizeof(str) - 1 && (c = fgetc(fp)) != '\n' && c != EOF)
// 			{
// 				printf("ind1: %d\n",ind1);
// 				ind1++;
// 				*q++ = c;
// 				printf("c: %d; %c\n",c,c);
// 				// cout<<"c: "<<c<<"; "<<char(c)<<endl;
// 			}
// 			printf("new str: %s\n",str);

// 			while (c != '\n' && c != EOF)
// 			{
// 				printf("ind2: %d\n",ind2);
// 				c = fgetc(fp);
// 				printf("c: %d; %c\n",c,c);
// 				ind2++;
// 			}

// 			if (c == EOF) {
// 				scanres = EOF;
// 				goto badread;
// 			}
// 			*q = 0;
// 			assert(strlen(str) < 8192);
// 			if (q - str > 1 && strcmp(str, " (null)") != 0) p->anno = strdup(str + 1); // skip leading space
// 			else p->anno = strdup("");

// 			printf("p->anno: %s\n",p->anno);
// 			// read the rest
// 			scanres = fscanf(fp, "%lld%d%d", &xx, &p->len, &p->n_ambs);
// 			if (scanres != 3) goto badread;
// 			p->offset = xx;
// 			printf("offset: %d\n",xx);

// 		}
// 		err_fclose(fp);
// 	}
// 	{ // read .amb
// 		int64_t l_pac;
// 		int32_t n_seqs;
// 		fp = xopen(fname = amb_filename, "r");
// 		scanres = fscanf(fp, "%lld%d%d", &xx, &n_seqs, &bns->n_holes);
// 		assert(bns->n_holes >= 0 && bns->n_holes <= INT_MAX);
// 		if (scanres != 3) goto badread;
// 		l_pac = xx;
// 		xassert(l_pac == bns->l_pac && n_seqs == bns->n_seqs, "inconsistent .ann and .amb files.");
//         if(bns->n_holes){
//             bns->ambs = (bntamb1_t*)calloc(bns->n_holes, sizeof(bntamb1_t));
//             assert(bns->ambs != NULL);
//         }
//         else{
//             bns->ambs = 0;
//         }
// 		for (i = 0; i < bns->n_holes; ++i) {
// 			bntamb1_t *p = bns->ambs + i;
// 			scanres = fscanf(fp, "%lld%d%8192s", &xx, &p->len, str);
// 			if (scanres != 3) goto badread;
// 			p->offset = xx;
// 			p->amb = str[0];
// 		}
// 		err_fclose(fp);
// 	}
// 	{ // open .pac
// 		bns->fp_pac = xopen(pac_filename, "rb");
// 	}
// 	return bns;

//  badread:
// 	if (EOF == scanres) {
// 		err_fatal(__func__, "Error reading %s : %s\n", fname, ferror(fp) ? strerror(errno) : "Unexpected end of file");
// 	}
// 	err_fatal(__func__, "Parse error reading %s\n", fname);
// }

// idx->bns = bns_restore(prefix);

bntseq_t *bns_restore(const char *prefix) // prefix: ../../input-datasets/fmi/broad
{
	auto bns_restore_time_start_1 = high_resolution_clock::now();
	printf("In bns_restore, prefix: %s\n", prefix);
	char ann_filename[PATH_MAX], amb_filename[PATH_MAX], pac_filename[PATH_MAX], alt_filename[PATH_MAX];
	FILE *fp;
	bntseq_t *bns;
	// assert(strlen(prefix) + 4 < 1024);
	strcpy_s(ann_filename, PATH_MAX, prefix);
	strcat_s(ann_filename, PATH_MAX, ".ann");
	strcpy_s(amb_filename, PATH_MAX, prefix);
	strcat_s(amb_filename, PATH_MAX, ".amb");
	strcpy_s(pac_filename, PATH_MAX, prefix);
	strcat_s(pac_filename, PATH_MAX, ".pac");

	auto bns_restore_time_end_1 = high_resolution_clock::now();

	printf("Going into bns_restore_core()\n");

	auto bns_restore_core_time_start = high_resolution_clock::now();
	bns = bns_restore_core(ann_filename, amb_filename, pac_filename);
	auto bns_restore_core_time_end = high_resolution_clock::now();
    auto duration_bns_restore_core = duration_cast<microseconds>(bns_restore_core_time_end - bns_restore_core_time_start);
    cout << "Time taken by bns_restore_core: "<<duration_bns_restore_core.count() <<" microseconds" << endl;

	auto bns_restore_time_start_2 = high_resolution_clock::now();

	bns->n_seqs_enc = encrypt_plaintext_integer_to_ciphertext(bns->n_seqs);

	if (bns == 0)
		return 0;
	strcpy_s(alt_filename, PATH_MAX, prefix);
	strcat_s(alt_filename, PATH_MAX, ".alt");
	// printf("alt_filename: %s\n", alt_filename);
	if ((fp = fopen(alt_filename, "r")) != 0)
	{					// read .alt file if present
		char str[1024]; // vecCT str_enc(1024,0);
		khash_t(str) * h;
		int c, i, absent;
		CT c_enc;
		khint_t k;
		h = kh_init(str);
		assert(h != NULL);
		for (i = 0; operate_and_decrypt(bns->n_seqs_enc, "-", i) > 0; ++i)
		{ 
			int32_t n_seqs;
			CT ct;
			if (Serial::DeserializeFromFile("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/packed_cts/p_name/"+to_string(i), ct, SerType::BINARY) == false) {
				std::cerr << "Could not read the ciphertext packed_cts/p_name/"+to_string(i) << std::endl;
				return NULL;
			}
			vector<int64_t> vv=decrypt_ciphertext_to_plaintext_vector(ct);
			string ss="";
			for(int kk=0;kk<vv.size();kk++)
			{
				// cout<<"vv["<<kk<<"]="<<vv[kk]<<endl;
				ss+=vv[kk];
				if(vv[kk]==0)
					break;
			}	
			const int length = ss.length();
  
			// declaring character array (+1 for null terminator)
			char* qq = new char[length + 1];
			strcpy(qq, ss.c_str());

			// printf("i: %d, qq: %s; strlen(qq): %d; bns->anns[i].name: %s, strlen(bns->anns[i].name): %d\n",i,qq,strlen(qq),bns->anns[i].name,strlen(bns->anns[i].name));
			strcpy(bns->anns[i].name,qq);
			// printf("qq: %s; strlen(qq): %d; bns->anns[i].name: %s, strlen(bns->anns[i].name): %d\n",qq,strlen(qq),bns->anns[i].name,strlen(bns->anns[i].name));

			// k = kh_put(str, h, bns->anns[i].name, &absent);
			k = kh_put(str, h, qq, &absent);
			kh_val(h, k) = i;
		}
		i = 0;
		while (compare_enc((c_enc = encrypt_plaintext_integer_to_ciphertext(fgetc(fp))), EOF) != 1)
		{
			// printf("c: %c; i: %d\n", decrypt_ciphertext_to_plaintext_vector(c_enc)[0], i);
			if (compare_enc(c_enc, '\t') || compare_enc(c_enc, '\n') || compare_enc(c_enc, '\r'))
			{
				str[i] = 0;
				if (str[0] != '@')
				{
					// printf("kh_val(h, k): %d\n", kh_val(h, k));
					k = kh_get(str, h, str);
					if (k != kh_end(h))
					{
						bns->anns[kh_val(h, k)].is_alt = 1;
						bns->anns[kh_val(h, k)].is_alt_enc = encrypt_plaintext_integer_to_ciphertext(1);
					}
				}
				while (compare_enc(c_enc, '\n') == 0 && compare_enc(c_enc, EOF) == 0)
					c_enc = encrypt_plaintext_integer_to_ciphertext(fgetc(fp));
				i = 0;
			}
			else
			{
				// printf("Decrypting c_enc\n");
				// printf("decrypt_ciphertext_to_plaintext_vector(c_enc)[0]: %c\n",decrypt_ciphertext_to_plaintext_vector(c_enc)[0]);
				str[i++] = decrypt_ciphertext_to_plaintext_vector(c_enc)[0]; // FIXME: potential segfault here
			}
		}
		printf("Out\n");
		
		// while ((c = fgetc(fp)) != EOF)
		// {
		// 	// printf("c: %c; i: %d\n",c,i);
		// 	if (c == '\t' || c == '\n' || c == '\r')
		// 	{
		// 		str[i] = 0;
		// 		if (str[0] != '@')
		// 		{
		// 			// printf("kh_val(h, k): %d\n",kh_val(h, k));
		// 			k = kh_get(str, h, str);
		// 			if (k != kh_end(h))
		// 			{
		// 				// printf("In=1\n");
		// 				// bns->anns[kh_val(h, k)].is_alt = 1;
		// 				bns->anns[kh_val(h, k)].is_alt_enc = encrypt_plaintext_integer_to_ciphertext(1);
		// 			}
		// 		}
		// 		while (c != '\n' && c != EOF)
		// 			c = fgetc(fp);
		// 		i = 0;
		// 	}
		// 	else
		// 		str[i++] = c; // FIXME: potential segfault here
		// }
		// printf("Out\n");

		kh_destroy(str, h);
		fclose(fp);
	}

	auto bns_restore_time_end_2 = high_resolution_clock::now();
	auto duration_bns_restore = duration_cast<microseconds>(bns_restore_time_end_1-bns_restore_time_start_1+bns_restore_time_end_2-bns_restore_time_start_2);
    cout << "Time taken by bns_restore: "<<duration_bns_restore.count()<< " microseconds" << endl;

	return bns;
}

// bntseq_t *bns_restore(const char *prefix)
// {
// 	printf("prefix: %s\n",prefix);
// 	char ann_filename[PATH_MAX], amb_filename[PATH_MAX], pac_filename[PATH_MAX], alt_filename[PATH_MAX];
// 	FILE *fp;
// 	bntseq_t *bns;
// 	//assert(strlen(prefix) + 4 < 1024);
// 	strcpy_s(ann_filename, PATH_MAX, prefix); strcat_s(ann_filename, PATH_MAX, ".ann");
// 	strcpy_s(amb_filename, PATH_MAX, prefix); strcat_s(amb_filename, PATH_MAX, ".amb");
// 	strcpy_s(pac_filename, PATH_MAX, prefix); strcat_s(pac_filename, PATH_MAX, ".pac");
// 	bns = bns_restore_core(ann_filename, amb_filename, pac_filename);
// 	if (bns == 0) return 0;
//     strcpy_s(alt_filename, PATH_MAX, prefix); strcat_s(alt_filename, PATH_MAX, ".alt");
// 	if ((fp = fopen(alt_filename, "r")) != 0) { // read .alt file if present
// 		char str[1024];
// 		khash_t(str) *h;									// #define khash_t(name) kh_##name##_t

// 		// typedef unsigned int khint32_t;
// 		// typedef khint32_t khint_t;
// 		// typedef khint_t khiter_t;

// 		int c, i, absent;
// 		khint_t k;
// 		h = kh_init(str);
//         assert(h != NULL);
// 		for (i = 0; i < bns->n_seqs; ++i) {
// 			k = kh_put(str, h, bns->anns[i].name, &absent);
// 			kh_val(h, k) = i;
// 		}
// 		i = 0;
// 		while ((c = fgetc(fp)) != EOF) {
// 			printf("c: %c; i: %d\n",c,i);
// 			if (c == '\t' || c == '\n' || c == '\r') {
// 				str[i] = 0;
// 				if (str[0] != '@') {
// 					printf("kh_val(h, k): %d\n",kh_val(h, k));
// 					k = kh_get(str, h, str);
// 					if (k != kh_end(h))
// 						bns->anns[kh_val(h, k)].is_alt = 1;
// 				}
// 				while (c != '\n' && c != EOF) c = fgetc(fp);
// 				i = 0;
// 			} else str[i++] = c; // FIXME: potential segfault here
// 		}
// 		kh_destroy(str, h);
// 		fclose(fp);
// 	}
// 	return bns;
// }

void bns_destroy(bntseq_t *bns)
{
	if (bns == 0)
		return;
	else
	{
		int i;
		if (bns->fp_pac)
			err_fclose(bns->fp_pac);
		free(bns->ambs);
		for (i = 0; i < bns->n_seqs; ++i)
		{
			free(bns->anns[i].name);
			free(bns->anns[i].anno);
		}
		free(bns->anns);
		free(bns);
	}
}

#define _set_pac(pac, l, c) ((pac)[(l) >> 2] |= (c) << ((~(l)&3) << 1))
#define _get_pac(pac, l) ((pac)[(l) >> 2] >> ((~(l)&3) << 1) & 3)

static uint8_t *add1(const kseq_t *seq, bntseq_t *bns, uint8_t *pac, int64_t *m_pac, int *m_seqs, int *m_holes, bntamb1_t **q)
{
	bntann1_t *p;
	int i, lasts;
	if (bns->n_seqs == *m_seqs)
	{
		*m_seqs <<= 1;
		bns->anns = (bntann1_t *)realloc(bns->anns, *m_seqs * sizeof(bntann1_t));
		assert(bns->anns != NULL);
	}
	p = bns->anns + bns->n_seqs;
	p->name = strdup((char *)seq->name.s);
	p->anno = seq->comment.l > 0 ? strdup((char *)seq->comment.s) : strdup("(null)");
	p->gi = 0;
	p->len = seq->seq.l;
	p->offset = (bns->n_seqs == 0) ? 0 : (p - 1)->offset + (p - 1)->len;
	p->n_ambs = 0;
	for (i = lasts = 0; i < seq->seq.l; ++i)
	{
		int c = nst_nt4_table[(int)seq->seq.s[i]];
		if (c >= 4)
		{ // N
			if (lasts == seq->seq.s[i])
			{ // contiguous N
				++(*q)->len;
			}
			else
			{
				if (bns->n_holes == *m_holes)
				{
					(*m_holes) <<= 1;
					bns->ambs = (bntamb1_t *)realloc(bns->ambs, (*m_holes) * sizeof(bntamb1_t));
				}
				*q = bns->ambs + bns->n_holes;
				(*q)->len = 1;
				(*q)->offset = p->offset + i;
				(*q)->amb = seq->seq.s[i];
				++p->n_ambs;
				++bns->n_holes;
			}
		}
		lasts = seq->seq.s[i];
		{ // fill buffer
			if (c >= 4)
				c = lrand48() & 3;
			if (bns->l_pac == *m_pac)
			{ // double the pac size
				*m_pac <<= 1;
				pac = (uint8_t *)realloc(pac, *m_pac / 4);
				memset(pac + bns->l_pac / 4, 0, (*m_pac - bns->l_pac) / 4);
			}
			_set_pac(pac, bns->l_pac, c);
			++bns->l_pac;
		}
	}
	++bns->n_seqs;
	return pac;
}

int64_t bns_fasta2bntseq(gzFile fp_fa, const char *prefix, int for_only)
{
	extern void seq_reverse(int len, ubyte_t *seq, int is_comp); // in bwaseqio.c
	kseq_t *seq;
	char name[PATH_MAX];
	bntseq_t *bns;
	uint8_t *pac = 0;
	int32_t m_seqs, m_holes;
	int64_t ret = -1, m_pac, l;
	bntamb1_t *q;
	FILE *fp;

	// initialization
	seq = kseq_init(fp_fa);
	bns = (bntseq_t *)calloc(1, sizeof(bntseq_t));
	assert(bns != NULL);
	bns->seed = 11; // fixed seed for random generator
	srand48(bns->seed);
	m_seqs = m_holes = 8;
	m_pac = 0x10000;
	bns->anns = (bntann1_t *)calloc(m_seqs, sizeof(bntann1_t));
	assert(bns->anns != NULL);
	bns->ambs = (bntamb1_t *)calloc(m_holes, sizeof(bntamb1_t));
	assert(bns->ambs != NULL);
	pac = (uint8_t *)calloc(m_pac / 4, 1);
	if (pac == NULL)
	{
		perror("Allocation of pac failed");
		exit(EXIT_FAILURE);
	}
	q = bns->ambs;
	// assert(strlen(prefix) + 4 < 1024);
	strcpy_s(name, PATH_MAX, prefix);
	strcat_s(name, PATH_MAX, ".pac");
	fp = xopen(name, "wb");
	// read sequences
	while (kseq_read(seq) >= 0)
		pac = add1(seq, bns, pac, &m_pac, &m_seqs, &m_holes, &q);
	if (!for_only)
	{ // add the reverse complemented sequence
		m_pac = (bns->l_pac * 2 + 3) / 4 * 4;
		pac = (uint8_t *)realloc(pac, m_pac / 4);
		if (pac == NULL)
		{
			perror("Reallocation of pac failed");
			exit(EXIT_FAILURE);
		}
		memset(pac + (bns->l_pac + 3) / 4, 0, (m_pac - (bns->l_pac + 3) / 4 * 4) / 4);
		for (l = bns->l_pac - 1; l >= 0; --l, ++bns->l_pac)
			_set_pac(pac, bns->l_pac, 3 - _get_pac(pac, l));
	}
	ret = bns->l_pac;
	{ // finalize .pac file
		ubyte_t ct;
		err_fwrite(pac, 1, (bns->l_pac >> 2) + ((bns->l_pac & 3) == 0 ? 0 : 1), fp);
		// the following codes make the pac file size always (l_pac/4+1+1)
		if (bns->l_pac % 4 == 0)
		{
			ct = 0;
			err_fwrite(&ct, 1, 1, fp);
		}
		ct = bns->l_pac % 4;
		err_fwrite(&ct, 1, 1, fp);
		// close .pac file
		err_fflush(fp);
		err_fclose(fp);
	}
	bns_dump(bns, prefix);
	bns_destroy(bns);
	kseq_destroy(seq);
	free(pac);
	return ret;
}

int bwa_fa2pac(int argc, char *argv[])
{
	int c, for_only = 0;
	gzFile fp;
	while ((c = getopt(argc, argv, "f")) >= 0)
	{
		switch (c)
		{
		case 'f':
			for_only = 1;
			break;
		}
	}
	if (argc == optind)
	{
		fprintf(stderr, "Usage: bwa fa2pac [-f] <in.fasta> [<out.prefix>]\n");
		return 1;
	}
	fp = xzopen(argv[optind], "r");
	bns_fasta2bntseq(fp, (optind + 1 < argc) ? argv[optind + 1] : argv[optind], for_only);
	err_gzclose(fp);
	return 0;
}

int bns_pos2rid(const bntseq_t *bns, int64_t pos_f)
{
	int left, mid, right;
	if (pos_f >= bns->l_pac)
		return -1;
	left = 0;
	mid = 0;
	right = bns->n_seqs;
	while (left < right)
	{ // binary search
		mid = (left + right) >> 1;
		if (pos_f >= bns->anns[mid].offset)
		{
			if (mid == bns->n_seqs - 1)
				break;
			if (pos_f < bns->anns[mid + 1].offset)
				break; // bracketed
			left = mid + 1;
		}
		else
			right = mid;
	}
	return mid;
}

int bns_intv2rid(const bntseq_t *bns, int64_t rb, int64_t re)
{
	int is_rev, rid_b, rid_e;
	if (rb < bns->l_pac && re > bns->l_pac)
		return -2;
	assert(rb <= re);
	rid_b = bns_pos2rid(bns, bns_depos(bns, rb, &is_rev));
	rid_e = rb < re ? bns_pos2rid(bns, bns_depos(bns, re - 1, &is_rev)) : rid_b;
	return rid_b == rid_e ? rid_b : -1;
}

int bns_cnt_ambi(const bntseq_t *bns, int64_t pos_f, int len, int *ref_id)
{
	int left, mid, right, nn;
	if (ref_id)
		*ref_id = bns_pos2rid(bns, pos_f);
	left = 0;
	right = bns->n_holes;
	nn = 0;
	while (left < right)
	{
		mid = (left + right) >> 1;
		if (pos_f >= bns->ambs[mid].offset + bns->ambs[mid].len)
			left = mid + 1;
		else if (pos_f + len <= bns->ambs[mid].offset)
			right = mid;
		else
		{ // overlap
			if (pos_f >= bns->ambs[mid].offset)
			{
				nn += bns->ambs[mid].offset + bns->ambs[mid].len < pos_f + len ? bns->ambs[mid].offset + bns->ambs[mid].len - pos_f : len;
			}
			else
			{
				nn += bns->ambs[mid].offset + bns->ambs[mid].len < pos_f + len ? bns->ambs[mid].len : len - (bns->ambs[mid].offset - pos_f);
			}
			break;
		}
	}
	return nn;
}

uint8_t *bns_get_seq(int64_t l_pac, const uint8_t *pac, int64_t beg, int64_t end, int64_t *len)
{
	uint8_t *seq = 0;
	if (end < beg)
		end ^= beg, beg ^= end, end ^= beg; // if end is smaller, swap
	if (end > l_pac << 1)
		end = l_pac << 1;
	if (beg < 0)
		beg = 0;
	if (beg >= l_pac || end <= l_pac)
	{
		int64_t k, l = 0;
		*len = end - beg;
		seq = (uint8_t *)malloc(end - beg + 64);
		assert(seq != NULL);
		if (beg >= l_pac)
		{ // reverse strand
			int64_t beg_f = (l_pac << 1) - 1 - end;
			int64_t end_f = (l_pac << 1) - 1 - beg;
			for (k = end_f; k > beg_f; --k)
			{
				seq[l++] = 3 - _get_pac(pac, k);
			}
		}
		else
		{ // forward strand
			for (k = beg; k < end; ++k)
			{
				seq[l++] = _get_pac(pac, k);
			}
		}
	}
	else
		*len = 0; // if bridging the forward-reverse boundary, return nothing
	return seq;
}

uint8_t *bns_fetch_seq(const bntseq_t *bns, const uint8_t *pac, int64_t *beg, int64_t mid, int64_t *end, int *rid)
{
	int64_t far_beg, far_end, len;
	int is_rev;
	uint8_t *seq;

	if (*end < *beg)
		*end ^= *beg, *beg ^= *end, *end ^= *beg; // if end is smaller, swap
	// printf("%ld %ld %ld\n", *beg, mid, *end);
	assert(*beg <= mid && mid < *end);

	*rid = bns_pos2rid(bns, bns_depos(bns, mid, &is_rev));
	far_beg = bns->anns[*rid].offset;
	far_end = far_beg + bns->anns[*rid].len;
	if (is_rev)
	{ // flip to the reverse strand
		int64_t tmp = far_beg;
		far_beg = (bns->l_pac << 1) - far_end;
		far_end = (bns->l_pac << 1) - tmp;
	}
	*beg = *beg > far_beg ? *beg : far_beg;
	*end = *end < far_end ? *end : far_end;

	seq = bns_get_seq(bns->l_pac, pac, *beg, *end, &len);

	if (seq == 0 || *end - *beg != len)
	{
		fprintf(stderr, "[E::%s] begin=%ld, mid=%ld, end=%ld, len=%ld, seq=%p, rid=%d, far_beg=%ld, far_end=%ld\n",
				__func__, (long)*beg, (long)mid, (long)*end, (long)len, seq, *rid, (long)far_beg, (long)far_end);
	}
	assert(seq && *end - *beg == len); // assertion failure should never happen
	return seq;
}
