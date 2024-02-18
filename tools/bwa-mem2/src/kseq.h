/* The MIT License

   Copyright (c) 2008, 2009, 2011 Attractive Chaos <attractor@live.co.uk>

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


#ifndef AC_KSEQ_H
#define AC_KSEQ_H

#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include "memcpy_bwamem.h"
#include "../../../palisade_header.h"

#include <ctime>

//# include "palisade.h"
// using namespace std;
// using namespace lbcrypto;

// #ifndef total_io_time_def
// # define total_io_time_def
// double total_io_time =  0;
// #endif 

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

#define KS_SEP_SPACE 0 // isspace(): \t, \n, \v, \f, \r
#define KS_SEP_TAB   1 // isspace() && !' '
#define KS_SEP_LINE  2 // line separator: "\n" (Unix) or "\r\n" (Windows)
#define KS_SEP_MAX   2

#define __KS_TYPE(type_t)						\
	typedef struct __kstream_t {				\
		unsigned char *buf;						\
		int begin, end, is_eof;					\
		type_t f;								\
		vecCT enc_buf;							\
	} kstream_t;

#define ks_eof(ks) ((ks)->is_eof && (ks)->begin >= (ks)->end)
#define ks_rewind(ks) ((ks)->is_eof = (ks)->begin = (ks)->end = 0)   

#define __KS_BASIC(type_t, __bufsize)								\
	static inline kstream_t *ks_init(type_t f)						\
	{																\
		kstream_t *ks = (kstream_t*)calloc(1, sizeof(kstream_t));	\
        assert(ks != NULL);                                         \
		ks->f = f;													\
		ks->buf = (unsigned char*)malloc(__bufsize);				\
		ks->enc_buf.resize(__bufsize);								\
		/*ks->enc_buf.shrink_to_fit();		*/						\
		/* printf("ks->enc_buf.size(): %d\n",ks->enc_buf.size());	*/						\
        assert(ks->buf != NULL);                                    \
		assert(ks->enc_buf.size()!=0);								\
		return ks;													\
	}																\
	static inline void ks_destroy(kstream_t *ks)					\
	{																\
		if (ks) {													\
			free(ks->buf);											\
			ks->enc_buf.clear();									\
			ks->enc_buf.shrink_to_fit();							\
			free(ks);												\
		}															\
	}

#define __KS_GETC(__read, __bufsize)						\
	/* static inline int ks_getc(kstream_t *ks)	*/			\
	static inline Ciphertext<DCRTPoly> ks_getc(kstream_t *ks)	\
	{														\
		clock_t begin, end;			\
		/* printf("Entered ks_getc, ks->begin: %d, ks->end: %d, ks->is_eof: %d\n", ks->begin, ks->end, ks->is_eof);		*/					\
		if (ks->is_eof && ks->begin >= ks->end) return encrypt_plaintext_integer_to_ciphertext(-1);	\
		if (ks->begin >= ks->end) {							\
			/* begin=clock();				*/					\
			ks->begin = 0;									\
			/* printf("1 Reading into ks->buf via __read()\n");	*/	\
			ks->end = __read(ks->f, ks->buf, __bufsize);	\
			if (ks->end == 0) { ks->is_eof = 1; 			\
			return encrypt_plaintext_integer_to_ciphertext(-1);}	\
			/* ks->enc_buf.clear();				*/			\
			/* printf("Assigning from ks->buf to ks->enc_buf; ks->end: %d\n",ks->end);	*/\
			/* printf("&ks->enc_buf[0]: %x\n",&ks->enc_buf[0]);	*/\
			for(int i=0; i<ks->end; i++)					\
			{												\
				int c = (int)(ks->buf[i]);					\
				/*printf("#####\n");*/						\
				/* printf("%c",c);		*/			\
				/* clock_t c_begin=clock();			*/					\
				Ciphertext<DCRTPoly> ciphertext = encrypt_plaintext_integer_to_ciphertext(c);		\
				/* clock_t c_end=clock();									\
				printf("time for one byte encryption: %lf\n",double(c_end - c_begin) / CLOCKS_PER_SEC);			*/								\
				ks->enc_buf[i] = ciphertext;				\
			}															\
			/* printf("Assigned!; sizeof(kstream ks->buf): %lu; sizeof(kstream ks->buf)/sizeof(ks->buf[0]): %lu\n",sizeof(ks->buf),sizeof(ks->buf)/sizeof(ks->buf[0]));			*/			\
			/* end = clock();						\
			double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;							\
			printf("time for one read: %lf\n",elapsed_secs);			*/								\
			/* total_io_time+=elapsed_secs;	*/															\
			/* printf("total_io_time: %lf\n",total_io_time);	*/											\
			/* printf("ks->end: %d\n",ks->end);		*/		\
		}													\
		/* return (int)ks->buf[ks->begin++];*/				\
		/* cout<<"ks_getc is returning: \n"; /*<<ks->enc_buf[ks->begin++]<<endl;*/	\
		/* printf("Leaving ks_getc, ks->begin: %d, ks->end: %d, ks->is_eof: %d\n", ks->begin, ks->end, ks->is_eof);				*/			\
		return ks->enc_buf[ks->begin++];					\
	}

#ifndef KSTRING_T
#define KSTRING_T kstring_t		
typedef struct __kstring_t {		
	size_t l, m;		
	char* s; vecCT enc_s;			
} kstring_t;		
#endif

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif
/*str->l = (i-ks->begin) is returned in case of comment; 1+(i-ks->begin) returned in main case; so str->l is total size of chars copied*/
#define __KS_GETUNTIL(__read, __bufsize)								\
	static int ks_getuntil2(kstream_t *ks, int delimiter, kstring_t *str, Ciphertext<DCRTPoly> *dret, int append) \
	{																	\
		/* printf("At start, ks->begin: %d, ks->end: %d\n",ks->begin,ks->end);		*/			\
		clock_t begin, end;			\
		double elapsed_secs;					\
		int gotany = 0;													\
		if (dret) *dret = 0;											\
		str->l = append? str->l : 0;									\
		for (;;) {														\
			int i;														\
			if (ks->begin >= ks->end) {									\
				if (!ks->is_eof) {										\
					/* begin = clock();		*/							\
					ks->begin = 0;										\
					/* printf("2 Reading into ks->buf via __read()\n");	*/	\
					ks->end = __read(ks->f, ks->buf, __bufsize);		\
					/* printf("Assigning from ks->buf to ks->enc_buf ; ks->end: %d\n",ks->end);	*/\
					int c=0;	CT c_ct;										\
					/* printf("&ks->enc_buf[0]: %x\n",&ks->enc_buf[0]);	*/\
					for(int h=0; h<ks->end; h++)						\
					{													\
						printf("%d ",h);											\
						c = (int)(ks->buf[h]);						\
						/* clock_t c_begin=clock();				*/					\
						c_ct = encrypt_plaintext_integer_to_ciphertext(c);		\
						/* printf("original ks->enc_buf[h]: %c\n",decrypt_ciphertext_to_plaintext_vector(ks->enc_buf[h])[0]);	*/	\
						/* clock_t c_end=clock();									\
						printf("time for one byte encryption: %lf\n",double(c_end - c_begin) / CLOCKS_PER_SEC);				*/							\
						ks->enc_buf[h]= c_ct;				\
					}													\
					/* printf("Assigned!; sizeof(kstream ks->buf): %lu; sizeof(kstream ks->buf)/sizeof(kstream ks->buf[0]): %lu\n",sizeof(ks->buf),sizeof(ks->buf)/sizeof(ks->buf[0]));		*/	\
					/* end = clock();						\
					elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;							\
					printf("time for one read: %lf\n",elapsed_secs);			*/								\
					/* total_io_time+=elapsed_secs;		*/														\
					/* printf("total_io_time: %lf\n",total_io_time);		*/										\
					if (ks->end == 0) { ks->is_eof = 1; break; }		\
					/* printf("Read new chars into ks->buf; now ks->begin: %d, ks->end: %d\n",ks->begin,ks->end);*/			\
				} else break;											\
			}															\
			if (delimiter == KS_SEP_LINE) { \
				for (i = ks->begin; i < ks->end; ++i) \
					if (decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(ks->enc_buf[i], encode_integer_to_plaintext('\n')))[0]==0) break; \
			} else if (delimiter > KS_SEP_MAX) {						\
				for (i = ks->begin; i < ks->end; ++i)					\
					if (decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(ks->enc_buf[i], encode_integer_to_plaintext(delimiter)))[0]==0) break;					\
			} else if (delimiter == KS_SEP_SPACE) {						\
				for (i = ks->begin; i < ks->end; ++i)					\
					if (decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(ks->enc_buf[i], encode_integer_to_plaintext('\n')))[0]==0 || \
					decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(ks->enc_buf[i], encode_integer_to_plaintext(' ')))[0]==0 ||	\
					decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(ks->enc_buf[i], encode_integer_to_plaintext('\f')))[0]==0 || \
					decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(ks->enc_buf[i], encode_integer_to_plaintext('\r')))[0]==0 || \
					decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(ks->enc_buf[i], encode_integer_to_plaintext('\t')))[0]==0 || \
					decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(ks->enc_buf[i], encode_integer_to_plaintext('\v')))[0]==0)	\
						break;						\
			} else if (delimiter == KS_SEP_TAB) {						\
				for (i = ks->begin; i < ks->end; ++i)					\
					if( (decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(ks->enc_buf[i], encode_integer_to_plaintext('\n')))[0]==0 || \
					decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(ks->enc_buf[i], encode_integer_to_plaintext(' ')))[0]==0 ||	\
					decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(ks->enc_buf[i], encode_integer_to_plaintext('\f')))[0]==0 || \
					decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(ks->enc_buf[i], encode_integer_to_plaintext('\r')))[0]==0 || \
					decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(ks->enc_buf[i], encode_integer_to_plaintext('\t')))[0]==0 || \
					decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(ks->enc_buf[i], encode_integer_to_plaintext('\v')))[0]==0) && decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(ks->enc_buf[i], encode_integer_to_plaintext(' ')))[0]!=0) break; \
			} else i = 0; /* never come to here! */						\
			/*printf("i: %d\n",i);		/*i is index of \n in ks->buf in case of KS_SEP_LINE*/	\
			/* Below runs if ks->buf[i]==delimiter (like \n) is 'not' encountered in the string read into ks->buf from file */	\
			/* printf("Earlier, str->m: %d, str->l: %d, i: %d, ks->begin: %d, str->enc_s.size(): %d\n",str->m,str->l,i,ks->begin,str->enc_s.size());	*/						\
			if (str->m - str->l < (size_t)(i - ks->begin + 1)) {		\
				str->m = str->l + (i - ks->begin) + 1;					\
				kroundup32(str->m);										\
				/* str->s = (char*)calloc(str->m,sizeof(char));	*/		\
				str->enc_s.resize(str->m);									\
			}															\
			gotany = 1;													\
			/* memcpy_bwamem(str->s + str->l, str->m - str->l, ks->buf + ks->begin, i - ks->begin, __FILE__, __LINE__); */ \
			memcpy_bwamem_enc(&str->enc_s, str->l, str->m - str->l, &ks->enc_buf, ks->begin, i - ks->begin, __FILE__, __LINE__); \
			/* printf("Out of memcpy_bwamem_enc!\n");		*/	\
			/* vector<CT> vv(str->enc_s.begin()+str->l,str->enc_s.end());		\
			char* enc_s_string=convert_ciphertext_vector_to_plaintext_string(vv);			\
			printf("enc_s_string: %s\n",enc_s_string);										\
			strcpy(str->s+str->l,enc_s_string);							*/	\
			/* printf("Out of memcpy_bwamem!\n");	*/		\
			/* (i-ks->begin) bytes are copied from (ks->buf+ks->begin) to (str->s+str->l); str->m-str->l = ks->end-ks->begin+1 */	\
			/* printf("After memcpy_bwamem, str->s+str->l: %s, str->l: %d, str->m: %d\n",str->s+str->l,str->l,str->m);	*/\
			str->l = str->l + (i - ks->begin);							\
			ks->begin = i + 1;											\
			/* printf("neww str->l: %d, ks->begin=i+1= %d\n",str->l,ks->begin);		*/	\
			/* If i==ks->end, then while reading from ks->buf, ks->begin reached ks->end, so, i=ks->begin=ks->end now. But we are yet to read the delimiter. */	\
			/*If before memcpy_bwamem, ks->begin: 16319, ks->end: 16384, */	\
			if (i < ks->end) {											\
				/*ks->buf[i] is DELIMITER*/								\
				/* printf("i: %d, ks->end: %d\n",i,ks->end);	*/			\
				/* printf("decrypt_ciphertext_to_plaintext_vector(ks->enc_buf[i])[0]: %c\n",(char)(decrypt_ciphertext_to_plaintext_vector(ks->enc_buf[i])[0]));	 */	\
				if (dret) *dret = ks->enc_buf[i];						\
				break;													\
			}															\
		}																\
		/* printf("Out of loop!\n");	*/							\
		if (!gotany && ks_eof(ks)) return -1;							\
		/* printf("str->s+str->l: %s, str->l: %d, str->m: %d\n",str->s+str->l,str->l,str->m);			\
		printf("str->enc_s.size(): %d\n",str->enc_s.size());				\
		vector<CT> vv(str->enc_s.begin()+str->l,str->enc_s.end());		\
		printf("convert_ciphertext_vector_to_plaintext_string(vv): %s\n",convert_ciphertext_vector_to_plaintext_string(vv));		*/	\
		if ((str->enc_s).size() == 0) {												\
			str->m = 1;													\
			/* str->s = (char*)calloc(1, 1);	*/						\
			(str->enc_s).resize(1);											\
			/* (str->enc_s).shrink_to_fit();	*/							\
            /* assert(str->s != NULL);    */                                 \
			assert((str->enc_s).size()!=0);									\
		} else if (delimiter == KS_SEP_LINE && str->l > 1 && decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(str->enc_s[str->l-1] , encode_integer_to_plaintext('\r')))[0]==0) --str->l; \
		/* printf("next!\n");		*/	\
		/* str->s[str->l] = '\0';		*/								\
		str->enc_s[str->l]=encrypt_plaintext_integer_to_ciphertext('\0');			\
		/* printf("At end, ks->begin: %d, str->l: %d\n\n",ks->begin,str->l); 	*/ \
		return str->l;													\
	}																	\
	static inline int ks_getuntil(kstream_t *ks, int delimiter, kstring_t *str, Ciphertext<DCRTPoly> *dret) \
	{ return ks_getuntil2(ks, delimiter, str, dret, 0); }

#define KSTREAM_INIT(type_t, __read, __bufsize) \
	__KS_TYPE(type_t)							\
	__KS_BASIC(type_t, __bufsize)				\
	__KS_GETC(__read, __bufsize)				\
	__KS_GETUNTIL(__read, __bufsize)

#define kseq_rewind(ks)  ((ks)->f->is_eof = (ks)->f->begin = (ks)->f->end = 0; (ks)->last_char = 0; (ks)->last_char = encrypt_plaintext_integer_to_ciphertext(0);) 

#define __KSEQ_BASIC(SCOPE, type_t)										\
	SCOPE kseq_t *kseq_init(type_t fd)									\
	{																	\
		kseq_t *s = (kseq_t*)calloc(1, sizeof(kseq_t));					\
		/* s->last_char = 0;		*/	\
		s->last_char = encrypt_plaintext_integer_to_ciphertext(0);		\
        assert(s != NULL);                                              \
		s->f = ks_init(fd);												\
		return s;														\
	}																	\
	SCOPE void kseq_destroy(kseq_t *ks)									\
	{																	\
		if (!ks) return;												\
		free(ks->name.s); free(ks->comment.s); free(ks->seq.s);	free(ks->qual.s);   \
		(ks->name.enc_s).clear(); (ks->comment.enc_s).clear(); (ks->seq.enc_s).clear(); (ks->qual.enc_s).clear(); 	\
		(ks->name.enc_s).shrink_to_fit(); (ks->comment.enc_s).shrink_to_fit(); (ks->seq.enc_s).shrink_to_fit(); (ks->qual.enc_s).shrink_to_fit(); \
		ks_destroy(ks->f);												\
		free(ks);														\
	}

/* Return value:
   >=0  length of the sequence (normal)
   -1   end-of-file
   -2   truncated quality string
 */
/* ks_getc() returns ks->enc_buf[ks->begin++] */
#define __KSEQ_READ(SCOPE) \
	SCOPE int64_t kseq_read(kseq_t *seq) \
	{ \
		/* printf("########################################################\n");			*/							\
		Ciphertext<DCRTPoly> c; \
		kstream_t *ks = seq->f; \
		if (compare_enc(seq->last_char,0)) { /* then jump to the next header line */ \
			while(!compare_enc((c = ks_getc(ks)),-1) && !compare_enc(c,'>') && !compare_enc(c,'@'));		\
			if (compare_enc(c,-1)) return -1; /* end of file */ \
			/* clock_t cmp_begin=clock();			\
			int cmp_ret=(decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(c, encode_integer_to_plaintext(-1)))[0]!=0);				\
			clock_t cmp_end=clock();				\
			double cmp_elapsed_secs = double(cmp_end - cmp_begin) / CLOCKS_PER_SEC;							\
			printf("time for one cmp comparison operation: %lf\n",cmp_elapsed_secs);										*/	\
			/* Ciphertext<DCRTPoly> c2=encode_integer_to_plaintext(2);		*/	\
			/* clock_t mul_begin=clock();	*/		\
			/* cc->EvalMult(c, c2);		*/		\
			/* clock_t mul_end=clock();   */				\
			/* double mul_elapsed_secs = double(mul_end - mul_begin) / CLOCKS_PER_SEC;							\
			printf("time for one mul operation: %lf\n\n",mul_elapsed_secs);			*/								\
			seq->last_char = c; \
			/* printf("decrypt_ciphertext_to_plaintext_vector(seq->last_char)[0]: %c\n",decrypt_ciphertext_to_plaintext_vector(seq->last_char)[0]);	*/	\
			/* seq->last_char=decrypt_ciphertext_to_plaintext_vector(c)[0];	*/	\
		} /* else: the first header char has been read in the previous call */ \
		seq->comment.l = seq->seq.l = seq->qual.l = 0; /* reset all members */ \
		if (ks_getuntil(ks, 0, &seq->name, &c) < 0) return -1; /*0 here means read till space; normal exit: EOF */ \
		/* printf("convert_ciphertext_vector_to_plaintext_string(seq->name.enc_s): %s\n",convert_ciphertext_vector_to_plaintext_string(seq->name.enc_s)); */\
		if (!compare_enc(c,'\n')) 						\
		{							\
			ks_getuntil(ks, KS_SEP_LINE, &seq->comment, 0); /* read FASTA/Q comment */ \
			/*printf("seq->comment.s: %s\n",seq->comment.s);*/	\
			vecInt comment=convert_ciphertext_vector_to_plaintext_vector(seq->comment.enc_s);	\
			for(int i=0;i<comment.size();i++)		\
				cout<<(char)(comment[i]);					\
			cout<<endl;								\
		}					\
		if ((seq->seq.enc_s).size()==0) { /* we can do this in the loop below, but that is slower */ \
			seq->seq.m = 256; \
			/* seq->seq.s = (char*)malloc(seq->seq.m);*/ /* seq.m is the total space for seq.s, and seq.l is current position of character read */ \
			vecCT v(256,0); \
			seq->seq.enc_s=v;									\
            assert((seq->seq.enc_s).size()!=0);             \
		} \
		/*The below loop will run >1 times; in 1st iteration, it will read the main sequence, and in next, it will read + and then end*/	\
		while (decrypt_ciphertext_to_plaintext_vector(cc->EvalSub((c = ks_getc(ks)), encode_integer_to_plaintext(-1)))[0]!=0 && decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(c, encode_integer_to_plaintext('>')))[0]!=0 && decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(c, encode_integer_to_plaintext('+')))[0]!=0 && decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(c, encode_integer_to_plaintext('@')))[0]!=0) { \
			if (decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(c, encode_integer_to_plaintext('\n')))[0]==0) continue; /* skip empty lines */ \
			seq->seq.enc_s[seq->seq.l++] = c; /* this is safe: we always have enough space for 1 char */ \
			ks_getuntil2(ks, KS_SEP_LINE, &seq->seq, 0, 1); /* read the rest of the line */ \
			/* printf("seq->seq.s: %s\n",seq->seq.s); */\
			/* printf("convert_ciphertext_vector_to_plaintext_string(seq->seq.enc_s): %s\n",convert_ciphertext_vector_to_plaintext_string(seq->seq.enc_s)); */\
		} \
		if (decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(c, encode_integer_to_plaintext('>')))[0]==0 || decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(c, encode_integer_to_plaintext('@')))[0]==0) {	\
			/* seq->last_char = decrypt_ciphertext_to_plaintext_vector(c)[0]; /* the first header char has been read */	\
			seq->last_char=c;	}\
		if (seq->seq.l + 1 >= seq->seq.m) { /* seq->seq.s[seq->seq.l] below may be out of boundary */ \
			seq->seq.m = seq->seq.l + 16;								\
			kroundup32(seq->seq.m); /* rounded to the next closest 2^k */ \
			(seq->seq.enc_s).resize(seq->seq.m);			\
			/* seq->seq.s = (char*)calloc(seq->seq.m,sizeof(char));	*/	\
            assert((seq->seq.enc_s).size()!=0);					\
		} \
		seq->seq.enc_s[seq->seq.l] = encrypt_plaintext_integer_to_ciphertext(0);	/* null terminated string */ \
		if (decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(c, encode_integer_to_plaintext('+')))[0]!=0) return seq->seq.l; /* FASTA */ \
		if (seq->qual.m < seq->seq.m) {	/* allocate memory for qual in case insufficient */ \
			seq->qual.m = seq->seq.m; \
			/* seq->qual.s = (char*)calloc(seq->qual.m,sizeof(char)); */ \
			(seq->qual.enc_s).resize(seq->qual.m);				\
		} \
		while (decrypt_ciphertext_to_plaintext_vector(cc->EvalSub((c = ks_getc(ks)), encode_integer_to_plaintext(-1)))[0]!=0 && decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(c, encode_integer_to_plaintext('\n')))[0]!=0); /* skip the rest of '+' line */ \
		if (decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(c, encode_integer_to_plaintext(-1)))[0]==0) return -2; /* error: no quality string */ \
		while (ks_getuntil2(ks, KS_SEP_LINE, &seq->qual, 0 , 1) >= 0 && seq->qual.l < seq->seq.l); \
		/* printf("convert_ciphertext_vector_to_plaintext_string(seq->qual.enc_s): %s\n",convert_ciphertext_vector_to_plaintext_string(seq->qual.enc_s)); */\
		/* seq->last_char=0;		*/	\
		seq->last_char = encrypt_plaintext_integer_to_ciphertext(0);	/* we have not come to the next header line */ \
		/* printf("seq->name.l: %d, seq->comment.l: %d, seq->seq.l: %d, seq->qual.l: %d\n\n",seq->name.l,seq->comment.l,seq->seq.l,seq->qual.l); */		\
		/* cout<<"seq->last_char: "<<decrypt_ciphertext_to_plaintext_vector(seq->last_char)[0]<<"; seq->seq.l: "<<seq->seq.l<<"; seq->qual.l: "<<seq->qual.l<<endl;	*/	\
		if (seq->seq.l != seq->qual.l) return -2; /* error: qual string is of a different length */ \
		return seq->seq.l; \
	}

#define __KSEQ_TYPE(type_t)						\
	typedef struct {							\
		kstring_t name, comment, seq, qual;		\
		Ciphertext<DCRTPoly> last_char; 					\
		/* Ciphertext<DCRTPoly>  enc_last_char;	*/	\
		kstream_t *f;							\
	} kseq_t;

#define KSEQ_INIT2(SCOPE, type_t, __read)		\
	KSTREAM_INIT(type_t, __read, 16384)			\
	__KSEQ_TYPE(type_t)							\
	__KSEQ_BASIC(SCOPE, type_t)					\
	__KSEQ_READ(SCOPE)

#define KSEQ_INIT(type_t, __read) KSEQ_INIT2(static, type_t, __read)

#define KSEQ_DECLARE(type_t)											\
	__KS_TYPE(type_t)													\
		__KSEQ_TYPE(type_t)												\
		extern kseq_t *kseq_init(type_t fd);							\
	void kseq_destroy(kseq_t *ks);										\
	int64_t kseq_read(kseq_t *seq);											\
	static int ks_getuntil2(kstream_t *ks, int delimiter, kstring_t *str, Ciphertext<DCRTPoly> *dret, int append);
	
#endif
