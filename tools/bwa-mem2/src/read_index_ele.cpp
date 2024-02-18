/*************************************************************************************
                           The MIT License

   BWA-MEM2  (Sequence alignment using Burrows-Wheeler Transform),
   Copyright (C) 2019  Intel Corporation, Heng Li.

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

#include "read_index_ele.h"
#ifdef __cplusplus
extern "C" {
#endif
#include "safe_mem_lib.h"
#include "safe_str_lib.h"
#include <chrono>
using namespace std::chrono;

#ifdef __cplusplus
}
#endif

indexEle::indexEle()
{
    idx = (bwaidx_fm_t*) calloc(1, sizeof(bwaidx_fm_t));
    assert(idx != NULL);
}

indexEle::~indexEle()
{
    if (idx == 0) return;
    if (idx->mem == 0)
    {
        if (idx->bns) bns_destroy(idx->bns);
        if (idx->pac) free(idx->pac);
    } else {
        free(idx->bns->anns); free(idx->bns);
        if (!idx->is_shm) free(idx->mem);
    }
    free(idx);
}

// #define BWA_IDX_BWT 0x1
// #define BWA_IDX_BNS 0x2
// #define BWA_IDX_PAC 0x4
// #define BWA_IDX_ALL 0x7

// bwa_idx_load_ele(ref_file_name, BWA_IDX_ALL);
void indexEle::bwa_idx_load_ele(const char* hint, int which)            
{
    auto bwa_idx_load_ele_start_1 = high_resolution_clock::now();

    // char *prefix;
    vecCT prefix;
    int l_hint = strlen(hint);
    // int l_hint=hint.size()-1;
    // prefix = (char *) malloc(l_hint + 3 + 4 + 1);
    prefix.resize(l_hint + 3 + 4 + 1);
    // assert(prefix != NULL);
    assert(prefix.size()!=0);
    // strcpy_s(prefix, l_hint + 3 + 4 + 1, hint);
    char* hh=(char*)(malloc(strlen(hint)*sizeof(char))); strcpy(hh,hint);
    assign_string_to_vecCT(prefix,hh,-1);

    fprintf(stderr, "* Index prefix: %s\n", convert_ciphertext_vector_to_plaintext_string(prefix));    // Index prefix: ../../input-datasets/fmi/broad
    int duration_2;
    auto bwa_idx_load_ele_end_1 = high_resolution_clock::now();

    // printf("which & BWA_IDX_BNS : %d\n",which & BWA_IDX_BNS);
    idx = (bwaidx_fm_t*) calloc(1, sizeof(bwaidx_fm_t));
    if (which & BWA_IDX_BNS) {
        // printf("Inside\n");
        int i, c;
        
        cout<<"Going into bns_restore!\n"<<endl;
        idx->bns = bns_restore(convert_ciphertext_vector_to_plaintext_string(prefix));
        // cout<<"Out of bns_restore!\n"<<endl;
        auto bwa_idx_load_ele_start_2 = high_resolution_clock::now();

        if (idx->bns == 0) {
            printf("Error!! : [%s] bns is NULL!!\n", __func__);
            exit(EXIT_FAILURE);
        }
        // printf("After checking idx->bns==0\n");
        // cout<<"idx->bns->n_seqs: "<<idx->bns->n_seqs<<endl;
        
        cout<<"decrypt_ciphertext_to_plaintext_vector(idx->bns->n_seqs_enc)[0]: "<<decrypt_ciphertext_to_plaintext_vector(idx->bns->n_seqs_enc)[0]<<endl;
        // for (i = c = 0; operate_and_decrypt(idx->bns->n_seqs_enc,"-",i)>0; ++i)
        // {
        //     cout<<"i: "<<i<<endl;
        //     if (idx->bns->anns[i].is_alt) 
        //     {
        //         ++c;
        //         cout<<"New c: "<<c<<endl;
        //     }
        // }

        for (i = c = 0; i<decrypt_ciphertext_to_plaintext_vector(idx->bns->n_seqs_enc)[0]; ++i)
        {
            // printf("i: %d\n",i);
            // cout<<"idx->bns->anns[i].is_alt: "<<idx->bns->anns[i].is_alt<<endl;
            // cout<<"idx->bns->anns[i].is_alt_enc: "<<decrypt_ciphertext_to_plaintext_vector(idx->bns->anns[i].is_alt_enc)[0]<<endl;
            if (decrypt_ciphertext_to_plaintext_vector(idx->bns->anns[i].is_alt_enc)[0]) 
                {
                    ++c;
                    // printf("new c: %d\n",c);
                }
        }
        
        fprintf(stderr, "* Read %d ALT contigs\n", c);
        
        // if (which & BWA_IDX_PAC)
        // {
        //     idx->pac = (uint8_t*) calloc(idx->bns->l_pac/4+1, 1); idx->pac_enc.resize(idx->bns->l_pac/4+1);
        //     assert(idx->pac != NULL);
        //     err_fread_noeof(idx->pac, 1, idx->bns->l_pac/4+1, idx->bns->fp_pac); // concatenated 2-bit encoded sequence
        //     err_fclose(idx->bns->fp_pac);

        //     for(int h=0;h<idx->bns->l_pac/4+1; h++)
        //     {
        //         idx->pac_enc[h]=encrypt_plaintext_integer_to_ciphertext(idx->pac[h]);
        //     }

        //     idx->bns->fp_pac = 0;
        // }

        if (which & BWA_IDX_PAC)
        {
            idx->pac = (uint8_t*) calloc(decrypt_ciphertext_to_plaintext_vector(idx->bns->l_pac_enc)[0]/4+1, 1);
            assert(idx->pac != NULL);
            err_fread_noeof(idx->pac, 1, decrypt_ciphertext_to_plaintext_vector(idx->bns->l_pac_enc)[0]/4+1, idx->bns->fp_pac); // concatenated 2-bit encoded sequence
            err_fclose(idx->bns->fp_pac);
            idx->bns->fp_pac = 0;
        }

        auto bwa_idx_load_ele_end_2 = high_resolution_clock::now();
        duration_2=duration_cast<microseconds>(bwa_idx_load_ele_end_2-bwa_idx_load_ele_start_2).count();
    }
    // free(prefix);
    prefix.clear();
    

    auto duration_bwa_idx_load_ele = duration_cast<microseconds>(bwa_idx_load_ele_end_1-bwa_idx_load_ele_start_1);
    cout << "Time taken by bwa_idx_load_ele: "<<duration_bwa_idx_load_ele.count()+duration_2<< " microseconds" << endl;
}

#include <sys/file.h>
char* indexEle::bwa_idx_infer_prefix(const char *hint)
{
    char *prefix;
    int l_hint;
    FILE *fp;
    l_hint = strlen(hint);
    prefix = (char *) malloc(l_hint + 3 + 4 + 1);
    assert(prefix != NULL);
    strcpy_s(prefix, l_hint + 3 + 4 + 1, hint);
    strcpy_s(prefix + l_hint, 8, ".64.bwt");
    if ((fp = fopen(prefix, "rb")) != 0)
    {
        fclose(fp);
        prefix[l_hint + 3] = 0;
        return prefix;
    } else {
        strcpy_s(prefix + l_hint, 8, ".bwt");
        if ((fp = fopen(prefix, "rb")) == 0)
        {
            free(prefix);
            return 0;
        } else {
            //flock(fileno(fp), 1);
            //flock(fileno(fp), 1);  // Unlock the file
            fclose(fp);
            prefix[l_hint] = 0;
            return prefix;
        }
    }
}

#if TEST
//int main(int argc, char* argv[])
//{
//  printf("Testing read_index_ele...\n");
//  indexEle *bwaEle = new indexEle();
//  
//  bwaEle->bwa_idx_load_ele("/projects/PCL-GBB/wasim/read_and_ref_data_1/hgaa.fa",
//                          BWA_IDX_ALL);
//
//  delete bwaEle;
//}
#endif
