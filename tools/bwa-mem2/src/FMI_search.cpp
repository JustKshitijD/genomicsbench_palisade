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

Authors: Sanchit Misra <sanchit.misra@intel.com>; Vasimuddin Md <vasimuddin.md@intel.com>;
*****************************************************************************************/

#include <stdio.h>
#include "sais.h"
#include "FMI_search.h"
#include "memcpy_bwamem.h"
#include "profiling.h"
#include "../../../palisade_header.h"

#include <chrono>
using namespace std::chrono;

#ifdef __cplusplus
extern "C"
{
#endif
#include "safe_str_lib.h"
#ifdef __cplusplus
}
#endif

FMI_search::FMI_search(const char *fname)
// FMI_search::FMI_search(vecCT fname)
{
    fprintf(stderr, "* Entering FMI_search\n");
    strcpy(file_name, fname);
    strcpy_s(file_name, PATH_MAX, fname);
    // file_name.resize(PATH_MAX);
    // strdup_enc(fname,file_name);

    reference_seq_len = 0;
    printf("Initializinf FMI_search!!!!\n");
    reference_seq_len_enc = encrypt_plaintext_integer_to_ciphertext(0);
    // reference_seq_len_read_plaintext=0;

    sentinel_index = 0;
    sentinel_index_enc = encrypt_plaintext_integer_to_ciphertext(0);
    // sentinel_index_read_plaintext = 0;

    index_alloc = 0;
    sa_ls_word = NULL;
    sa_ms_byte = NULL;
    cp_occ = NULL;
    one_hot_mask_array = NULL;
    // map<pair<int,int>,pair<pair<int,int>,pair<int,int> > > cp_occ_cp_count_map;        //For cp_occ[i].cp_count[j]-  ((i,j),((ciphertext_1_no,index_in_ciphertext_1),(ciphertext_2_no,index_in_ciphertext_2)))
    // int cp_occ_cp_count_ciphertext_count=0;
}

FMI_search::~FMI_search()
{
    if (sa_ms_byte)
    {
        _mm_free(sa_ms_byte);
        sa_ms_byte_enc.clear();
    }
    if (sa_ls_word)
    {
        _mm_free(sa_ls_word);
        sa_ls_word_enc.clear();
    }
    if (cp_occ)
    {
        _mm_free(cp_occ);
    }
    if (one_hot_mask_array)
        _mm_free(one_hot_mask_array);
}

int64_t FMI_search::pac_seq_len(const char *fn_pac)
{
    FILE *fp;
    int64_t pac_len;
    uint8_t c;
    fp = xopen(fn_pac, "rb");
    err_fseek(fp, -1, SEEK_END);
    pac_len = err_ftell(fp);
    err_fread_noeof(&c, 1, 1, fp);
    err_fclose(fp);
    return (pac_len - 1) * 4 + (int)c;
}

void FMI_search::pac2nt(const char *fn_pac, std::string &reference_seq)
{
    uint8_t *buf2;
    int64_t i, pac_size, seq_len;
    FILE *fp;

    // initialization
    seq_len = pac_seq_len(fn_pac);
    assert(seq_len > 0);
    assert(seq_len <= 0x7fffffffffL);
    fp = xopen(fn_pac, "rb");

    // prepare sequence
    pac_size = (seq_len >> 2) + ((seq_len & 3) == 0 ? 0 : 1);
    buf2 = (uint8_t *)calloc(pac_size, 1);
    assert(buf2 != NULL);
    err_fread_noeof(buf2, 1, pac_size, fp);
    err_fclose(fp);
    for (i = 0; i < seq_len; ++i)
    {
        int nt = buf2[i >> 2] >> ((3 - (i & 3)) << 1) & 3;
        switch (nt)
        {
        case 0:
            reference_seq += "A";
            break;
        case 1:
            reference_seq += "C";
            break;
        case 2:
            reference_seq += "G";
            break;
        case 3:
            reference_seq += "T";
            break;
        default:
            fprintf(stderr, "ERROR! Value of nt is not in 0,1,2,3!");
            exit(EXIT_FAILURE);
        }
    }
    for (i = seq_len - 1; i >= 0; i--)
    {
        char c = reference_seq[i];
        switch (c)
        {
        case 'A':
            reference_seq += "T";
            break;
        case 'C':
            reference_seq += "G";
            break;
        case 'G':
            reference_seq += "C";
            break;
        case 'T':
            reference_seq += "A";
            break;
        }
    }
    free(buf2);
}

int FMI_search::build_fm_index(const char *ref_file_name, char *binary_seq, int64_t ref_seq_len, int64_t *sa_bwt, int64_t *count)
{
    // printf("ref_seq_len = %ld\n", ref_seq_len);
    fflush(stdout);

    char outname[PATH_MAX];

    strcpy_s(outname, PATH_MAX, ref_file_name);
    strcat_s(outname, PATH_MAX, CP_FILENAME_SUFFIX);
    // sprintf(outname, "%s.bwt.2bit.%d", ref_file_name, CP_BLOCK_SIZE);

    std::fstream outstream(outname, std::ios::out | std::ios::binary);
    outstream.seekg(0);

    // printf("count = %ld, %ld, %ld, %ld, %ld\n", count[0], count[1], count[2], count[3], count[4]);
    fflush(stdout);

    uint8_t *bwt;

    ref_seq_len++;
    outstream.write((char *)(&ref_seq_len), 1 * sizeof(int64_t));
    outstream.write((char *)count, 5 * sizeof(int64_t));

    int64_t i;
    int64_t ref_seq_len_aligned = ((ref_seq_len + CP_BLOCK_SIZE - 1) / CP_BLOCK_SIZE) * CP_BLOCK_SIZE;
    int64_t size = ref_seq_len_aligned * sizeof(uint8_t);
    bwt = (uint8_t *)_mm_malloc(size, 64);
    assert_not_null(bwt, size, index_alloc);

    int64_t sentinel_index = -1;
    for (i = 0; i < ref_seq_len; i++)
    {
        if (sa_bwt[i] == 0)
        {
            bwt[i] = 4;
            // printf("BWT[%ld] = 4\n", i);
            sentinel_index = i;
        }
        else
        {
            char c = binary_seq[sa_bwt[i] - 1];
            switch (c)
            {
            case 0:
                bwt[i] = 0;
                break;
            case 1:
                bwt[i] = 1;
                break;
            case 2:
                bwt[i] = 2;
                break;
            case 3:
                bwt[i] = 3;
                break;
            default:
                fprintf(stderr, "ERROR! i = %ld, c = %c\n", i, c);
                exit(EXIT_FAILURE);
            }
        }
    }
    for (i = ref_seq_len; i < ref_seq_len_aligned; i++)
        bwt[i] = DUMMY_CHAR;

    // printf("CP_SHIFT = %d, CP_MASK = %d\n", CP_SHIFT, CP_MASK);
    // printf("sizeof CP_OCC = %ld\n", sizeof(CP_OCC));
    fflush(stdout);
    // create checkpointed occ
    int64_t cp_occ_size = (ref_seq_len >> CP_SHIFT) + 1;
    CP_OCC *cp_occ = NULL;

    size = cp_occ_size * sizeof(CP_OCC);
    cp_occ = (CP_OCC *)_mm_malloc(size, 64);
    assert_not_null(cp_occ, size, index_alloc);
    memset(cp_occ, 0, cp_occ_size * sizeof(CP_OCC));
    int64_t cp_count[16];

    memset(cp_count, 0, 16 * sizeof(int64_t));
    for (i = 0; i < ref_seq_len; i++)
    {
        if ((i & CP_MASK) == 0)
        {
            CP_OCC cpo;
            // cpo.cp_count_enc[0] = encrypt_plaintext_integer_to_ciphertext(cp_count[0]);
            // cpo.cp_count_enc[1] = encrypt_plaintext_integer_to_ciphertext(cp_count[1]);
            // cpo.cp_count_enc[2] = encrypt_plaintext_integer_to_ciphertext(cp_count[2]);
            // cpo.cp_count_enc[3] = encrypt_plaintext_integer_to_ciphertext(cp_count[3]);

            cpo.cp_count[0] = cp_count[0];
            cpo.cp_count[1] = cp_count[1];
            cpo.cp_count[2] = cp_count[2];
            cpo.cp_count[3] = cp_count[3];

            int32_t j;
            // cpo.one_hot_bwt_str_enc[0] = encrypt_plaintext_integer_to_ciphertext(0);
            // cpo.one_hot_bwt_str_enc[1] = encrypt_plaintext_integer_to_ciphertext(0);
            // cpo.one_hot_bwt_str_enc[2] = encrypt_plaintext_integer_to_ciphertext(0);
            // cpo.one_hot_bwt_str_enc[3] = encrypt_plaintext_integer_to_ciphertext(0);

            cpo.one_hot_bwt_str[0] = 0;
            cpo.one_hot_bwt_str[1] = 0;
            cpo.one_hot_bwt_str[2] = 0;
            cpo.one_hot_bwt_str[3] = 0;

            for (j = 0; j < CP_BLOCK_SIZE; j++)
            {
                // cpo.one_hot_bwt_str_enc[0] = shift_left(cpo.one_hot_bwt_str_enc[0],1);
                // cpo.one_hot_bwt_str_enc[1] = shift_left(cpo.one_hot_bwt_str_enc[1],1);
                // cpo.one_hot_bwt_str_enc[2] = shift_left(cpo.one_hot_bwt_str_enc[2],1);
                // cpo.one_hot_bwt_str_enc[3] = shift_left(cpo.one_hot_bwt_str_enc[3],1);

                cpo.one_hot_bwt_str[0] = cpo.one_hot_bwt_str[0] << 1;
                cpo.one_hot_bwt_str[1] = cpo.one_hot_bwt_str[1] << 1;
                cpo.one_hot_bwt_str[2] = cpo.one_hot_bwt_str[2] << 1;
                cpo.one_hot_bwt_str[3] = cpo.one_hot_bwt_str[3] << 1;

                uint8_t c = bwt[i + j];
                // printf("c = %d\n", c);
                if (c < 4)
                {
                    // cpo.one_hot_bwt_str_enc[c] = cc->EvalAdd(cpo.one_hot_bwt_str_enc[c], 1);
                    cpo.one_hot_bwt_str[c] = cpo.one_hot_bwt_str[c] + 1;
                }
            }

            cp_occ[i >> CP_SHIFT] = cpo;
            // set_cp_occ_cp_count_at_index(i >> CP_SHIFT,0,cpo.cp_count[0]);
            // set_cp_occ_cp_count_at_index(i >> CP_SHIFT,1,cpo.cp_count[1]);
            // set_cp_occ_cp_count_at_index(i >> CP_SHIFT,2,cpo.cp_count[2]);
            // set_cp_occ_cp_count_at_index(i >> CP_SHIFT,3,cpo.cp_count[3]);
        }
        cp_count[bwt[i]]++;
    }
    outstream.write((char *)cp_occ, cp_occ_size * sizeof(CP_OCC));
    _mm_free(cp_occ);
    _mm_free(bwt);

#if SA_COMPRESSION

    size = ((ref_seq_len >> SA_COMPX) + 1) * sizeof(uint32_t);
    uint32_t *sa_ls_word = (uint32_t *)_mm_malloc(size, 64);
    assert_not_null(sa_ls_word, size, index_alloc);
    size = ((ref_seq_len >> SA_COMPX) + 1) * sizeof(int8_t);
    int8_t *sa_ms_byte = (int8_t *)_mm_malloc(size, 64);
    assert_not_null(sa_ms_byte, size, index_alloc);
    int64_t pos = 0;
    for (i = 0; i < ref_seq_len; i++)
    {
        if ((i & SA_COMPX_MASK) == 0)
        {
            sa_ls_word[pos] = sa_bwt[i] & 0xffffffff;
            sa_ms_byte[pos] = (sa_bwt[i] >> 32) & 0xff;
            pos++;
        }
    }
    fprintf(stderr, "pos: %d, ref_seq_len__: %ld\n", pos, ref_seq_len >> SA_COMPX);
    outstream.write((char *)sa_ms_byte, ((ref_seq_len >> SA_COMPX) + 1) * sizeof(int8_t));
    outstream.write((char *)sa_ls_word, ((ref_seq_len >> SA_COMPX) + 1) * sizeof(uint32_t));

#else

    size = ref_seq_len * sizeof(uint32_t);
    uint32_t *sa_ls_word = (uint32_t *)_mm_malloc(size, 64);
    assert_not_null(sa_ls_word, size, index_alloc);
    size = ref_seq_len * sizeof(int8_t);
    int8_t *sa_ms_byte = (int8_t *)_mm_malloc(size, 64);
    assert_not_null(sa_ms_byte, size, index_alloc);
    for (i = 0; i < ref_seq_len; i++)
    {
        sa_ls_word[i] = sa_bwt[i] & 0xffffffff;
        sa_ms_byte[i] = (sa_bwt[i] >> 32) & 0xff;
    }
    outstream.write((char *)sa_ms_byte, ref_seq_len * sizeof(int8_t));
    outstream.write((char *)sa_ls_word, ref_seq_len * sizeof(uint32_t));

#endif

    outstream.write((char *)(&sentinel_index), 1 * sizeof(int64_t));
    outstream.close();
    // printf("max_occ_ind = %ld\n", i >> CP_SHIFT);
    fflush(stdout);

    _mm_free(sa_ms_byte);
    _mm_free(sa_ls_word);
    return 0;
}

int FMI_search::build_index()
{

    char *prefix = file_name;
    uint64_t startTick;
    startTick = __rdtsc();
    index_alloc = 0;

    std::string reference_seq;
    char pac_file_name[PATH_MAX];
    strcpy_s(pac_file_name, PATH_MAX, prefix);
    strcat_s(pac_file_name, PATH_MAX, ".pac");
    // sprintf(pac_file_name, "%s.pac", prefix);
    pac2nt(pac_file_name, reference_seq);
    int64_t pac_len = reference_seq.length();
    int status;
    int64_t size = pac_len * sizeof(char);
    char *binary_ref_seq = (char *)_mm_malloc(size, 64);
    index_alloc += size;
    assert_not_null(binary_ref_seq, size, index_alloc);
    char binary_ref_name[PATH_MAX];
    strcpy_s(binary_ref_name, PATH_MAX, prefix);
    strcat_s(binary_ref_name, PATH_MAX, ".0123");
    // sprintf(binary_ref_name, "%s.0123", prefix);
    std::fstream binary_ref_stream(binary_ref_name, std::ios::out | std::ios::binary);
    binary_ref_stream.seekg(0);
    fprintf(stderr, "init ticks = %llu\n", __rdtsc() - startTick);
    startTick = __rdtsc();
    int64_t i, count[16];
    memset(count, 0, sizeof(int64_t) * 16);
    for (i = 0; i < pac_len; i++)
    {
        switch (reference_seq[i])
        {
        case 'A':
            binary_ref_seq[i] = 0, ++count[0];
            break;
        case 'C':
            binary_ref_seq[i] = 1, ++count[1];
            break;
        case 'G':
            binary_ref_seq[i] = 2, ++count[2];
            break;
        case 'T':
            binary_ref_seq[i] = 3, ++count[3];
            break;
        default:
            binary_ref_seq[i] = 4;
        }
    }
    count[4] = count[0] + count[1] + count[2] + count[3];
    count[3] = count[0] + count[1] + count[2];
    count[2] = count[0] + count[1];
    count[1] = count[0];
    count[0] = 0;
    fprintf(stderr, "ref seq len = %ld\n", pac_len);
    binary_ref_stream.write(binary_ref_seq, pac_len * sizeof(char));
    fprintf(stderr, "binary seq ticks = %llu\n", __rdtsc() - startTick);
    startTick = __rdtsc();

    size = (pac_len + 2) * sizeof(int64_t);
    int64_t *suffix_array = (int64_t *)_mm_malloc(size, 64);
    index_alloc += size;
    assert_not_null(suffix_array, size, index_alloc);
    startTick = __rdtsc();
    // status = saisxx<const char *, int64_t *, int64_t>(reference_seq.c_str(), suffix_array + 1, pac_len, 4);
    status = saisxx(reference_seq.c_str(), suffix_array + 1, pac_len);
    suffix_array[0] = pac_len;
    fprintf(stderr, "build suffix-array ticks = %llu\n", __rdtsc() - startTick);
    startTick = __rdtsc();

    build_fm_index(prefix, binary_ref_seq, pac_len, suffix_array, count);
    fprintf(stderr, "build fm-index ticks = %llu\n", __rdtsc() - startTick);
    _mm_free(binary_ref_seq);
    _mm_free(suffix_array);
    return 0;
}

// int FMI_search::return_cp_occ_cp_count_at_index(int i, int j)
// {
//     pair<int,int> p= cp_occ_cp_count_map[make_pair(i,j)];
//     int cp_occ_cp_count_ciphertext_count=p.first;
//     int cp_occ_cp_count_index=p.second;

//     CT c=deserialize_ciphertext_from_file("ciphertext_cp_occ_cp_count_"+to_string(cp_occ_cp_count_ciphertext_count));
//     vector<int64_t> vv=decrypt_ciphertext_to_plaintext_vector(c);

//     return vv[j];
// }

// void FMI_search::set_cp_occ_cp_count_at_index(int i, int j, int val)
// {
//     pair<int,int> p= cp_occ_cp_count_map[make_pair(i,j)];
//     int cp_occ_cp_count_ciphertext_count=p.first;
//     int cp_occ_cp_count_index=p.second;

//     CT c=deserialize_ciphertext_from_file("ciphertext_cp_occ_cp_count_"+to_string(cp_occ_cp_count_ciphertext_count));
//     vector<int64_t> vv=decrypt_ciphertext_to_plaintext_vector(c);

//     vv[j]=val;
//     Plaintext plaintext = cc->MakePackedPlaintext(vv);
//     auto c = cc->Encrypt(kp.publicKey, plaintext);

//     if (!Serial::SerializeToFile("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/serializer_folder/ciphertext_cp_occ_cp_count_"+to_string(cp_occ_cp_count_ciphertext_count)+".txt" ,c, SerType::BINARY)) {
//         std::cerr
//             << "Error writing serialization of ciphertext_cp_occ_cp_count_"+to_string(cp_occ_cp_count_ciphertext_count)
//             << std::endl;
//         return;
//     }
// }

void FMI_search::load_index()
{
    auto load_index_time_start = high_resolution_clock::now();

    one_hot_mask_array = (uint64_t *)_mm_malloc(64 * sizeof(uint64_t), 64);
    one_hot_mask_array[0] = 0;
    uint64_t base = 0x8000000000000000L;
    one_hot_mask_array[1] = base;
    int64_t i = 0;
    for (i = 2; i < 64; i++)
    {
        one_hot_mask_array[i] = (one_hot_mask_array[i - 1] >> 1) | base;
    }
    /*
    one_hot_mask_array[0]=0; one_hot_mask_array[1]=100...0(63 zeroes); one_hot_mask_array[2]=11000..0(62 zeroes); .....; one_hot_mask_array[63]=1111...1(0 zeroes)
    */

    char *ref_file_name = file_name;
    // vecCT ref_file_name; strdup_enc(file_name,ref_file_name);
    // beCalls = 0;
    char cp_file_name[PATH_MAX];
    // vecCT cp_file_name; cp_file_name.resize(PATH_MAX);
    // strdup_enc(ref_file_name,cp_file_name);
    strcpy_s(cp_file_name, PATH_MAX, ref_file_name); // Now, cp_file_name is $INPUTS_DIR/fmi/broad

    // strcat_enc(cp_file_name,CP_FILENAME_SUFFIX,-1);
    strcat_s(cp_file_name, PATH_MAX, CP_FILENAME_SUFFIX); // Now, cp_file_name is $INPUTS_DIR/fmi/broad.bwt.2bit.64

    // printf("cp_file_name: %s\n",convert_ciphertext_vector_to_plaintext_string(cp_file_name));

    // Read the BWT and FM index of the reference sequence
    FILE *cpstream = NULL;
    cpstream = fopen(cp_file_name, "rb"); // cp_file_name is ../input-datasets/fmi/broad.bwt.2bit.64
    if (cpstream == NULL)
    {
        fprintf(stderr, "ERROR! Unable to open the file: %s\n", cp_file_name);
        exit(EXIT_FAILURE);
    }
    else
    {
        fprintf(stderr, "* Index file found. Loading index from %s\n", cp_file_name);
    }

    // err_fread_noeof(&reference_seq_len_read_plaintext, sizeof(int64_t), 1, cpstream);
    // reference_seq_len=encrypt_plaintext_integer_to_ciphertext(reference_seq_len_read_plaintext);
    // reference_seq_len_encrypted_bit_vector=get_encrypted_bits_vector(reference_seq_len_read_plaintext);

    // assert(decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(reference_seq_len, encrypt_plaintext_integer_to_ciphertext(0)))[0]>0);
    // assert(decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(reference_seq_len, encrypt_plaintext_integer_to_ciphertext(0x7fffffffffL)))[0]<=0);

    // fprintf(stderr, "* Reference seq len for bi-index = %ld\n", decrypt_ciphertext_to_plaintext_vector(reference_seq_len)[0]);

    // // create checkpointed occ
    // CT cp_occ_size_enc = cc->EvalAdd(shift_encrypted_bit_vector_and_return_integer(reference_seq_len_encrypted_bit_vector, CP_SHIFT), encrypt_plaintext_integer_to_ciphertext(1));  //CP_SHIFT=6
    // int64_t cp_occ_size = decrypt_ciphertext_to_plaintext_vector(cp_occ_size_enc)[0];
    // printf("cp_occ_size: %ld\n",cp_occ_size);

    err_fread_noeof(&reference_seq_len, sizeof(int64_t), 1, cpstream);
    reference_seq_len_enc = encrypt_plaintext_integer_to_ciphertext(reference_seq_len);
    reference_seq_len_encrypted_bit_vector = get_encrypted_bits_vector(reference_seq_len);

    // reference_seq_len: 6434693835; type(reference_seq_len): l
    // cout<<"reference_seq_len: "<<reference_seq_len<<"; type(reference_seq_len): "<<typeid(reference_seq_len).name()<<endl;
    // decrypt_ciphertext_to_plaintext_vector(reference_seq_len_enc)[0]: -744758
    // cout<<"decrypt_ciphertext_to_plaintext_vector(reference_seq_len_enc)[0]: "<<decrypt_ciphertext_to_plaintext_vector(reference_seq_len_enc)[0]<<endl;
    fprintf(stderr, "* Reference seq len for bi-index = %ld\n", decrypt_ciphertext_to_plaintext_vector(reference_seq_len_enc)[0]);

    assert(decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(reference_seq_len_enc, encode_integer_to_plaintext(0)))[0] > 0);
    // assert(decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(reference_seq_len_enc, encode_integer_to_plaintext(0x7fffffffffL)))[0]<=0);

    // create checkpointed occ
    CT cp_occ_size_enc = cc->EvalAdd(shift_encrypted_bit_vector_and_return_integer(reference_seq_len_encrypted_bit_vector, -1 * CP_SHIFT), encrypt_plaintext_integer_to_ciphertext(1)); // CP_SHIFT=6
    int64_t cp_occ_size_read = decrypt_ciphertext_to_plaintext_vector(cp_occ_size_enc)[0];                                                                                              // 10^8
    // printf("cp_occ_size_read: %ld\n", cp_occ_size_read);

    // fprintf(stderr, "* Reference seq len for bi-index = %ld\n", reference_seq_len);

    cp_occ = NULL;
    cp_occ_read = NULL;

    // int64_t count[5];
    // the following statement requests a 64-byte aligned memory block for 8 floating point elements-  farray = (float *)__mm_malloc(8*sizeof(float), 64);
    // https://stackoverflow.com/questions/3994035/what-is-aligned-memory-allocation
    // Reading into count
    err_fread_noeof(&count[0], sizeof(int64_t), 5, cpstream);
    for (int ii = 0; ii < 5; ii++)
    {
        count_enc[ii] = encrypt_plaintext_integer_to_ciphertext(count[ii]);
    }

    // cp_occ_size_read = reference_seq_len >> CP_SHIFT(6) + 1
    // if ((cp_occ_read = (CP_OCC_READ *)_mm_malloc(cp_occ_size_read * sizeof(CP_OCC_READ), 64)) == NULL) {
    //     fprintf(stderr, "ERROR! unable to allocated cp_occ_read memory\n");
    //     exit(EXIT_FAILURE);
    // }
    // malloc array gives access errors; even if one element is a CT or even if one element if a struct with a CT as member.  Calloc is error free.
    // Earlier, we used vector<CT> to replace all malloc(int/char). But here, it is like malloc(struct), so we can't convert one struct to CT.
    // So, we'll use calloc instead of vector<CT>.

    if ((cp_occ = (CP_OCC *)calloc(cp_occ_size_read, sizeof(CP_OCC))) == NULL)
    {
        fprintf(stderr, "ERROR! unable to allocated cp_occ memory\n");
        exit(EXIT_FAILURE);
    }

    // printf("Created cp_occ\n");

    // if ((cp_occ_read = (CP_OCC_READ *)calloc(cp_occ_size_read, sizeof(CP_OCC_READ))) == NULL)
    // {
    //     fprintf(stderr, "ERROR! unable to allocated cp_occ memory\n");
    //     exit(EXIT_FAILURE);
    // }

    // cp_occ = (CP_OCC*)(calloc(1,sizeof(CP_OCC)));
    // cp_occ->cp_count_enc.resize(tot_cp_occ_ct_count);
    // cp_occ->one_hot_bwt_str=(uint64_t*)(calloc(tot_cp_occ_ct_count,sizeof(uint64_t)));

    // cout<<"Reading into cp_occ_read!\n";

    err_fread_noeof(cp_occ, sizeof(CP_OCC), cp_occ_size_read, cpstream);

    //     FILE* fp2 = fopen("cp_occ_cp_count_file.txt", "w");
    //     if (fp2 == NULL) {
    //         printf("Failed to open file\n");
    //         return;
    //     }
    //     for (int k = 0; k < cp_occ_size_read; k++) {
    //         fprintf(fp2, "%d\n", cp_occ[k].cp_count[0]);
    //         fprintf(fp2, "%d\n", cp_occ[k].cp_count[1]);
    //         fprintf(fp2, "%d\n", cp_occ[k].cp_count[2]);
    //         fprintf(fp2, "%d\n\n", cp_occ[k].cp_count[3]);
    //     }
    //    fclose(fp2);

    //    fp2 = fopen("cp_occ_one_hot_bwt_str_file.txt", "w");
    //     if (fp2 == NULL) {
    //         printf("Failed to open file\n");
    //         return;
    //     }
    //     for (int k = 0; k < cp_occ_size_read; k++) {
    //         fprintf(fp2, "%d\n", cp_occ[k].one_hot_bwt_str[0]);
    //         fprintf(fp2, "%d\n", cp_occ[k].one_hot_bwt_str[1]);
    //         fprintf(fp2, "%d\n", cp_occ[k].one_hot_bwt_str[2]);
    //         fprintf(fp2, "%d\n\n", cp_occ[k].one_hot_bwt_str[3]);
    //     }
    //    fclose(fp2);

    // printf("Read into cp_occ\n");

    // vector<int64_t> vv;
    // for(int i=0;i<cp_occ_size_read;i++)
    // {
    //     for(int j=0;j<4;j++)
    //     {
    //         int64_t zz=cp_occ[i].cp_count[j];
    //         cout<<"i: "<<i<<"; j: "<<j<<"; zz: "<<zz<<endl;
    //         if(zz<p)
    //         {
    //             if(cp_occ_cp_count_index==16384)                      // we are going to write to cp_occ_cp_count_index; but if it is 16384, then vv is already full and we need to encrypt vv to CT and write current in fresh vv
    //             {
    //                 if(vv.size()!=0)
    //                 {
    //                     CT c=encrypt_plaintext_vector_to_ciphertext(vv);
    //                     if (!Serial::SerializeToFile("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/serializer_folder/ciphertext_cp_occ_cp_count_"+to_string(cp_occ_cp_count_ciphertext_count)+".txt" ,c, SerType::BINARY)) {
    //                         std::cerr
    //                             << "Error writing serialization of ciphertext_cp_occ_cp_count_"+to_string(cp_occ_cp_count_ciphertext_count)
    //                             << std::endl;
    //                         return;
    //                     }
    //                     cp_occ_cp_count_ciphertext_count++;
    //                 }

    //                 vv.resize(16384,0);
    //                 cp_occ_cp_count_index=0;
    //             }

    //             if(vv.size()==0)
    //                 vv.resize(16384,0);

    //             vv[cp_occ_cp_count_index]=zz;
    //             // cp_occ_cp_count_map.insert(make_pair(make_pair(i,j), make_pair(make_pair(cp_occ_cp_count_ciphertext_count,cp_occ_cp_count_index),make_pair(cp_occ_cp_count_ciphertext_count,cp_occ_cp_count_index)))); // cp_occ_cp_count_ciphertext_count starts from 0
    //             cp_occ_cp_count_map.insert(make_pair(make_pair(i,j), make_pair(cp_occ_cp_count_ciphertext_count,cp_occ_cp_count_index))); // cp_occ_cp_count_ciphertext_count starts from 0

    //             cout<<"Single- cp_occ_cp_count_ciphertext_count: "<<cp_occ_cp_count_ciphertext_count<<"; cp_occ_cp_count_index: "<<cp_occ_cp_count_index<<endl;

    //             cp_occ_cp_count_index++;
    //         }
    //         else
    //         {
    //             int64_t index1, index2, ct_count1, ct_count2;

    //             while(zz>0)
    //             {
    //                 if(cp_occ_cp_count_index==16384)
    //                 {
    //                     if(vv.size()!=0)
    //                     {
    //                         CT c=encrypt_plaintext_vector_to_ciphertext(vv);
    //                         if (!Serial::SerializeToFile("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/serializer_folder/ciphertext_cp_occ_cp_count_"+to_string(cp_occ_cp_count_ciphertext_count)+".txt" ,c, SerType::BINARY)) {
    //                             std::cerr
    //                                 << "Error writing serialization of ciphertext_cp_occ_cp_count_"+to_string(cp_occ_cp_count_ciphertext_count)
    //                                 << std::endl;
    //                             return;
    //                         }
    //                         cp_occ_cp_count_ciphertext_count++;
    //                     }

    //                     vv.resize(16384,0);
    //                     cp_occ_cp_count_index=0;
    //                 }
    //                 if(vv.size()==0)
    //                     vv.resize(16384,0);

    //                 vv[cp_occ_cp_count_index]=zz%p; // every time we write to vv[cp_occ_cp_count_index], we make sure that cp_occ_cp_count_index is not 16384; if it is, we create CT out of it
    //                 zz=(zz-zz%p)/p;
    //                 cout<<"zz2: "<<zz<<endl;

    //                 if(zz!=0)
    //                 {
    //                     index1=cp_occ_cp_count_index; ct_count1=cp_occ_cp_count_ciphertext_count;
    //                 }
    //                 else{
    //                     index2=cp_occ_cp_count_index; ct_count2=cp_occ_cp_count_ciphertext_count;
    //                 }

    //                 cout<<"cp_occ_cp_count_index: "<<cp_occ_cp_count_index<<"; vv[cp_occ_cp_count_index]= "<<vv[cp_occ_cp_count_index]<<endl;
    //                 cp_occ_cp_count_index++;

    //             }

    //             // as zz was too big, so 2 components of that have been formed and written to vv; note that the above while-loop will run for only 2 iterations
    //             // this means that cp_occ[i].cp_count[j] is represented by decrypt_ciphertext_to_vector(ciphertext_cp_occ_cp_count_ct_count1)[index1] and decrypt_ciphertext_to_vector(ciphertext_cp_occ_cp_count_ct_count2)[index2]
    //             cp_occ_cp_count_map.insert(make_pair(make_pair(i,j), make_pair(make_pair(ct_count1,index1),make_pair(ct_count2,index2)))); // cp_occ_cp_count_ciphertext_count starts from 0
    //             cout<<"Two times- ct_count1: "<<ct_count1<<"; index1: "<<index1<<"; ct_count2: "<<ct_count2<<"; index2: "<<index2<<endl;
    //         }

    //     }
    // }

    int64_t tot_cp_occ_ct_count;
    tot_cp_occ_ct_count = cp_occ_size_read / 4096;
    if (cp_occ_size_read % 4096 == 0)
        tot_cp_occ_ct_count++;

    // for(int64_t i=0;i<tot_cp_occ_ct_count;i++)
    // {
    //     vector<int64_t> vv(16384); // uint64_t one_hot_bwt_str_copy[16384];
    //     int j=0;

    //     for(j=0;j<16384 && i*16384+j<4*cp_occ_size_read;j++)            // j goes till 16384; index i in cp_occ[i] goes till 4095;
    //     {
    //         cout<<i<<"*16384+"<<j<<" = "<<i*16384+j<<"; "<<endl;
    //         int64_t zz = cp_occ[(i*16384+j)/4].cp_count[(i*16384+j)%4];
    //         cout<<"zz: "<<zz<<endl;
    //         if(zz<p)
    //         {
    //             vv[j]=zz;
    //         }
    //         else
    //         {
    //             vv[j]=-1*(p-1);
    //             vector<int64_t> vv_2(16384,0); int k=0;
    //             while(zz>0)
    //             {
    //                 vv_2[k]=zz%p;
    //                 zz=(zz-zz%p)/p;
    //                 cout<<"k: "<<k<<"; vv_2[k]= "<<vv_2[k]<<endl;
    //                 k++;
    //             }
    //             CT c2=encrypt_plaintext_vector_to_ciphertext(vv_2);

    //             if (!Serial::SerializeToFile("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/serializer_folder/ciphertext_cp_occ_cp_count_"+to_string(i)+"_"+to_string(j)+".txt" ,c2, SerType::BINARY)) {
    //                 std::cerr
    //                     << "Error writing serialization of ciphertext_cp_occ_"+to_string(i)+"_"+to_string(j)
    //                     << std::endl;
    //                 return;
    //             }
    //             cout<<endl;
    //         }
    //         // cout<<"vv["<<j<<"]= "<<vv[j]<<endl;
    //         // cp_occ->one_hot_bwt_str[j]=cp_occ_read[(i*16384+j)/4].one_hot_bwt_str[(i*16384+j)%4];
    //     }
    //     // PCT p1=make_pair(make_pair(encrypt_plaintext_vector_to_ciphertext(vv),0),j+1);
    //     CT c1=encrypt_plaintext_vector_to_ciphertext(vv);
    //     if (!Serial::SerializeToFile("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/serializer_folder/ciphertext_cp_occ_cp_count_"+to_string(i)+".txt" ,c1, SerType::BINARY)) {
    //         std::cerr
    //             << "Error writing serialization of ciphertext_cp_occ_"+to_string(i)
    //             << std::endl;
    //         return;
    //     }
    //     cout<<endl;
    // }

    int64_t ii = 0;
    for (ii = 0; ii < 5; ii++) // update read count structure
    {
        count[ii] = count[ii] + 1;
        count_enc[ii] = cc->EvalAdd(count_enc[ii], encode_integer_to_plaintext(1));
        // cout << "decrypt_ciphertext_to_plaintext_vector(count_enc[" << ii << "])[0]: " << decrypt_ciphertext_to_plaintext_vector(count_enc[ii])[0] << endl;
    }

    // #define SA_COMPRESSION 1
    // #define SA_COMPX 03 // (= power of 2)
    // #define SA_COMPX_MASK 0x7    // 0x7 or 0x3 or 0x1

#if SA_COMPRESSION

    int64_t reference_seq_len_ = (reference_seq_len >> SA_COMPX) + 1;
    // cout << "reference_seq_len_: " << reference_seq_len_ << endl; // 8*10^8
    CT reference_seq_len_enc_ = cc->EvalAdd(shift_encrypted_bit_vector_and_return_integer(reference_seq_len_encrypted_bit_vector, -1 * SA_COMPX), encrypt_plaintext_integer_to_ciphertext(1));
    // cout << "decrypt_ciphertext_to_plaintext_vector(reference_seq_len_enc_)[0]: " << decrypt_ciphertext_to_plaintext_vector(reference_seq_len_enc_)[0] << endl;

    sa_ms_byte = (int8_t *)_mm_malloc(reference_seq_len_ * sizeof(int8_t), 64);
    sa_ls_word = (uint32_t *)_mm_malloc(reference_seq_len_ * sizeof(uint32_t), 64);
    // sa_ms_byte_enc.resize(reference_seq_len_); // sa_ms_byte_enc.shrink_to_fit();
    // sa_ls_word_enc.resize(reference_seq_len_); // sa_ls_word_enc.shrink_to_fit();
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
    // int sz=0;
    // while(decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(reference_seq_len_enc, encode_integer_to_plaintext(sz)))[0]!=0)
    // {
    //     sa_ms_byte_enc.resize(reference_seq_len_enc_,0);
    // sa_ls_word_enc.resize(reference_seq_len_enc_,0);
    //     sz++;
    // }

    // While decrypting, sa_ms_byte, type-cast it to int8_t. For sa_ls_word, type-cast it to uint32_t

    err_fread_noeof(sa_ms_byte, sizeof(int8_t), reference_seq_len_, cpstream);
    err_fread_noeof(sa_ls_word, sizeof(uint32_t), reference_seq_len_, cpstream);

    //     FILE *fp;
    //     fp = fopen("sa_ms_byte_file.txt", "w");
    //     if (fp == NULL) {
    //         printf("Failed to open file\n");
    //         return;
    //     }
    //     for (int k = 0; k < reference_seq_len_; k++) {
    //         fprintf(fp, "%d\n", sa_ms_byte[k]);
    //     }
    //    fclose(fp);

    //    fp = fopen("sa_ls_word_file.txt", "w");
    //     if (fp == NULL) {
    //         printf("Failed to open file\n");
    //         return;
    //     }
    //     for (int k = 0; k < reference_seq_len_; k++) {
    //         fprintf(fp, "%d\n", sa_ls_word[k]);
    //     }
    //    fclose(fp);

    // for(int i=0;i<reference_seq_len_;i++)
    // {
    //     cout<<"sa_ls_word["<<i<<"] = "<<sa_ls_word[i]<<endl;
    //     cout<<"sa_ms_byte["<<i<<"] = "<<sa_ms_byte[i]<<endl;
    //     sa_ms_byte_enc[i]=encrypt_plaintext_integer_to_ciphertext(sa_ms_byte[i]);
    //     sa_ls_word_enc[i]=encrypt_plaintext_integer_to_ciphertext(sa_ls_word[i]);
    // }

#else

    sa_ms_byte = (int8_t *)_mm_malloc(reference_seq_len * sizeof(int8_t), 64);
    sa_ls_word = (uint32_t *)_mm_malloc(reference_seq_len * sizeof(uint32_t), 64);
    // sa_ms_byte_enc.resize(reference_seq_len, 0);
    // sa_ls_word_enc.resize(reference_seq_len, 0);

    err_fread_noeof(sa_ms_byte, sizeof(int8_t), reference_seq_len, cpstream);
    err_fread_noeof(sa_ls_word, sizeof(uint32_t), reference_seq_len, cpstream);
    // for (int i = 0; i < reference_seq_len; i++)
    // {
    //     sa_ms_byte_enc[i] = encrypt_plaintext_integer_to_ciphertext(sa_ms_byte[i]);
    //     sa_ls_word_enc[i] = encrypt_plaintext_integer_to_ciphertext(sa_ls_word[i]);
    // }

#endif

    sentinel_index = -1;
    sentinel_index_enc = encrypt_plaintext_integer_to_ciphertext(-1);

#if SA_COMPRESSION
    err_fread_noeof(&sentinel_index, sizeof(int64_t), 1, cpstream);
    // fprintf(stderr, "* sentinel-index: %ld\n", sentinel_index);
    sentinel_index_enc = encrypt_plaintext_integer_to_ciphertext(sentinel_index);

    // fprintf(stderr, "* sentinel-index: %ld\n", sentinel_index);
    fprintf(stderr, "* sentinel-index-enc: %ld\n", decrypt_ciphertext_to_plaintext_vector(sentinel_index_enc)[0]);

#endif
    fclose(cpstream);

    int64_t x;
#if !SA_COMPRESSION
    for (x = 0; x < reference_seq_len; x++) // x is a position
    {
// fprintf(stderr, "x: %ld\n", x);
#if SA_COMPRESSION
        if (decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(get_sa_entry_compressed(x), encode_integer_to_plaintext(0)))[0] == 0)
        {
            sentinel_index = x;
            sentinel_index_enc = encrypt_plaintext_integer_to_ciphertext(x);
            break;
        }
#else
        if (decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(get_sa_entry(x), encode_integer_to_plaintext(0)))[0] == 0)
        {
            sentinel_index = x;
            sentinel_index_enc = encrypt_plaintext_integer_to_ciphertext(x);
            break;
        }
#endif
    }
    fprintf(stderr, "\nsentinel_index: %ld\n", x);
    fprintf(stderr, "* sentinel-index-enc: %ld\n", decrypt_ciphertext_to_plaintext_vector(sentinel_index_enc)[0]);
#endif

    // fprintf(stderr, "* Count:\n");
    // for (x = 0; x < 5; x++)
    // {
    //     fprintf(stderr, "%ld,\t%lu\n", x, (unsigned long)count[x]);
    // }
    fprintf(stderr, "* count_enc:\n");
    for (x = 0; x < 5; x++)
    {
        fprintf(stderr, "%ld,\t%lu\n", x, (unsigned long)decrypt_ciphertext_to_plaintext_vector(count_enc[x])[0]);
    }
    fprintf(stderr, "\n");

    // ref_file_name is  ../../input-datasets/fmi/broad
    fprintf(stderr, "* Reading other elements of the index from files %s\n",
            ref_file_name);

    auto load_index_time_end = high_resolution_clock::now();
    auto duration_load_index_before_bwa_idx_load_ele = duration_cast<microseconds>(load_index_time_end - load_index_time_start);

    cout << "Time taken by load_index_before_bwa_idx_load_ele: "<<duration_load_index_before_bwa_idx_load_ele.count() <<" microseconds" << endl;

    // In bwa.h, #define BWA_IDX_ALL 0x7
    bwa_idx_load_ele(ref_file_name, BWA_IDX_ALL);
    // printf("Out of bwa_idx_load_ele!\n");

    fprintf(stderr, "* Done reading Index!!\n");
}

void FMI_search::getSMEMsOnePosOneThread(uint8_t *enc_qdb,
                                         int16_t *query_pos_array,
                                         int32_t *min_intv_array,
                                         int32_t *rid_array,
                                         int32_t numReads,
                                         int32_t batch_size,
                                         const bseq1_t *seq_,
                                         int32_t *query_cum_len_ar,
                                         int32_t max_readlength,
                                         int32_t minSeedLen,
                                         SMEM *matchArray,
                                         int64_t *__numTotalSmem)
{
    // printf("In getSMEMsOnePosOneThread\n");
    int64_t numTotalSmem = *__numTotalSmem;
    SMEM prevArray[max_readlength];

    uint32_t i;
    // Perform SMEM for original reads
    for (i = 0; i < numReads; i++)
    {
        // printf("i: %d\n", i);
        int x = query_pos_array[i];
        int32_t rid = rid_array[i];
        int next_x = x + 1;

        int readlength = seq_[rid].l_seq;
        int offset = query_cum_len_ar[rid];
        // uint8_t a = enc_qdb[rid * readlength + x];
        uint8_t a = enc_qdb[offset + x];
        // printf("offset+x: %d, a: %d\n", offset + x, a);

        if (a < 4)
        {
            SMEM smem;
            smem.rid = rid;
            smem.m = x;
            smem.n = x;
            smem.k = decrypt_ciphertext_to_plaintext_vector(count_enc[a])[0];
            smem.l = decrypt_ciphertext_to_plaintext_vector(count_enc[3 - a])[0];
            smem.s = decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(count_enc[a + 1], count_enc[a]))[0];
            int numPrev = 0;

            // printf("x: %d, readlength: %d, rid: %d\n", x, readlength, rid);
            // printf("smem.m: %d, smem.n: %d, smem.k: %d, smem.l: %d, smem.s: %d\n", smem.m, smem.n, smem.k, smem.l, smem.s);

            int j;
            // printf("x: %d, readlength: %d, rid: %d\n", x, readlength, rid);

            for (j = x + 1; j < readlength; j++)
            {
                // a = enc_qdb[rid * readlength + j];
                a = enc_qdb[offset + j];
                // printf("j: %d, a: %d\n", j, a);
                next_x = j + 1;
                if (a < 4)
                {
                    // printf("min_intv_array[i]: %d\n", min_intv_array[i]);
                    SMEM smem_ = smem;

                    // Forward extension is backward extension with the BWT of reverse complement
                    smem_.k = smem.l;
                    smem_.l = smem.k;

                    // printf("smem_.k: %d, smem_.l: %d\n", smem_.k, smem_.l);
                    SMEM newSmem_ = backwardExt(smem_, 3 - a);
                    // SMEM newSmem_ = forwardExt(smem_, 3 - a);
                    SMEM newSmem = newSmem_;
                    newSmem.k = newSmem_.l;
                    newSmem.l = newSmem_.k;
                    newSmem.n = j;

                    // printf("newSmem.k: %d, newSmem.l: %d, newSmem.n: %d, newSmem.s: %d\n", newSmem.k, newSmem.l, newSmem.n, newSmem.s);

                    int32_t s_neq_mask = newSmem.s != smem.s;

                    prevArray[numPrev] = smem;
                    numPrev += s_neq_mask;
                    if (newSmem.s < min_intv_array[i])
                    {
                        next_x = j;
                        // printf("newSmem.s: %d, min_intv_array[%d]: %d\n", newSmem.s, i, min_intv_array[i]);
                        break;
                    }
                    smem = newSmem;
                    // #ifdef ENABLE_PREFETCH
                    //                     _mm_prefetch((const CT *)( &cp_occ[(smem.k) >> CP_SHIFT]), _MM_HINT_T0);
                    //                     _mm_prefetch((const CT *)(&cp_occ[(smem.l) >> CP_SHIFT]), _MM_HINT_T0);

                    // #endif
                }
                else
                {
                    // printf("a>4; break!\n");
                    break;
                }
            }
            if (smem.s >= min_intv_array[i])
            {

                prevArray[numPrev] = smem;
                numPrev++;
            }

            SMEM *prev;
            prev = prevArray;

            int p;
            for (p = 0; p < (numPrev / 2); p++)
            {
                SMEM temp = prev[p];
                prev[p] = prev[numPrev - p - 1];
                prev[numPrev - p - 1] = temp;
            }

            // Backward search
            int cur_j = readlength;
            for (j = x - 1; j >= 0; j--)
            {
                //ncout << "j2: " << j << endl;
                int numCurr = 0;
                int curr_s = -1;
                // a = enc_qdb[rid * readlength + j];
                a = enc_qdb[offset + j];

                if (a > 3)
                {
                    break;
                }
                for (p = 0; p < numPrev; p++)
                {
                    SMEM smem = prev[p];
                    SMEM newSmem = backwardExt(smem, a);
                    newSmem.m = j;

                    if ((newSmem.s < min_intv_array[i]) && ((smem.n - smem.m + 1) >= minSeedLen))
                    {
                        cur_j = j;

                        matchArray[numTotalSmem++] = smem;

                        // printf("1 numTotalSmem after matchArray[numTotalSmem++]: %d\n", numTotalSmem);
                        break;
                    }
                    if ((newSmem.s >= min_intv_array[i]) && (newSmem.s != curr_s))
                    {
                        curr_s = newSmem.s;
                        prev[numCurr++] = newSmem;
                        // #ifdef ENABLE_PREFETCH
                        //                         _mm_prefetch((const CT *)(&cp_occ[(newSmem.k) >> CP_SHIFT]), _MM_HINT_T0);
                        //                         _mm_prefetch((const CT *)(&cp_occ[(newSmem.k + newSmem.s) >> CP_SHIFT]), _MM_HINT_T0);
                        // #endif
                        break;
                    }
                }
                p++;
                for (; p < numPrev; p++)
                {
                    SMEM smem = prev[p];

                    SMEM newSmem = backwardExt(smem, a);
                    newSmem.m = j;

                    if ((newSmem.s >= min_intv_array[i]) && (newSmem.s != curr_s))
                    {
                        curr_s = newSmem.s;
                        prev[numCurr++] = newSmem;
                        // #ifdef ENABLE_PREFETCH
                        //                         _mm_prefetch((const CT *)(&cp_occ[(newSmem.k) >> CP_SHIFT]), _MM_HINT_T0);
                        //                         _mm_prefetch((const CT *)(&cp_occ[(newSmem.k + newSmem.s) >> CP_SHIFT]), _MM_HINT_T0);
                        // #endif
                    }
                }
                numPrev = numCurr;
                if (numCurr == 0)
                {
                    break;
                }
            }
            if (numPrev != 0)
            {
                SMEM smem = prev[0];
                if (((smem.n - smem.m + 1) >= minSeedLen))
                {

                    matchArray[numTotalSmem++] = smem;
                    // printf("2 numTotalSmem after matchArray[numTotalSmem++]: %d\n", numTotalSmem);
                }
                numPrev = 0;
            }
        }
        query_pos_array[i] = next_x;
    }
    (*__numTotalSmem) = numTotalSmem;
}

void FMI_search::getSMEMsAllPosOneThread(uint8_t *enc_qdb,
                                         int32_t *min_intv_array,
                                         int32_t *rid_array,
                                         int32_t numReads,
                                         int32_t batch_size,
                                         const bseq1_t *seq_,
                                         int32_t *query_cum_len_ar,
                                         int32_t max_readlength,
                                         int32_t minSeedLen,
                                         SMEM *matchArray,
                                         int64_t *__numTotalSmem)
{
    // printf("In getSMEMsAllPosOneThread!\n");
    int16_t *query_pos_array = (int16_t *)_mm_malloc(numReads * sizeof(int16_t), 64);

    int32_t i;
    for (i = 0; i < numReads; i++)
        query_pos_array[i] = 0;

    int32_t numActive = numReads;
    (*__numTotalSmem) = 0;

    do
    {
        int32_t head = 0;
        int32_t tail = 0;
        // printf("numActive: %d\n", numActive);
        for (head = 0; head < numActive; head++)
        {
            // printf("head: %d\n", head);
            int readlength = seq_[rid_array[head]].l_seq;
            // printf("query_pos_array[%d]:%d ;readlength: %d\n", head, query_pos_array[head], readlength);
            if (query_pos_array[head] < readlength)
            {
                rid_array[tail] = rid_array[head];
                query_pos_array[tail] = query_pos_array[head];
                min_intv_array[tail] = min_intv_array[head];
                // printf("min_intv_array[%d]= %d\n", tail, min_intv_array[tail]);

                tail++;
            }
        }
        // printf("tail: %d\n", tail);
        getSMEMsOnePosOneThread(enc_qdb,
                                query_pos_array,
                                min_intv_array,
                                rid_array,
                                tail,
                                batch_size,
                                seq_,
                                query_cum_len_ar,
                                max_readlength,
                                minSeedLen,
                                matchArray,
                                __numTotalSmem);
        numActive = tail;
    } while (numActive > 0);

    _mm_free(query_pos_array);
}

int64_t FMI_search::bwtSeedStrategyAllPosOneThread(uint8_t *enc_qdb,
                                                   int32_t *max_intv_array,
                                                   int32_t numReads,
                                                   const bseq1_t *seq_,
                                                   int32_t *query_cum_len_ar,
                                                   int32_t minSeedLen,
                                                   SMEM *matchArray)
{
    int32_t i;

    int64_t numTotalSeed = 0;

    for (i = 0; i < numReads; i++)
    {
        int readlength = seq_[i].l_seq;
        int16_t x = 0;
        while (x < readlength)
        {
            int next_x = x + 1;

            // Forward search
            SMEM smem;
            smem.rid = i;
            smem.m = x;
            smem.n = x;

            int offset = query_cum_len_ar[i];
            uint8_t a = enc_qdb[offset + x];
            // uint8_t a = enc_qdb[i * readlength + x];

            if (a < 4)
            {
                smem.k = decrypt_ciphertext_to_plaintext_vector(count_enc[a])[0];
                smem.l = decrypt_ciphertext_to_plaintext_vector(count_enc[3 - a])[0];
                smem.s = decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(count_enc[a + 1], count_enc[a]))[0];

                int j;
                for (j = x + 1; j < readlength; j++)
                {
                    next_x = j + 1;
                    // a = enc_qdb[i * readlength + j];
                    a = enc_qdb[offset + j];
                    if (a < 4)
                    {
                        SMEM smem_ = smem;

                        // Forward extension is backward extension with the BWT of reverse complement
                        smem_.k = smem.l;
                        smem_.l = smem.k;
                        SMEM newSmem_ = backwardExt(smem_, 3 - a);
                        // SMEM smem = backwardExt(smem, 3 - a);
                        // smem.n = j;
                        SMEM newSmem = newSmem_;
                        newSmem.k = newSmem_.l;
                        newSmem.l = newSmem_.k;
                        newSmem.n = j;
                        smem = newSmem;
                        // #ifdef ENABLE_PREFETCH
                        //                         _mm_prefetch((const CT *)(&cp_occ[(smem.k) >> CP_SHIFT]), _MM_HINT_T0);
                        //                         _mm_prefetch((const CT *)(&cp_occ[(smem.l) >> CP_SHIFT]), _MM_HINT_T0);
                        // #endif

                        if ((smem.s < max_intv_array[i]) && ((smem.n - smem.m + 1) >= minSeedLen))
                        {

                            if (smem.s > 0)
                            {
                                matchArray[numTotalSeed++] = smem;
                            }
                            break;
                        }
                    }
                    else
                    {

                        break;
                    }
                }
            }
            x = next_x;
        }
    }
    return numTotalSeed;
}

void FMI_search::getSMEMs(uint8_t *enc_qdb,
                          int32_t numReads,
                          int32_t batch_size,
                          int32_t readlength,
                          int32_t minSeedLen,
                          int32_t nthreads,
                          SMEM *matchArray,
                          int64_t *numTotalSmem)
{
    SMEM *prevArray = (SMEM *)_mm_malloc(nthreads * readlength * sizeof(SMEM), 64);
    SMEM *currArray = (SMEM *)_mm_malloc(nthreads * readlength * sizeof(SMEM), 64);

    // #pragma omp parallel num_threads(nthreads)
    {
        int tid = 0; // omp_get_thread_num();   // removed omp
        numTotalSmem[tid] = 0;
        SMEM *myPrevArray = prevArray + tid * readlength;
        SMEM *myCurrArray = prevArray + tid * readlength;

        int32_t perThreadQuota = (numReads + (nthreads - 1)) / nthreads;
        int32_t first = tid * perThreadQuota;
        int32_t last = (tid + 1) * perThreadQuota;
        if (last > numReads)
            last = numReads;
        SMEM *myMatchArray = matchArray + first * readlength;

        uint32_t i;
        // Perform SMEM for original reads
        for (i = first; i < last; i++)
        {
            int x = readlength - 1;
            int numPrev = 0;
            int numSmem = 0;

            while (x >= 0)
            {
                // Forward search
                SMEM smem;
                smem.rid = i;
                smem.m = x;
                smem.n = x;
                uint8_t a = enc_qdb[i * readlength + x];

                if (a > 3)
                {
                    x--;
                    continue;
                }
                smem.k = count[a];
                smem.l = count[3 - a];
                smem.s = count[a + 1] - count[a];

                int j;
                for (j = x + 1; j < readlength; j++)
                {
                    a = enc_qdb[i * readlength + j];
                    if (a < 4)
                    {
                        SMEM smem_ = smem;

                        // Forward extension is backward extension with the BWT of reverse complement
                        smem_.k = smem.l;
                        smem_.l = smem.k;
                        SMEM newSmem_ = backwardExt(smem_, 3 - a);
                        SMEM newSmem = newSmem_;
                        newSmem.k = newSmem_.l;
                        newSmem.l = newSmem_.k;
                        newSmem.n = j;

                        if (newSmem.s != smem.s)
                        {
                            myPrevArray[numPrev] = smem;
                            numPrev++;
                        }
                        smem = newSmem;
                        if (newSmem.s == 0)
                        {
                            break;
                        }
                    }
                    else
                    {
                        myPrevArray[numPrev] = smem;
                        numPrev++;
                        break;
                    }
                }
                if (smem.s != 0)
                {
                    myPrevArray[numPrev++] = smem;
                }

                SMEM *curr, *prev;
                prev = myPrevArray;
                curr = myCurrArray;

                int p;
                for (p = 0; p < (numPrev / 2); p++)
                {
                    SMEM temp = prev[p];
                    prev[p] = prev[numPrev - p - 1];
                    prev[numPrev - p - 1] = temp;
                }

                int next_x = x - 1;

                // Backward search
                int cur_j = readlength;
                for (j = x - 1; j >= 0; j--)
                {
                    int numCurr = 0;
                    int curr_s = -1;
                    a = enc_qdb[i * readlength + j];
                    // printf("a = %d\n", a);
                    if (a > 3)
                    {
                        next_x = j - 1;
                        break;
                    }
                    for (p = 0; p < numPrev; p++)
                    {
                        SMEM smem = prev[p];
                        SMEM newSmem = backwardExt(smem, a);
                        newSmem.m = j;

                        if (newSmem.s == 0)
                        {
                            if ((numCurr == 0) && (j < cur_j))
                            {
                                cur_j = j;
                                if ((smem.n - smem.m + 1) >= minSeedLen)
                                    myMatchArray[numTotalSmem[tid] + numSmem++] = smem;
                            }
                        }
                        if ((newSmem.s != 0) && (newSmem.s != curr_s))
                        {
                            curr_s = newSmem.s;
                            curr[numCurr++] = newSmem;
                        }
                    }
                    SMEM *temp = prev;
                    prev = curr;
                    curr = temp;
                    numPrev = numCurr;
                    if (numCurr == 0)
                    {
                        next_x = j;
                        break;
                    }
                    else
                    {
                        next_x = j - 1;
                    }
                }
                if (numPrev != 0)
                {
                    SMEM smem = prev[0];
                    if ((smem.n - smem.m + 1) >= minSeedLen)
                        myMatchArray[numTotalSmem[tid] + numSmem++] = smem;
                    numPrev = 0;
                }
                x = next_x;
            }
            numTotalSmem[tid] += numSmem;
        }
    }

    _mm_free(prevArray);
    _mm_free(currArray);
}

int compare_smem(const void *a, const void *b)
{
    SMEM *pa = (SMEM *)a;
    SMEM *pb = (SMEM *)b;

    if (pa->rid < pb->rid)
        return -1;
    if (pa->rid > pb->rid)
        return 1;

    if (pa->m < pb->m)
        return -1;
    if (pa->m > pb->m)
        return 1;
    if (pa->n > pb->n)
        return -1;
    if (pa->n < pb->n)
        return 1;
    return 0;
}

void FMI_search::sortSMEMs(SMEM *matchArray,
                           int64_t numTotalSmem[],
                           int32_t numReads,
                           int32_t readlength,
                           int nthreads)
{
    int tid;
    int32_t perThreadQuota = (numReads + (nthreads - 1)) / nthreads;
    for (tid = 0; tid < nthreads; tid++)
    {
        int32_t first = tid * perThreadQuota;
        SMEM *myMatchArray = matchArray + first * readlength;
        qsort(myMatchArray, numTotalSmem[tid], sizeof(SMEM), compare_smem);
    }
}

SMEM FMI_search::backwardExt(SMEM smem, uint8_t a)
{
    // beCalls++;
    uint8_t b;

    int64_t k[4], l[4], s[4];
    for (b = 0; b < 4; b++)
    {
        int64_t sp = (int64_t)(smem.k);
        int64_t ep = (int64_t)(smem.k) + (int64_t)(smem.s);
        GET_OCC(sp, b, occ_id_sp, y_sp, occ_sp, one_hot_bwt_str_c_sp, match_mask_sp);
        GET_OCC(ep, b, occ_id_ep, y_ep, occ_ep, one_hot_bwt_str_c_ep, match_mask_ep);
        // k[b] = count[b] + occ_sp;
        k[b] = decrypt_ciphertext_to_plaintext_vector(count_enc[b])[0]+occ_sp;
        s[b] = occ_ep - occ_sp;
    }

    int64_t sentinel_offset = 0;
    if ((smem.k <= decrypt_ciphertext_to_plaintext_vector(sentinel_index_enc)[0]) && ((smem.k + smem.s) > decrypt_ciphertext_to_plaintext_vector(sentinel_index_enc)[0]))
        sentinel_offset = 1;
    l[3] = smem.l + sentinel_offset;
    l[2] = l[3] + s[3];
    l[1] = l[2] + s[2];
    l[0] = l[1] + s[1];

    smem.k = k[a];
    smem.l = l[a];
    smem.s = s[a];
    return smem;
}

// sa_entry = sa_ms_byte[pos]*2^32 + sa_ls_word[pos]
int64_t FMI_search::get_sa_entry(int64_t pos)
{
    int64_t sa_entry = sa_ms_byte[pos];
    // CT sa_entry_enc = sa_ms_byte_enc[pos];      // static_cast<int16_t>(sa_ms_byte[i]) values are either 0/1 ;  1 << 32 =1*2^32 < p/2 (p=12869861377)

    sa_entry = sa_entry << 32;
    // sa_entry_enc=shift_left(sa_entry_enc,32);

    sa_entry = sa_entry + sa_ls_word[pos];
    // sa_entry_enc = cc->EvalAdd(sa_entry_enc,sa_ls_word_enc[pos]);
    return sa_entry;
}

void FMI_search::get_sa_entries(int64_t *posArray, int64_t *coordArray, uint32_t count, int32_t nthreads)
{
    uint32_t i;
    // #pragma omp parallel for num_threads(nthreads)
    for (i = 0; i < count; i++)
    {
        int64_t pos = posArray[i];
        int64_t sa_entry = sa_ms_byte[pos];
        sa_entry = sa_entry << 32;
        sa_entry = sa_entry + sa_ls_word[pos];
        //_mm_prefetch((const char *)(sa_ms_byte + pos + SAL_PFD), _MM_HINT_T0);
        coordArray[i] = sa_entry;
    }
}

void FMI_search::get_sa_entries(SMEM *smemArray, int64_t *coordArray, int32_t *coordCountArray, uint32_t count, int32_t max_occ)
{
    uint32_t i;
    int32_t totalCoordCount = 0;
    for (i = 0; i < count; i++)
    {
        int32_t c = 0;
        SMEM smem = smemArray[i];
        int64_t hi = smem.k + smem.s;
        int64_t step = (smem.s > max_occ) ? smem.s / max_occ : 1;
        int64_t j;
        for (j = smem.k; (j < hi) && (c < max_occ); j += step, c++)
        {
            int64_t pos = j;
            int64_t sa_entry = sa_ms_byte[pos];
            sa_entry = sa_entry << 32;
            sa_entry = sa_entry + sa_ls_word[pos];
            //_mm_prefetch((const char *)(sa_ms_byte + pos + SAL_PFD * step), _MM_HINT_T0);
            coordArray[totalCoordCount + c] = sa_entry;
        }
        coordCountArray[i] = c;
        totalCoordCount += c;
    }
}

// sa_compression
// In compression, int64_t reference_seq_len_ = (reference_seq_len >> SA_COMPX(03) ) + 1;   and sa_ls_word and sa_ms_byte were made to be of size reference_seq_len_ using malloc
// pos can go till reference_seq_len
// if ((pos & SA_COMPX_MASK(0x7) ) == 0):  sa_entry = sa_ms_byte[pos >> SA_COMPX]*2^32 + sa_ls_word[pos >> SA_COMPX]  (pos value exceeds 2^SA_COMPX_MASK or referece_seq_len_, which is size of sa_ls_word)
// else: sa_entry = sa_ms_byte[sp >> SA_COMPX]*2^32 + sa_ls_word[sp >> SA_COMPX] + offset

// CT FMI_search::get_sa_entry_compressed(int64_t pos, int tid)
// {
//     if ((pos & SA_COMPX_MASK) == 0) {                                   // SA_COMPX_MASK is 0x7

//         #if  SA_COMPRESSION
//         int64_t sa_entry = sa_ms_byte[pos >> SA_COMPX];
//         CT sa_entry_enc = sa_ms_byte_enc[pos >> SA_COMPX];
//         #else
//         int64_t sa_entry = sa_ms_byte[pos];     // simulation
//         CT sa_entry_enc = sa_ms_byte_enc[pos];
//         #endif

//         sa_entry = sa_entry << 32;
//         sa_entry_enc = shift_left(sa_entry_enc,32);                     // max value 1*2^32

//         #if  SA_COMPRESSION
//         sa_entry = sa_entry + sa_ls_word[pos >> SA_COMPX];
//         sa_entry_enc = cc->EvalAdd(sa_entry_enc,sa_ls_word_enc[pos >> SA_COMPX]);
//         #else
//         sa_entry = sa_entry + sa_ls_word[pos];   // simulation
//         sa_entry_enc = cc->EvalAdd(sa_entry_enc,sa_ls_word_enc[pos]);
//         #endif

//         return sa_entry_enc;
//     }
//     else {
//         // tprof[MEM_CHAIN][tid] ++;
//         int64_t offset = 0;
//         int64_t sp = pos;
//         while(true)
//         {
//             int64_t occ_id_pp_ = sp >> CP_SHIFT;
//             int64_t y_pp_ = CP_BLOCK_SIZE - (sp & CP_MASK) - 1;
//             uint64_t *one_hot_bwt_str = cp_occ[occ_id_pp_].one_hot_bwt_str;
//             uint8_t b;

//             if((one_hot_bwt_str[0] >> y_pp_) & 1)
//                 b = 0;
//             else if((one_hot_bwt_str[1] >> y_pp_) & 1)
//                 b = 1;
//             else if((one_hot_bwt_str[2] >> y_pp_) & 1)
//                 b = 2;
//             else if((one_hot_bwt_str[3] >> y_pp_) & 1)
//                 b = 3;
//             else
//                 b = 4;

//             if (b == 4) {
//                 return encrypt_plaintext_integer_to_ciphertext(offset);
//             }

//             GET_OCC(sp, b, occ_id_sp, y_sp, occ_sp, one_hot_bwt_str_c_sp, match_mask_sp);

//             sp = count[b] + occ_sp;

//             offset ++;
//             // tprof[ALIGN1][tid] ++;
//             if ((sp & SA_COMPX_MASK) == 0) break;
//         }
//         // assert((reference_seq_len >> SA_COMPX) - 1 >= (sp >> SA_COMPX));
//         #if  SA_COMPRESSION
//         int64_t sa_entry = sa_ms_byte[sp >> SA_COMPX];
//         CT sa_entry_enc = sa_ms_byte_enc[sp >> SA_COMPX];
//         #else
//         int64_t sa_entry = sa_ms_byte[sp];      // simultion
//         CT sa_entry_enc = sa_ms_byte_enc[sp];
//         #endif

//         sa_entry = sa_entry << 32;
//         sa_entry_enc = shift_left(sa_entry_enc,32);

//         #if  SA_COMPRESSION
//         sa_entry = sa_entry + sa_ls_word[sp >> SA_COMPX];
//         sa_entry_enc = cc->EvalAdd(sa_entry_enc, sa_ls_word_enc[sp >> SA_COMPX]);
//         #else
//         sa_entry = sa_entry + sa_ls_word[sp];      // simulation
//         sa_entry_enc = cc->EvalAdd(sa_entry_enc, sa_ls_word_enc[sp]);
//         #endif

//         sa_entry += offset;
//         sa_entry_enc=cc->EvalAdd(sa_entry_enc,encode_integer_to_plaintext(offset));

//         return sa_entry_enc;
//     }
// }

// sa_compression
int64_t FMI_search::get_sa_entry_compressed(int64_t pos, int tid)
{
    if ((pos & SA_COMPX_MASK) == 0)
    {

#if SA_COMPRESSION
        int64_t sa_entry = sa_ms_byte[pos >> SA_COMPX];
#else
        int64_t sa_entry = sa_ms_byte[pos];    // simulation
#endif

        sa_entry = sa_entry << 32;

#if SA_COMPRESSION
        sa_entry = sa_entry + sa_ls_word[pos >> SA_COMPX];
#else
        sa_entry = sa_entry + sa_ls_word[pos]; // simulation
#endif

        return sa_entry;
    }
    else
    {
        // tprof[MEM_CHAIN][tid] ++;
        int64_t offset = 0;
        int64_t sp = pos;
        while (true)
        {
            int64_t occ_id_pp_ = sp >> CP_SHIFT;
            int64_t y_pp_ = CP_BLOCK_SIZE - (sp & CP_MASK) - 1;
            uint64_t *one_hot_bwt_str = cp_occ[occ_id_pp_].one_hot_bwt_str;

            int64_t one_hot_bwt_str_0 = decrypt_ciphertext_to_plaintext_vector(cp_occ_one_hot_bwt_str_i(occ_id_pp_, 0))[0];
            int64_t one_hot_bwt_str_1 = decrypt_ciphertext_to_plaintext_vector(cp_occ_one_hot_bwt_str_i(occ_id_pp_, 1))[0];
            int64_t one_hot_bwt_str_2 = decrypt_ciphertext_to_plaintext_vector(cp_occ_one_hot_bwt_str_i(occ_id_pp_, 2))[0];
            int64_t one_hot_bwt_str_3 = decrypt_ciphertext_to_plaintext_vector(cp_occ_one_hot_bwt_str_i(occ_id_pp_, 3))[0];

            uint8_t b;

            // if ((one_hot_bwt_str[0] >> y_pp_) & 1)
            //     b = 0;
            // else if ((one_hot_bwt_str[1] >> y_pp_) & 1)
            //     b = 1;
            // else if ((one_hot_bwt_str[2] >> y_pp_) & 1)
            //     b = 2;
            // else if ((one_hot_bwt_str[3] >> y_pp_) & 1)
            //     b = 3;
            // else
            //     b = 4;

            if ((one_hot_bwt_str_0 >> y_pp_) & 1)
                b = 0;
            else if ((one_hot_bwt_str_1 >> y_pp_) & 1)
                b = 1;
            else if ((one_hot_bwt_str_2 >> y_pp_) & 1)
                b = 2;
            else if ((one_hot_bwt_str_3 >> y_pp_) & 1)
                b = 3;
            else
                b = 4;

            if (b == 4)
            {
                return offset;
            }

            GET_OCC(sp, b, occ_id_sp, y_sp, occ_sp, one_hot_bwt_str_c_sp, match_mask_sp);

            sp = count[b] + occ_sp;

            offset++;
            // tprof[ALIGN1][tid] ++;
            if ((sp & SA_COMPX_MASK) == 0)
                break;
        }
// assert((reference_seq_len >> SA_COMPX) - 1 >= (sp >> SA_COMPX));
#if SA_COMPRESSION
        int64_t sa_entry = sa_ms_byte[sp >> SA_COMPX];
#else
        int64_t sa_entry = sa_ms_byte[sp];     // simultion
#endif

        sa_entry = sa_entry << 32;

#if SA_COMPRESSION
        sa_entry = sa_entry + sa_ls_word[sp >> SA_COMPX];
#else
        sa_entry = sa_entry + sa_ls_word[sp];  // simulation
#endif

        sa_entry += offset;
        return sa_entry;
    }
}

// void FMI_search::get_sa_entries(SMEM *smemArray, int64_t *coordArray, int32_t *coordCountArray, uint32_t count, int32_t max_occ, int tid)
// {

//     uint32_t i;
//     int32_t totalCoordCount = 0;
//     for (i = 0; i < count; i++)
//     {
//         int32_t c = 0;
//         SMEM smem = smemArray[i];
//         int64_t hi = smem.k + smem.s;
//         int64_t step = (smem.s > max_occ) ? smem.s / max_occ : 1;
//         int64_t j;
//         for (j = smem.k; (j < hi) && (c < max_occ); j += step, c++)
//         {
//             int64_t pos = j;
//             int64_t sa_entry = _compressed(pos, tid);
//             coordArray[totalCoordCount + c] = sa_entry;
//         }
//         // coordCountArray[i] = c;
//         *coordCountArray += c;
//         totalCoordCount += c;
//     }
// }

// SA_COPMRESSION w/ PREFETCH
int64_t FMI_search::call_one_step(int64_t pos, int64_t &sa_entry, int64_t &offset)
{
    if ((pos & SA_COMPX_MASK) == 0)
    {
        sa_entry = sa_ms_byte[pos >> SA_COMPX];
        sa_entry = sa_entry << 32;
        sa_entry = sa_entry + sa_ls_word[pos >> SA_COMPX];
        // return sa_entry;
        return 1;
    }
    else
    {
        // int64_t offset = 0;
        int64_t sp = pos;

        int64_t occ_id_pp_ = sp >> CP_SHIFT;
        int64_t y_pp_ = CP_BLOCK_SIZE - (sp & CP_MASK) - 1;
        // vector<int64_t> vv; vv = convert_ciphertext_vector_to_plaintext_vector(cp_occ[occ_id_pp_].one_hot_bwt_str_enc);
        uint64_t *one_hot_bwt_str = cp_occ[occ_id_pp_].one_hot_bwt_str;
        uint8_t b;

        if ((one_hot_bwt_str[0] >> y_pp_) & 1)
            b = 0;
        else if ((one_hot_bwt_str[1] >> y_pp_) & 1)
            b = 1;
        else if ((one_hot_bwt_str[2] >> y_pp_) & 1)
            b = 2;
        else if ((one_hot_bwt_str[3] >> y_pp_) & 1)
            b = 3;
        else
            b = 4;
        if (b == 4)
        {
            sa_entry = 0;
            return 1;
        }

        GET_OCC(sp, b, occ_id_sp, y_sp, occ_sp, one_hot_bwt_str_c_sp, match_mask_sp);

        sp = count[b] + occ_sp;

        offset++;
        if ((sp & SA_COMPX_MASK) == 0)
        {

            sa_entry = sa_ms_byte[sp >> SA_COMPX];
            sa_entry = sa_entry << 32;
            sa_entry = sa_entry + sa_ls_word[sp >> SA_COMPX];

            sa_entry += offset;
            // return sa_entry;
            return 1;
        }
        else
        {
            sa_entry = sp;
            return 0;
        }
    } // else
}

void FMI_search::get_sa_entries_prefetch(SMEM *smemArray, int64_t *coordArray,
                                         int64_t *coordCountArray, int64_t count,
                                         const int32_t max_occ, int tid, int64_t &id_)
{

    // uint32_t i;
    int32_t totalCoordCount = 0;
    int32_t mem_lim = 0, id = 0;

    for (int i = 0; i < count; i++)
    {
        int32_t c = 0;
        SMEM smem = smemArray[i];
        mem_lim += smem.s;
    }

    int64_t *pos_ar = (int64_t *)_mm_malloc(mem_lim * sizeof(int64_t), 64);
    int64_t *map_ar = (int64_t *)_mm_malloc(mem_lim * sizeof(int64_t), 64);

    for (int i = 0; i < count; i++)
    {
        int32_t c = 0;
        SMEM smem = smemArray[i];
        int64_t hi = smem.k + smem.s;
        int64_t step = (smem.s > max_occ) ? smem.s / max_occ : 1;
        int64_t j;
        for (j = smem.k; (j < hi) && (c < max_occ); j += step, c++)
        {
            int64_t pos = j;
            pos_ar[id] = pos;
            map_ar[id++] = totalCoordCount + c;
            // int64_t sa_entry = _compressed(pos, tid);
            // coordArray[totalCoordCount + c] = sa_entry;
        }
        // coordCountArray[i] = c;
        *coordCountArray += c;
        totalCoordCount += c;
    }

    id_ += id;

    const int32_t sa_batch_size = 20;
    int64_t working_set[sa_batch_size], map_pos[sa_batch_size];
    ;
    int64_t offset[sa_batch_size] = {-1};

    int i = 0, j = 0;
    while (i < id && j < sa_batch_size)
    {
        int64_t pos = pos_ar[i];
        working_set[j] = pos;
        map_pos[j] = map_ar[i];
        offset[j] = 0;

        if (pos & SA_COMPX_MASK == 0)
        {
            _mm_prefetch(&sa_ms_byte[pos >> SA_COMPX], _MM_HINT_T0);
            _mm_prefetch(&sa_ls_word[pos >> SA_COMPX], _MM_HINT_T0);
        }
        else
        {
            int64_t occ_id_pp_ = pos >> CP_SHIFT;
            _mm_prefetch(&cp_occ[occ_id_pp_], _MM_HINT_T0);
        }
        i++;
        j++;
    }

    int lim = j, all_quit = 0;
    while (all_quit < id)
    {

        for (int k = 0; k < lim; k++)
        {
            int64_t sp = 0, pos = 0;
            bool quit;
            if (offset[k] >= 0)
            {
                quit = call_one_step(working_set[k], sp, offset[k]);
            }
            else
                continue;

            if (quit)
            {
                coordArray[map_pos[k]] = sp;
                all_quit++;

                if (i < id)
                {
                    pos = pos_ar[i];
                    working_set[k] = pos;
                    map_pos[k] = map_ar[i++];
                    offset[k] = 0;

                    if (pos & SA_COMPX_MASK == 0)
                    {
                        _mm_prefetch(&sa_ms_byte[pos >> SA_COMPX], _MM_HINT_T0);
                        _mm_prefetch(&sa_ls_word[pos >> SA_COMPX], _MM_HINT_T0);
                    }
                    else
                    {
                        int64_t occ_id_pp_ = pos >> CP_SHIFT;
                        _mm_prefetch(&cp_occ[occ_id_pp_], _MM_HINT_T0);
                    }
                }
                else
                    offset[k] = -1;
            }
            else
            {
                working_set[k] = sp;
                if (sp & SA_COMPX_MASK == 0)
                {
                    _mm_prefetch(&sa_ms_byte[sp >> SA_COMPX], _MM_HINT_T0);
                    _mm_prefetch(&sa_ls_word[sp >> SA_COMPX], _MM_HINT_T0);
                }
                else
                {
                    int64_t occ_id_pp_ = sp >> CP_SHIFT;
                    _mm_prefetch(&cp_occ[occ_id_pp_], _MM_HINT_T0);
                }
            }
        }
    }

    _mm_free(pos_ar);
    _mm_free(map_ar);
}
