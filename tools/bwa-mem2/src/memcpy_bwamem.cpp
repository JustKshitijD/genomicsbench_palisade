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

#include "memcpy_bwamem.h"
# include "../../../palisade_header.h"

errno_t memcpy_bwamem(void *dest, rsize_t dmax, const void *src, rsize_t smax, char *file_name, int line_num)
{
    errno_t ret;
    int64_t bytes_copied;
    for(bytes_copied = 0; bytes_copied < smax; bytes_copied += RSIZE_MAX_MEM)
    {
        //cout<<"bytes_copied: "<<bytes_copied<<endl;
        int64_t bytes_remaining = smax - bytes_copied;
        int64_t bytes_to_copy = (bytes_remaining > RSIZE_MAX_MEM) ? RSIZE_MAX_MEM : bytes_remaining;
        int64_t dest_bytes_remaining = dmax - bytes_copied;
        int64_t dest_bytes = (dest_bytes_remaining < bytes_to_copy) ? dest_bytes_remaining : bytes_to_copy;
        //cout<<"bytes_remaining: "<<bytes_remaining<<"; bytes_to_copy: "<<bytes_to_copy<<"; dest_bytes_remaining: "<<dest_bytes_remaining<<"; dest_bytes: "<<dest_bytes<<endl;
        //printf("(const char *)src + bytes_copied: %s\n",(const char *)src + bytes_copied);
        //printf("(char *)dest + bytes_copied: %s\n",(char *)dest + bytes_copied); 
        if((ret = memcpy_s((char *)dest + bytes_copied, dest_bytes, (const char *)src + bytes_copied, bytes_to_copy)) != 0)
        {
            fprintf(stderr, "[%s: %d] memcpy_s returned %d\n", file_name, line_num, ret);
            exit(EXIT_FAILURE);
        }
    }
    //cout<<"Returning from memcpy_bwamem!"<<endl;
    return 0;
}

// memcpy_bwamem_enc(str->enc_s, str->l, str->m - str->l, ks->enc_buf, ks->begin, i - ks->begin, __FILE__, __LINE__);
errno_t memcpy_bwamem_enc(vecCT* dest, int l, rsize_t dmax, vecCT* src, int begin, rsize_t smax, char *file_name, int line_num)
{
    //printf("l: %d; begin: %d\n",l,begin);
    /* printf("In memcpy_bwamem_enc\n");*/
    errno_t ret;
    int64_t bytes_copied;
    for(bytes_copied = 0; bytes_copied < smax; bytes_copied += RSIZE_MAX_MEM)
    {
        //printf("bytes_copied: %ld\n",bytes_copied);

        int64_t bytes_remaining = smax - bytes_copied;
        int64_t bytes_to_copy = (bytes_remaining > RSIZE_MAX_MEM) ? RSIZE_MAX_MEM : bytes_remaining;
        int64_t dest_bytes_remaining = dmax - bytes_copied;
        int64_t dest_bytes = (dest_bytes_remaining < bytes_to_copy) ? dest_bytes_remaining : bytes_to_copy;

        try{
            if((*dest).size()==0 || (*src).size()==0)
                throw(1);
            else if(dest_bytes>RSIZE_MAX_MEM || bytes_to_copy>RSIZE_MAX_MEM)
                throw(2);
            else if(bytes_to_copy>dest_bytes)
                throw(3);
            else
            {
                for(int i=0;i<bytes_to_copy;i++)
                {
                    //printf("i: %d\n",i);
                    (*dest)[l+bytes_copied+i]=(*src)[begin+bytes_copied+i];
                    //printf("decrypt_ciphertext_to_plaintext_vector((*dest)[l+bytes_copied+i])[0]: %c\n",(char)(decrypt_ciphertext_to_plaintext_vector((*dest)[l+bytes_copied+i])[0]));
                }
            }             
        }
        catch(int x){
            if(x==1)
            {
                printf("dest or src is NULL!\n");
            }
            else if(x==2)
            {
                printf("dest_bytes>RSIZE_MAX_MEM || bytes_to_copy>RSIZE_MAX_MEM\n");
            }
            else
            {
                printf("bytes_to_copy>dest_bytes\n");
            }

            if((*dest).size()!=0)
            {
                for(int i=l+bytes_copied;i<l+bytes_copied+dest_bytes;i++)
                    (*dest)[i]=0;
            }
            exit(EXIT_FAILURE);
        }

        // if((ret = memcpy_s((char *)dest + bytes_copied, dest_bytes, (const char *)src + bytes_copied, bytes_to_copy)) != 0)
        // {
        //     fprintf(stderr, "[%s: %d] memcpy_s returned %d\n", file_name, line_num, ret);
        //     exit(EXIT_FAILURE);
        // }
    }
    //printf("Returning 0!\n");
    return 0;
}

// // memcpy_bwamem_enc(str->enc_s, str->l, str->m - str->l, ks->enc_buf, ks->begin, i - ks->begin, __FILE__, __LINE__);
// errno_t memcpy_bwamem_enc(vecCT* dest, int l, rsize_t dmax, vecCT* src, int begin, rsize_t smax, char *file_name, int line_num)
// {
//     printf("l: %d; begin: %d\n",l,begin);
//     /* printf("In memcpy_bwamem_enc\n");*/
//     errno_t ret;
//     int64_t bytes_copied;
//     for(bytes_copied = 0; bytes_copied < smax; bytes_copied += RSIZE_MAX_MEM)
//     {
//         printf("bytes_copied: %ld\n",bytes_copied);

//         int64_t bytes_remaining = smax - bytes_copied;
//         int64_t bytes_to_copy = (bytes_remaining > RSIZE_MAX_MEM) ? RSIZE_MAX_MEM : bytes_remaining;
//         int64_t dest_bytes_remaining = dmax - bytes_copied;
//         int64_t dest_bytes = (dest_bytes_remaining < bytes_to_copy) ? dest_bytes_remaining : bytes_to_copy;

//         try{
//             if(dest.size()==0 || src.size()==0)
//                 throw(1);
//             else if(dest_bytes>RSIZE_MAX_MEM || bytes_to_copy>RSIZE_MAX_MEM)
//                 throw(2);
//             else if(bytes_to_copy>dest_bytes)
//                 throw(3);
//             else
//             {
//                 for(int i=0;i<bytes_to_copy;i++)
//                 {
//                     printf("i: %d\n",i);
//                     dest[l+bytes_copied+i]=src[begin+bytes_copied+i];
//                     printf("decrypt_ciphertext_to_plaintext_vector(dest[l+bytes_copied+i])[0]: %c\n",(char)(decrypt_ciphertext_to_plaintext_vector(dest[l+bytes_copied+i])[0]));
//                 }
//             }             
//         }
//         catch(int x){
//             if(x==1)
//             {
//                 printf("dest or src is NULL!\n");
//             }
//             else if(x==2)
//             {
//                 printf("dest_bytes>RSIZE_MAX_MEM || bytes_to_copy>RSIZE_MAX_MEM\n");
//             }
//             else
//             {
//                 printf("bytes_to_copy>dest_bytes\n");
//             }

//             if(dest.size()!=0)
//             {
//                 for(int i=l+bytes_copied;i<l+bytes_copied+dest_bytes;i++)
//                     dest[i]=0;
//             }
//             exit(EXIT_FAILURE);
//         }

//         // if((ret = memcpy_s((char *)dest + bytes_copied, dest_bytes, (const char *)src + bytes_copied, bytes_to_copy)) != 0)
//         // {
//         //     fprintf(stderr, "[%s: %d] memcpy_s returned %d\n", file_name, line_num, ret);
//         //     exit(EXIT_FAILURE);
//         // }
//     }
//     printf("Returning 0!\n");
//     return 0;
// }


