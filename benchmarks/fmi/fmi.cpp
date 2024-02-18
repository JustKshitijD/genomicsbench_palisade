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

#include<stdio.h>
#include<stdlib.h>
#include<stdint.h>
#include <omp.h>
#include <string.h>
#include "bwa.h"            
#include "FMI_search.h"
#include <chrono>

using namespace std::chrono;

#include <sys/types.h>
#include <unistd.h>

#include <ctime>

// #define VTUNE_ANALYSIS 1

#ifdef VTUNE_ANALYSIS
    #include <ittnotify.h>
#endif

#define CLMUL 8 // 8 x 64 bit cache line
// #define QUERY_DB_SIZE 1280000000
#define QUERY_DB_SIZE 2560000000L
int myrank, num_ranks;

int main(int argc, char **argv) {
    auto fmi_start_1 = high_resolution_clock::now();

    pid_t pid = getpid();
    printf("pid: %lu\n", pid);

    clock_t begin = clock();

#ifdef VTUNE_ANALYSIS
    __itt_pause();
#endif

    if(argc!=6)
    {
        printf("Need five arguments : ref_file query_set batch_size minSeedLen n_threads\n");
        return 1;
    }

    int32_t numReads = 0;
    int64_t total_size = 0;
    gzFile fp = gzopen(argv[2], "r");
	if (fp == 0)
	{
		fprintf(stderr, "[E::%s] fail to open file `%s'.\n", __func__, argv[2]);
        exit(EXIT_FAILURE);
	}
    
    // ../benchmarks/fmi/fmi $INPUTS_DIR/fmi/broad $INPUTS_DIR/fmi/small/SRR7733443_1m_1.fastq 512 19 1

    auto fmi_end_1 = high_resolution_clock::now();

    auto bseq_start = high_resolution_clock::now();
    // printf("before reading sequences\n");
    bseq1_t *seqs = bseq_read_one_fasta_file(QUERY_DB_SIZE, &numReads, fp, &total_size); // seqs array now has all input  
    auto bseq_end = high_resolution_clock::now();

    auto duration = duration_cast<microseconds>(bseq_end - bseq_start);

    cout << "Time taken by bseq_read_one_fasta_file: "
         << duration.count() << " microseconds" << endl;

    if(seqs == NULL)
    {
        printf("ERROR! seqs = NULL\n");
        exit(EXIT_FAILURE);
    }
    int32_t *query_cum_len_ar = (int32_t *)_mm_malloc(numReads * sizeof(int32_t), 64);

    // vecCT ref_file_name_enc;
    // assign_string_to_vecCT(ref_file_name_enc,argv[1],-1);
    FMI_search *fmiSearch = new FMI_search(argv[1]); // argv[1] is $INPUTS_DIR/fmi/broad
    fmiSearch->load_index();
    // printf("Out of load index!\n");

    auto fmi_start_2 = high_resolution_clock::now();

    int max_readlength = seqs[0].l_seq;
    int min_readlength = seqs[0].l_seq;
    for(int i = 1; i < numReads; i++)
    {
        if(max_readlength < seqs[i].l_seq)
            max_readlength = seqs[i].l_seq;
        if(min_readlength > seqs[i].l_seq)
            min_readlength = seqs[i].l_seq;
    }
    assert(max_readlength > 0);
    assert(max_readlength < 10000);
    assert(numReads > 0);
    assert(numReads * max_readlength < QUERY_DB_SIZE); // #define QUERY_DB_SIZE 2560000000L
    printf("numReads = %d, max_readlength = %d, min_readlength = %d\n", numReads, max_readlength, min_readlength);
    uint8_t *enc_qdb=(uint8_t *)malloc((int64_t)numReads * max_readlength * sizeof(uint8_t));

    int64_t cind,st;
#if 0
    printf("Priting query\n");
    for(st = 0; st < max_readlength; st++)
    {
        printf("%c", seqs[0].seq[st]);
    }
    printf("\n");
#endif

    // int32_t *query_cum_len_ar = (int32_t *)_mm_malloc(numReads * sizeof(int32_t), 64);
    uint64_t r;
    for (st=0; st < numReads; st++) {
        // printf("st: %d, max_readlength: %d\n",st,max_readlength);
        query_cum_len_ar[st] = st * max_readlength;
        cind=st*max_readlength;
        for(r = 0; r < max_readlength; ++r) {
            // printf("decrypt_ciphertext_to_plaintext_vector(seqs[st].enc_seq[%d])[0]: %d\n",r,decrypt_ciphertext_to_plaintext_vector(seqs[st].enc_seq[r])[0]);
            // printf("r+cind: %d\n",r+cind);
            // printf("decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(seqs[st].enc_seq[r], encode_integer_to_plaintext('C')))[0]: %d\n",decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(seqs[st].enc_seq[r], encode_integer_to_plaintext('C')))[0]);

            if(decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(seqs[st].enc_seq[r], encode_integer_to_plaintext('A')))[0]==0)
                enc_qdb[r+cind]=0;
            else if(decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(seqs[st].enc_seq[r], encode_integer_to_plaintext('C')))[0]==0)
                enc_qdb[r+cind]=1;
            else if(decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(seqs[st].enc_seq[r], encode_integer_to_plaintext('G')))[0]==0)
                enc_qdb[r+cind]=2;
            else if(decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(seqs[st].enc_seq[r], encode_integer_to_plaintext('T')))[0]==0)
                enc_qdb[r+cind]=3;     
            else 
                enc_qdb[r+cind]=4;

            // switch(seqs[st].seq[r])
            // {
            //     case 'A': enc_qdb[r+cind]=0;
            //               break;
            //     case 'C': enc_qdb[r+cind]=1;
            //               break;
            //     case 'G': enc_qdb[r+cind]=2;
            //               break;
            //     case 'T': enc_qdb[r+cind]=3;
            //               break;
            //     default: enc_qdb[r+cind]=4;
            // }

            // printf("%c %d\n", decrypt_ciphertext_to_plaintext_vector(seqs[st].enc_seq[r])[0], enc_qdb[r + cind]);
            // printf("enc_qdb[%d]= %d\n",r+cind, enc_qdb[r + cind]);
        }
    }

    int batch_size=0;
    batch_size=atoi(argv[3]);
    assert(batch_size > 0);
    assert(batch_size <= numReads);

    int32_t minSeedLen = atoi(argv[4]);         // 19- length of the substring of a 'read' from the fastq file, whose position is found in reference sequence using FMI kernel
    int numthreads=atoi(argv[5]);               // 1

    const int splitWidth = 10;
    const int maxMemIntv = 20;
    const double splitFactor = 1.5;

//     typedef struct smem_struct
// {
// #ifdef DEBUG
//     uint64_t info; // for debug
// #endif
//     uint32_t rid;
//     uint32_t m, n;
//     int64_t k, l, s;
// }SMEM;

// #define CLMUL 8 // 8 x 64 bit cache line

    assert(numthreads > 0);
    assert(numthreads <= omp_get_max_threads());                                 //   # define CLMUL 8
    SMEM *matchArray[CLMUL * numthreads];
    int32_t *min_intv_array[CLMUL * numthreads];
    int32_t *rid_array[CLMUL * numthreads];
    int16_t *query_pos_array[CLMUL * numthreads];

    int64_t num_batches = (numReads + batch_size - 1 ) / batch_size;
    int64_t *numTotalSmem = (int64_t *)_mm_malloc(num_batches * sizeof(int64_t), 64);;
    SMEM **batchStart = (SMEM **)_mm_malloc(num_batches * sizeof(SMEM *), 64);;
    
    uint64_t tim = __rdtsc();
    sleep(1);
    uint64_t proc_freq = __rdtsc() - tim;

#pragma omp parallel num_threads(numthreads)
    {
        int tid = omp_get_thread_num();

        if(tid == 0)
            printf("Running %d threads\n", omp_get_num_threads());
    }

    int64_t i;

    int64_t startTick, endTick;
#ifdef VTUNE_ANALYSIS
    __itt_resume();
#endif
    startTick = __rdtsc();
    memset(numTotalSmem, 0, num_batches * sizeof(int64_t));
    memset(batchStart, 0, num_batches * sizeof(int64_t));
    int64_t workTicks[CLMUL * numthreads];
    memset(workTicks, 0, CLMUL * numthreads * sizeof(int64_t));
    int64_t perThreadQuota = numReads / numthreads;

    int split_len = (int)(minSeedLen * splitFactor + .499);

#pragma omp parallel num_threads(numthreads)
    {
        int32_t tid = omp_get_thread_num();
        int64_t matchArrayAlloc = perThreadQuota * 20;
        matchArray[CLMUL * tid] = (SMEM *)malloc(matchArrayAlloc * sizeof(SMEM));
        min_intv_array[CLMUL * tid] = (int32_t *)malloc(matchArrayAlloc * sizeof(int32_t));
        rid_array[CLMUL * tid] = (int32_t *)malloc(matchArrayAlloc * sizeof(int32_t));
        query_pos_array[CLMUL * tid] = (int16_t *)malloc(matchArrayAlloc * sizeof(int16_t));

        int64_t myTotalSmems = 0;
        int64_t startTick = __rdtsc();

#pragma omp for schedule(dynamic, 1)
        for(i = 0; i < numReads; i += batch_size)
        {
            // printf("i: %d\n",i);
            int64_t st1 = __rdtsc();
            int32_t batch_count = batch_size;
            if((i + batch_count) > numReads) batch_count = numReads - i;
            int32_t j;
            for(j = 0; j < batch_count; j++)
            {
                min_intv_array[CLMUL * tid][j] = 1;
                rid_array[CLMUL * tid][j] = j;
            }
            int32_t batch_id = i/batch_size;
            printf("%d] i = %d, batch_count = %d, batch_size = %d, batch_id= %d\n", tid, i, batch_count, batch_size,batch_id);
            //fflush(stdout);
            if((matchArrayAlloc - myTotalSmems) < (batch_size * max_readlength))
            {
                printf("%d] realloc\n", tid);
                fflush(stdout);
                matchArrayAlloc *= 2;
                matchArray[CLMUL * tid] = (SMEM *)realloc(matchArray[CLMUL * tid], matchArrayAlloc * sizeof(SMEM)); 
                min_intv_array[CLMUL * tid] = (int32_t *)realloc(min_intv_array[CLMUL * tid], matchArrayAlloc * sizeof(int32_t)); 
                rid_array[CLMUL * tid] = (int32_t *)realloc(rid_array[CLMUL * tid], matchArrayAlloc * sizeof(int32_t));
                query_pos_array[CLMUL * tid] = (int16_t *)realloc(query_pos_array[CLMUL * tid], matchArrayAlloc * sizeof(int16_t));
            }
            int64_t num_smem1 = 0, num_smem2 = 0, num_smem3 = 0;
            fmiSearch->getSMEMsAllPosOneThread(enc_qdb + i * max_readlength,
                    min_intv_array[CLMUL * tid],
                    rid_array[CLMUL * tid],
                    batch_count,
                    batch_size,
                    seqs + i,
                    query_cum_len_ar,
                    max_readlength,
                    minSeedLen,
                    &matchArray[CLMUL * tid][myTotalSmems],
                    &num_smem1);
            // printf("num_smem1: %d\n",num_smem1);
            int64_t pos = 0;
            for (j = 0; j < num_smem1; j++) {
                SMEM *p = &matchArray[CLMUL * tid][myTotalSmems + j];
                int start = p->m, end = p->n +1;
                if (end - start < split_len || p->s > splitWidth) continue;
                // printf("%d:%d:%ld\n", start, end, p->s);
                rid_array[CLMUL * tid][pos] = p->rid;
                query_pos_array[CLMUL * tid][pos] = (end + start)>>1;
                min_intv_array[CLMUL * tid][pos] = p->s + 1;
                pos++;
            }
            
            // Reseed
            fmiSearch->getSMEMsOnePosOneThread(enc_qdb + i * max_readlength,
                    query_pos_array[CLMUL * tid],
                    min_intv_array[CLMUL * tid],
                    rid_array[CLMUL * tid],
                    pos,
                    pos,
                    seqs + i,
                    query_cum_len_ar,
                    max_readlength,
                    minSeedLen,
                    &matchArray[CLMUL * tid][myTotalSmems + num_smem1],
                    &num_smem2);

            // printf("num_smem2: %d\n",num_smem2);        
            // LAST
            for(j = 0; j < batch_count; j++)
            {
                min_intv_array[CLMUL * tid][j] = maxMemIntv;
            }
            num_smem3 = fmiSearch->bwtSeedStrategyAllPosOneThread(enc_qdb + i * max_readlength,
                    min_intv_array[CLMUL * tid],
                    batch_count,
                    seqs + i,
                    query_cum_len_ar,
                    minSeedLen + 1,
                    &matchArray[CLMUL * tid][myTotalSmems + num_smem1 + num_smem2]);
            int64_t totalSmem = num_smem1 + num_smem2 + num_smem3; 
            printf("num_smem1: %d, num_smem2: %d, num_smem3: %d\n",num_smem1,num_smem2,num_smem3);
            numTotalSmem[batch_id] = totalSmem;
            printf("numTotalSmem[batch_id]: %d\n",numTotalSmem[batch_id]);
            batchStart[batch_id] = matchArray[CLMUL * tid] + myTotalSmems;
            for(j = 0; j < totalSmem; j++)
            {
                matchArray[CLMUL * tid][myTotalSmems + j].rid += i;
            }
            fmiSearch->sortSMEMs(matchArray[CLMUL * tid] + myTotalSmems,
                    numTotalSmem + batch_id,
                    batch_count,
                    max_readlength,
                    1);
            myTotalSmems += totalSmem; 
            int64_t et1 = __rdtsc();
            workTicks[CLMUL * tid] += (et1 - st1);
        }

        int64_t endTick = __rdtsc();
        printf("%d] %ld ticks, workTicks = %ld\n", tid, endTick - startTick, workTicks[CLMUL * tid]);
    }

    auto fmi_end_2 = high_resolution_clock::now();
    auto duration_fmi = duration_cast<microseconds>(fmi_end_1-fmi_start_1+fmi_end_2-fmi_start_2);
    cout << "Time taken by fmi: "<<duration_fmi.count()<< " microseconds" << endl;

    endTick = __rdtsc();
#ifdef VTUNE_ANALYSIS
    __itt_pause();
#endif
    int64_t sumTicks = 0;
    int64_t maxTicks = 0;
    for(i = 0; i < numthreads; i++)
    {
        sumTicks += workTicks[CLMUL * i];
        if(workTicks[CLMUL * i] > maxTicks) maxTicks = workTicks[CLMUL * i];
    }
    double avgTicks = (sumTicks * 1.0) / numthreads;
    printf("avgTicks = %lf, maxTicks = %ld, load imbalance = %lf\n", avgTicks, maxTicks, maxTicks/avgTicks);

    printf("Consumed: %ld cycles, %0.4lf sec\n", endTick - startTick, (endTick - startTick) * 1.0 / proc_freq);

    int64_t totalSmem = 0;
    int32_t batch_id = 0;
    for(batch_id = 0; batch_id < num_batches; batch_id++)
    {
        printf("batch_id: %d, numTotalSmem[batch_id]: %d\n",batch_id,numTotalSmem[batch_id]);
        totalSmem += numTotalSmem[batch_id];
    }
    printf("totalSmems = %ld\n", totalSmem);

#ifdef PRINT_OUTPUT
    int32_t prevRid = -1;
    for(batch_id = 0; batch_id < num_batches; batch_id++)
    {
        SMEM *myMatchArray = batchStart[batch_id];
        int64_t i;
        for(i = 0; i < numTotalSmem[batch_id]; i++)
        {
            SMEM smem = myMatchArray[i];
            if(smem.rid != prevRid)
            {
                int32_t j;
                for(j = prevRid + 1; j <= smem.rid; j++)
                    printf("%u:\n", j);
            }
            prevRid = smem.rid;
            printf("[%u,%u]", smem.m, smem.n + 1);
            // printf("%u, %u]", smem.k, smem.s);
#if 0
            printf(" ["); 
            int64_t u1, u2, u3;
            u1 = smem.k;
            u2 = smem.k + smem.s;
            for(u3 = u1; u3 < u2; u3++)
            {
                printf("%ld,", fmiSearch->get_sa_entry(u3));
            }
            printf("]"); 
#endif
            printf("\n");
        }
    }
#endif
    _mm_free(query_cum_len_ar);
    free(enc_qdb);
    for(int tid = 0; tid < numthreads; tid++)
    {
        free(matchArray[CLMUL * tid]);
        free(min_intv_array[CLMUL * tid]);
        free(query_pos_array[CLMUL * tid]);
        free(rid_array[CLMUL * tid]);
    }
    _mm_free(numTotalSmem);
    _mm_free(batchStart);
    delete fmiSearch;

    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    printf("elapsed secs: %lf\n",elapsed_secs);

    return 0;
}

