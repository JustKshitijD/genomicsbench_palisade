// Copyright(c) 2016-2017, Intel Corporation
//
// Redistribution  and  use  in source  and  binary  forms,  with  or  without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of  source code  must retain the  above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
// * Neither the name  of Intel Corporation  nor the names of its contributors
//   may be used to  endorse or promote  products derived  from this  software
//   without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,  BUT NOT LIMITED TO,  THE
// IMPLIED WARRANTIES OF  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED.  IN NO EVENT  SHALL THE COPYRIGHT OWNER  OR CONTRIBUTORS BE
// LIABLE  FOR  ANY  DIRECT,  INDIRECT,  INCIDENTAL,  SPECIAL,  EXEMPLARY,  OR
// CONSEQUENTIAL  DAMAGES  (INCLUDING,  BUT  NOT LIMITED  TO,  PROCUREMENT  OF
// SUBSTITUTE GOODS OR SERVICES;  LOSS OF USE,  DATA, OR PROFITS;  OR BUSINESS
// INTERRUPTION)  HOWEVER CAUSED  AND ON ANY THEORY  OF LIABILITY,  WHETHER IN
// CONTRACT,  STRICT LIABILITY,  OR TORT  (INCLUDING NEGLIGENCE  OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,  EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#include <omp.h>
#include <getopt.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string.h>
#include <xmmintrin.h>
#include <algorithm>
#include <functional>
#include "pairhmm_common.h"

#include "PairHMMUnitTest.h"
#include "shacc_pairhmm.h"
#include <sys/time.h>
#include <sys/resource.h>
#include <time.h>

// #define VTUNE_ANALYSIS 1

#ifdef VTUNE_ANALYSIS
    #include <ittnotify.h>
#endif

using namespace std;
using namespace shacc_pairhmm;

// #define PRINT_OUTPUT
#define MAX_BATCHES 1000000
#define MAX_BATCH_SIZE 50000

// #define ENABLE_SORT 1

long g_total_cells = 0;

static const char *USAGE_MESSAGE =
"  -f, --testfile                       name of test file\n"
"  -l, --loop                           number of loops\n"
"  -t  --threads                        number of threads\n";

namespace opt
{
    const char* testfile;
    int loop = 1;
    int nThreads = 1;
}

static const char* shortopts = "f:l:t:";

enum { OPT_HELP = 1 };

static const struct option longopts[] = {
    { "testfile",       required_argument, NULL, 'f' },
    { "loop",           required_argument, NULL, 'l' },
    { "threads",        required_argument, NULL, 't' },
    { NULL, 0, NULL, 0 }
};

void initPairHMM();
void computelikelihoodsfloat(testcase *testcases, float *expected_result);
void computelikelihoodsboth(testcase *testcases, double *expected_result, int batch_size);

void normalize(string& str, int min_value=0) {
    for (int i = 0; i < str.length(); i++) {
        str[i] = max(min_value, str[i] - 33);
    }
}

// Every batch has batch.num_reads and batch.num_haps
// For batch.num_reads times, strings bases, q, i, d ,c are read into batch.reads[batch.num_reads]
// For batch.num_haps times, string is read into batch.haps[batch.num_haps]
Batch read_batch(istream &is) {
    Batch batch;

    is >> batch.num_reads >> batch.num_haps >> ws;
    // cout<<"-----------------------------------------------------------"<<endl;
    printf("batch.num_reads: %d, batch.num_haps: %d\n",batch.num_reads,batch.num_haps);

    batch.reads = (Read*)_mm_malloc(batch.num_reads * sizeof(Read), 64);
    batch.haps = (Haplotype*)_mm_malloc(batch.num_haps * sizeof(Haplotype), 64);
    batch.results = (double*)_mm_malloc(sizeof(double) * batch.num_reads * batch.num_haps, 64);

    long total_read_length = 0;
    long total_hap_length = 0;

    for (int r = 0; r < batch.num_reads; r++) {
        string bases, q, i, d, c;
        is >> bases >> q >> i >> d >> c >> ws;
        cout<<"bases: "<<bases<<endl;
        cout<<"q: "<<q<<endl;
        cout<<"i: "<<i<<endl;
        cout<<"d: "<<d<<endl;
        cout<<"c: "<<c<<endl;
        
        normalize(q, 6);
        normalize(i);
        normalize(d);
        normalize(c);

        printf("new q: ");
        for(int kk=0;kk<q.length();kk++)
            printf("%d ",q[kk]);
        printf("\n");
        printf("new i: ");
        for(int kk=0;kk<i.length();kk++)
            printf("%d ",i[kk]);
        printf("\n");
        printf("new d: ");
        for(int kk=0;kk<d.length();kk++)
            printf("%d ",d[kk]);
        printf("\n");
        printf("new c: ");
        for(int kk=0;kk<c.length();kk++)
            printf("%d ",c[kk]);
        printf("\n");

        int length = bases.size();
        total_read_length += length;

        Read* read = &batch.reads[r];
        read->length = length;
        read->bases = strndup(bases.c_str(), length);
        read->q = strndup(q.c_str(), length);
        read->i = strndup(i.c_str(), length);
        read->d = strndup(d.c_str(), length);
        read->c = strndup(c.c_str(), length);
    }

    for (int h = 0; h < batch.num_haps; h++) {
        string bases;
        is >> bases >> ws;
        cout<<"h: "<<h<<"; bases: "<<bases<<endl;
        int length = bases.size();
        total_hap_length += length;

        Haplotype* hap = &batch.haps[h];
        hap->length = length;
        hap->bases = strndup(bases.c_str(), length);
    }

    batch.num_cells = total_read_length * total_hap_length;

    return batch;
}

// struct Batch {
//     int id;
//     int num_reads;
//     int num_haps;
//     long num_cells;
//     Read* reads;
//     Haplotype* haps;
//     double* results;

//     bool operator < (const Batch& b) const
//     {
//         return (num_cells < b.num_cells);
//     }
//   };

vector<Batch> read_testfile(string filename) {
    istream *is;
    ifstream ifs;
    if (filename == "") {
        printf("Reading test data from stdin");
        is = &std::cin;
    }
    else {
        printf("Reading test data from file: %s\n", filename.c_str());
        ifs.open(filename.c_str());
        if (!ifs.is_open()) {
            printf("Cannot open file : %s", filename.c_str());
            exit(0);
        }
        is = &ifs;
    }

    vector<Batch> batches;
    fprintf(stderr, "size of 1 batch: %d\n", sizeof(Batch));
    int count = 0;
    while (!is->eof()) {
        Batch batch = read_batch(*is);
        batch.id = count++;
        batches.push_back(batch);
    }

    return batches;
}

int main(int argc, char** argv) {
#ifdef VTUNE_ANALYSIS
    __itt_pause();
#endif
    if (argc == 1) {
        std::cout << USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
        printf("c: %c, optarg: %s\n",c,optarg);
        switch (c) {
            case 'f': opt::testfile = optarg; break;                    // ../../input-datasets/phmm/small/5m.in
            case 'l': opt::loop = stoi(optarg); break;
            case 't': opt::nThreads = stoi(optarg); break;              // 1
            case OPT_HELP:
                std::cout << USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }

    // disable stdout buffering
    setbuf(stdout, NULL);

    initPairHMM();

    ConvertChar::init();

    vector<Batch> batches;
    batches = read_testfile(opt::testfile);
    printf("Num Batches %d, Num threads %d\n", batches.size(), opt::nThreads);

    testcase* testcases = (testcase *)_mm_malloc(MAX_BATCH_SIZE * opt::nThreads * sizeof(testcase), 64);

    #pragma omp parallel num_threads(opt::nThreads)
    {
        int tid = omp_get_thread_num();
        if (tid == 0) {
            printf("Running %d threads\n", opt::nThreads);
        }
    }

    struct timeval start_time, end_time;
    double runtime = 0;

    gettimeofday(&start_time, NULL);

#ifdef VTUNE_ANALYSIS
    __itt_resume();
#endif
    // cout<<"opt::loop: "<<opt::loop<<endl;

    while (opt::loop) {
        #ifdef ENABLE_SORT
            std::sort(batches.begin(), batches.end(), SortByCells());
        #endif

        #pragma omp parallel num_threads(opt::nThreads)
        {
            int tid = omp_get_thread_num();
            // Due to pragma below, each thread will act parallely on a different batch
            #pragma omp for schedule(dynamic)           // The #pragma omp for directive in OpenMP is used to parallelize loops in a shared-memory parallel programming model. It allows the iterations of a loop to be executed in parallel by multiple threads.
                for (int i = 0; i < batches.size(); i++) {
                    printf("----------------\n");
                    printf("i: %d\n",i);
                    int batch_size = batches[i].num_reads * batches[i].num_haps;
                    assert(batch_size <= MAX_BATCH_SIZE);
                    testcase *tc = &testcases[tid * MAX_BATCH_SIZE];
                    // For each batch, we will create testcases out of each of the different combinations of batch.reads[r] and batch.haps[h]
                    // Given batch i, a particular thread tid will work for all its num_reads and num_haps
                    for (int r = 0; r < batches[i].num_reads; r++) {
                        for (int h = 0; h < batches[i].num_haps; h++) {
                            // printf("r: %d, h: %d\n",r,h);
                            tc->rslen = batches[i].reads[r].length;
                            tc->haplen = batches[i].haps[h].length;
                            tc->hap = batches[i].haps[h].bases;
                            tc->rs = batches[i].reads[r].bases;
                            tc->q = batches[i].reads[r].q;
                            tc->i = batches[i].reads[r].i;
                            tc->d = batches[i].reads[r].d;
                            tc->c = batches[i].reads[r].c;
                            tc++;
                        }
                    }
                    // Use a separate testcase for each call
                    computelikelihoodsboth(&testcases[tid * MAX_BATCH_SIZE], batches[i].results, batch_size);
                }
        }
        opt::loop--;

        #ifdef ENABLE_SORT
            std::sort(batches.begin(), batches.end(), SortById());
        #endif
    }
#ifdef VTUNE_ANALYSIS
    __itt_pause();
#endif

    gettimeofday(&end_time, NULL);
    runtime += (end_time.tv_sec - start_time.tv_sec)*1e6 + end_time.tv_usec - start_time.tv_usec;

    for (int i = 0; i < batches.size(); i++) {
#ifdef PRINT_OUTPUT
        int batch_size = batches[i].num_reads * batches[i].num_haps;
        for (int j = 0; j < batch_size; j++) {
            printf("%lf\n", batches[i].results[j]);
        }
#endif
        _mm_free(batches[i].reads);
        _mm_free(batches[i].haps);
        _mm_free(batches[i].results);
    }
    _mm_free(testcases);
    printf("\nPairHMM completed. Kernel runtime: %.2f sec\n", runtime*1e-6);

}
