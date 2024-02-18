// phmm:- max_base_size: 250; max_q_size: 250; max_i_size: 250; max_d_size: 250; max_c_size: 250; max_num_reads: 110
// count goes from [0,549]; and for each count, a batch is read

// struct Batch {
//   int id;
//   int num_reads;					
//   int num_haps;
//   long num_cells;
//   Read* reads;             // array of size batch.num_reads(=110 for small data; 1193 for large data)
//   Haplotype* haps;				 // array of size batch.num_haps(=37 for small data; 128 for large data)
//   double* results;

//   bool operator < (const Batch& b) const
//   {
//       return (num_cells < b.num_cells);
//   }
// };

// struct Read {
//   int length;
//   const char* bases;				// each is of max length 250
//   const char* q;
//   const char* i;
//   const char* d;
//   const char* c;
// };

// struct Haplotype {
//   int length;
//   const char* bases;		// max haps_base_size is 473 for large and 302 for small
// };

// 1 read has (250*5+1)=1251 ints; step(16384/1251)=13 reads in 1 ciphertext; ceil(batch.num_reads/13) CTs required to cover all reads;
// 1 haplotype has (473+1)=474 ints; step(16384/474)=35 haplotypes in 1 ciphertext; ceil(batch.num_haps/35) CTs required to cover all haplotypes;


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
// #include "../../palisade_header.h"

double runtime = 0;
struct timeval start_time, end_time;

// #define VTUNE_ANALYSIS 1

#ifdef VTUNE_ANALYSIS
#include <ittnotify.h>
#endif

using namespace std;
using namespace shacc_pairhmm;

// #define PRINT_OUTPUT
#define MAX_BATCHES 1000000
#define MAX_BATCH_SIZE 50000

testcase *testcases;
int max_base_size = -100, max_q_size = -100, max_i_size = -100, max_d_size = -100, max_c_size = -100, max_num_reads = -100, max_num_haps = -100, max_haps_bases_size = -100;

// #define ENABLE_SORT 1

long g_total_cells = 0;

static const char *USAGE_MESSAGE =
    "  -f, --testfile                       name of test file\n"
    "  -l, --loop                           number of loops\n"
    "  -t  --threads                        number of threads\n";

namespace opt
{
    const char *testfile;
    int loop = 1;
    int nThreads = 1;
}

static const char *shortopts = "f:l:t:";

enum
{
    OPT_HELP = 1
};

static const struct option longopts[] = {
    {"testfile", required_argument, NULL, 'f'},
    {"loop", required_argument, NULL, 'l'},
    {"threads", required_argument, NULL, 't'},
    {NULL, 0, NULL, 0}};

void initPairHMM();
void computelikelihoodsfloat(testcase *testcases, float *expected_result);
void computelikelihoodsboth(testcase *testcases, double *expected_result, int batch_size);

void normalize(string &str, int min_value = 0)
{
    for (int i = 0; i < str.length(); i++)
    {
        str[i] = max(min_value, str[i] - 33);
    }
}

// Every batch has batch.num_reads and batch.num_haps
// For batch.num_reads times, strings bases, q, i, d ,c are read into batch.reads[batch.num_reads]
// For batch.num_haps times, string is read into batch.haps[batch.num_haps]
void read_batch(istream &is, int count, int nThreads, int loop)
{

    // #pragma omp parallel num_threads(nThreads)
    // {
    //     int tid = omp_get_thread_num();
    //     if (tid == 0) {
    //         printf("Running %d threads\n", nThreads);
    //     }
    // }
    // printf("Outside omp parallel\n");

    Batch batch;
    batch.id = count;

    vecInt int_ct_vector(16384, 0);
    vecInt reads_ct_vector(16384, 0);
    vecInt haps_ct_vector(16384, 0);

    is >> batch.num_reads >> batch.num_haps >> ws;
    CT num_reads_ct, num_haps_ct;
    num_reads_ct = encrypt_plaintext_integer_to_ciphertext(batch.num_reads);
    num_haps_ct = encrypt_plaintext_integer_to_ciphertext(batch.num_haps);

    // cout<<"-----------------------------------------------------------"<<endl;
    printf("batch.num_reads: %d, batch.num_haps: %d\n", batch.num_reads, batch.num_haps);
    max_num_reads = max(max_num_reads, batch.num_reads);
    max_num_haps = max(max_num_haps, batch.num_haps);

    // batch.reads array is of size batch.num_reads
    batch.reads = (Read *)_mm_malloc(batch.num_reads * sizeof(Read), 64);
    batch.haps = (Haplotype *)_mm_malloc(batch.num_haps * sizeof(Haplotype), 64);
    batch.results = (double *)_mm_malloc(sizeof(double) * batch.num_reads * batch.num_haps, 64);
    // batch.reads_data_in_ct.resize(ceil(batch.num_reads*1.0/13));
    // batch.haps_data_in_ct.resize(ceil(batch.num_haps*1.0/35));

    long total_read_length = 0;
    long total_hap_length = 0;

    for (int r = 0; r < batch.num_reads; r++)           // reading into &batches.reads[r]
    {
        printf("r: %d\n", r);
        string bases, q, i, d, c;
        is >> bases >> q >> i >> d >> c >> ws;
        // cout<<"bases: "<<bases<<endl;
        // cout<<"q: "<<q<<endl;
        // cout<<"i: "<<i<<endl;
        // cout<<"d: "<<d<<endl;
        // cout<<"c: "<<c<<endl;

        // max_base_size = max(max_base_size, (int)(bases.size()));
        // max_q_size = max(max_q_size, (int)(q.size()));
        // max_i_size = max(max_i_size, (int)(i.size()));
        // max_d_size = max(max_d_size, (int)(d.size()));
        // max_c_size = max(max_c_size, (int)(c.size()));

        normalize(q, 6);
        normalize(i);
        normalize(d);
        normalize(c);

        // printf("new q: ");
        // for (int kk = 0; kk < q.length(); kk++)
        //     printf("%d ", q[kk]);
        // printf("\n");
        // printf("new i: ");
        // for (int kk = 0; kk < i.length(); kk++)
        //     printf("%d ", i[kk]);
        // printf("\n");
        // printf("new d: ");
        // for (int kk = 0; kk < d.length(); kk++)
        //     printf("%d ", d[kk]);
        // printf("\n");
        // printf("new c: ");
        // for (int kk = 0; kk < c.length(); kk++)
        //     printf("%d ", c[kk]);
        // printf("\n");
        // printf("new bases: ");
        // for (int kk = 0; kk < bases.length(); kk++)
        //     printf("%d ", bases[kk]);
        // printf("\n\n");

        int length = bases.size();
        total_read_length += length;

        Read *read = &batch.reads[r];
        read->length = length;
        read->bases = strndup(bases.c_str(), length);
        read->q = strndup(q.c_str(), length);
        read->i = strndup(i.c_str(), length);
        read->d = strndup(d.c_str(), length);
        read->c = strndup(c.c_str(), length);
    }

    printf("batch.num_reads: %d\n", batch.num_reads);

    CT c1;
    vecInt v3;

    for (int r = 0; r < batch.num_reads; r++)
    {
        printf("------------------------------------------\n");
        printf("r: %d\n",r);
        if (r % 13 == 0 && r != 0)          // 13 reads in 1 CT
        {
            printf("r%13==0 for r: %d\n", r);
            v3.resize(16384, 0);
            c1 = encrypt_plaintext_vector_to_ciphertext(reads_ct_vector);       // size of reads_ct_vector is 16384

            cout<<"reads_ct_vector: "<<endl;
            for (int ii = 0; ii < v3.size(); ii++)
            {
                cout <<reads_ct_vector[ii]<<" ";
            }
            cout<<endl;

            printf("encrypted!\n");
            v3 = decrypt_ciphertext_to_plaintext_vector(c1);
            printf("decrypted!\n");

            cout<<"v3: "<<endl;
            for (int ii = 0; ii < v3.size(); ii++)
            {
                cout <<v3[ii]<<" ";
            }
            cout<<endl;

            batch.reads_data_in_ct.push_back(c1);
            printf("batch.reads_data_in_ct.size(): %d\n", batch.reads_data_in_ct.size());
            for(int kk=0;kk<16384;kk++)
                reads_ct_vector[kk]=0;
            // reads_ct_vector.resize(16384, 0);
            // reads_ct_vector=[0]*16384;
            printf("Zeroed read_ct_vector\n");
            for(int kk=0;kk<reads_ct_vector.size();kk++){
                cout<<reads_ct_vector[kk]<<" ";
            }
            cout<<endl;
            
        }

        // char q[strlen(batch.reads[r].q)+1];
        // char i[strlen(batch.reads[r].i)+1];
        // char d[strlen(batch.reads[r].d)+1];
        // char c[strlen(batch.reads[r].c)+1];
        // char bases[strlen(batch.reads[r].bases)];

        // strcpy(q,batch.reads[r].q); // copy batch.reads[r].q to q
        // strcpy(i,batch.reads[r].i);
        // strcpy(d,batch.reads[r].d);
        // strcpy(c,batch.reads[r].c);
        // strcpy(bases,batch.reads[r].bases);

        // const char *q = batch.reads[r].q;
        // const char *i = batch.reads[r].i;
        // const char *d = batch.reads[r].d;
        // const char *c = batch.reads[r].c;
        // const char *bases = batch.reads[r].bases;

        // int length = batch.reads[r].length;

        // if((r%13)==8){
        //     cout<<"strlen(q): "<<strlen(q)<<endl;
        //     int k=0;
        //     cout<<"q: "<<endl;
        //     while(k<strlen(q))
        //     {
        //         printf("%d ",q[k]);
        //         // cout<<"q["<<k<<"]="<<q[k]<<endl;
        //         k++;
        //     }
        //     cout<<endl;
        // }

        int st = 1251 * (r%13);
        int k = 0;
        while (k < 250 && k < strlen(batch.reads[r].q) && k<batch.reads[r].length )
        {
            reads_ct_vector[st+k] = (int64_t)(batch.reads[r].q[k]);
            k++;
        }
        // reads_ct_vector[st+k]=(int64_t)('\0');
        printf("q: \n");
        for (int kk = 0; kk < strlen(batch.reads[r].q) && kk<batch.reads[r].length ; kk++)
            printf("%d ", kk<batch.reads[r].q[kk]);
        printf("\n");

        printf("After q\n");
        for(int kk=0;kk<reads_ct_vector.size();kk++){
            cout<<reads_ct_vector[kk]<<" ";
        }
        cout<<endl<<endl;

        st = 1251 * (r%13) + 250;
        k = 0;
        while (k <250 && k < strlen(batch.reads[r].i) && k<batch.reads[r].length )
        {
            reads_ct_vector[st+k] = (int64_t)(batch.reads[r].i[k]);
            k++;
        }

        printf("i: \n");
        for (int kk = 0; kk < strlen(batch.reads[r].i) && kk<batch.reads[r].length ; kk++)
            printf("%d ", batch.reads[r].i[kk]);
        printf("\n");
        // reads_ct_vector[st+k]=(int64_t)('\0');
        printf("After i\n");
        for(int kk=0;kk<reads_ct_vector.size();kk++){
            cout<<reads_ct_vector[kk]<<" ";
        }
        cout<<endl<<endl;

        st = 1251 * (r%13) + 500;
        k = 0;
        while (k < 250 && k < strlen(batch.reads[r].d) && k<batch.reads[r].length)
        {
            reads_ct_vector[st+k] = (int64_t)(batch.reads[r].d[k]);
            k++;
        }
        // reads_ct_vector[st+k]=(int64_t)('\0');
        printf("d: \n");
        for (int kk = 0; kk < strlen(batch.reads[r].d) && kk<batch.reads[r].length; kk++)
            printf("%d ", batch.reads[r].d[kk]);
        printf("\n");
        printf("After d\n");
        for(int kk=0;kk<reads_ct_vector.size();kk++){
            cout<<reads_ct_vector[kk]<<" ";
        }
        cout<<endl<<endl;

        st = 1251 * (r%13) + 750;
        k = 0;
        while (k < 250 && k < strlen(batch.reads[r].c) && k<batch.reads[r].length)
        {
            reads_ct_vector[st+k] = (int64_t)(batch.reads[r].c[k]);
            k++;
        }
        // reads_ct_vector[st+k]=(int64_t)('\0');
        printf("c: \n");
        for (int kk = 0; kk < strlen(batch.reads[r].c) && kk<batch.reads[r].length; kk++)
            printf("%d ", batch.reads[r].c[kk]);
        printf("\n");
        printf("After c\n");
        for(int kk=0;kk<reads_ct_vector.size();kk++){
            cout<<reads_ct_vector[kk]<<" ";
        }
        cout<<endl<<endl;

        st = 1251 * (r%13) + 1000;
        k = 0;
        while (k < 250 && k < strlen(batch.reads[r].bases) && k<batch.reads[r].length)
        {
            reads_ct_vector[st+k] = (int64_t)(batch.reads[r].bases[k]);
            k++;
        }

        printf("bases: \n");
        for (int kk = 0; kk < strlen(batch.reads[r].bases) && kk<batch.reads[r].length; kk++)
            printf("%d ", batch.reads[r].bases[kk]);
        printf("\n");
        // reads_ct_vector[st+k]=(int64_t)('\0');
        printf("After bases\n");
        for(int kk=0;kk<reads_ct_vector.size();kk++){
            cout<<reads_ct_vector[kk]<<" ";
        }
        cout<<endl;

        reads_ct_vector[1251 * (r%13) + 1250] = (int64_t)(batch.reads[r].length);
    }
    //reads_data_in_ct is vecCT whose one element is CT that encrypts 13 reads
    batch.reads_data_in_ct.push_back(encrypt_plaintext_vector_to_ciphertext(reads_ct_vector));

    printf("batch.reads_data_in_ct.size(): %d\n", batch.reads_data_in_ct.size());
    for (int r = 0; r < batch.num_reads; r++)
    {
        printf("--------------------------\n");
        printf("r: %d\n", r);

        const char *q = batch.reads[r].q;
        const char *i = batch.reads[r].i;
        const char *d = batch.reads[r].d;
        const char *c = batch.reads[r].c;
        const char *bases = batch.reads[r].bases;
        int length = batch.reads[r].length;
        printf("length: %d\n", length);
        int ct_count = r/13;
        printf("ct_count: %d\n", ct_count);
        CT c2 = batch.reads_data_in_ct[ct_count];
        vecInt v2 = decrypt_ciphertext_to_plaintext_vector(c2);

        printf("q: \n");
        for (int s = 0; s < strlen(q); s++)
            printf("%d ", q[s]);
        printf("\n");
        for (int s = 0; s < 250 && v2[1251*(r%13)+s] != 0; s++)
        {
            printf("%lld ", v2[1251*(r%13)+s]);
        }
        printf("\n\n");

        printf("i: \n");
        for (int s = 0; s < strlen(i); s++)
            printf("%d ", i[s]);
        printf("\n");
        for (int s = 250; s < 500 && v2[1251*(r%13)+s] != 0; s++)
        {
            printf("%lld ", v2[1251*(r%13)+s]);
        }
        printf("\n\n");

        printf("d: \n");
        for (int s = 0; s < strlen(d); s++)
            printf("%d ", d[s]);
        printf("\n");
        for (int s = 500; s < 750 && v2[1251*(r%13)+s] != 0; s++)
        {
            printf("%lld ", v2[1251*(r%13)+s]);
        }
        printf("\n\n");

        printf("c: \n");
        for (int s = 0; s < strlen(c); s++)
            printf("%d ", c[s]);
        printf("\n");
        for (int s = 750; s < 1000 && v2[1251*(r%13)+s] != 0; s++)
        {
            printf("%lld ", v2[1251*(r%13)+s]);
        }
        printf("\n\n");

        printf("bases: \n");
        for (int s = 0; s < strlen(bases); s++)
            printf("%d ", bases[s]);
        printf("\n");
        for (int s = 1000; s < 1250 && v2[1251*(r%13)+s] != 0; s++)
        {
            printf("%lld ", v2[1251*(r%13)+s]);
        }
        printf("\n\n");

        printf("length: \n");
        printf("%d\n", length);
        printf("%lld\n", v2[1251*(r%13)+1250]);
    }

    // printf("batch.reads_data_in_ct.size(): %d\n",batch.reads_data_in_ct.size());
    // for(int kk=0;kk<batch.reads_data_in_ct.size();kk++){
    //     vecInt vv=batch.reads_data_in_ct[kk];
    // }

    // max strlen(batch.haps[r]) across all batches is 473; In this max length case, bases[0] to bases[472] will be copied 
    
    for (int h = 0; h < batch.num_haps; h++)
    {
        string bases;
        is >> bases >> ws;
        cout << "h: " << h << "; bases: " << bases << endl;
        max_haps_bases_size = max(max_haps_bases_size, (int)(bases.size()));
        int length = bases.size();
        total_hap_length += length;

        Haplotype *hap = &batch.haps[h];
        hap->length = length;
        hap->bases = strndup(bases.c_str(), length);    
        // c_str() converts string to array of chars and terminates with \0
    }
    
    for (int r = 0; r < batch.num_haps; r++)
    {
        // const char *bases = batch.haps[r].bases;
        // int length = batch.haps[r].length;

        if(r%35==0 && r!=0){
            printf("r%35==0 for r: %d\n", r);
            v3.resize(16384, 0);
            c1 = encrypt_plaintext_vector_to_ciphertext(haps_ct_vector);       // size of reads_ct_vector is 16384

            cout<<"haps_ct_vector: "<<endl;
            for (int ii = 0; ii < v3.size(); ii++)
            {
                cout <<haps_ct_vector[ii]<<" ";
            }
            cout<<endl;

            printf("encrypted!\n");
            v3 = decrypt_ciphertext_to_plaintext_vector(c1);
            printf("decrypted!\n");

            cout<<"v3: "<<endl;
            for (int ii = 0; ii < v3.size(); ii++)
            {
                cout <<v3[ii]<<" ";
            }
            cout<<endl;

            batch.haps_data_in_ct.push_back(c1);
            printf("batch.haps_data_in_ct.size(): %d\n", batch.haps_data_in_ct.size());
            for(int kk=0;kk<16384;kk++)
                haps_ct_vector[kk]=0;
            // reads_ct_vector.resize(16384, 0);
            // reads_ct_vector=[0]*16384;
            printf("Zeroed haps_ct_vector\n");
            for(int kk=0;kk<haps_ct_vector.size();kk++){
                cout<<haps_ct_vector[kk]<<" ";
            }
            cout<<endl;
        }

        int st = 474 * (r%35);  // 473 for bases; 1 for length; 35 haps in 1 CT
        int k = 0;

        while (k < 473 && k<strlen(batch.haps[r].bases) && k<batch.haps[r].length)
        {
            haps_ct_vector[st + k] = (int64_t)(batch.haps[r].bases[k]);
            k++;
        }

        haps_ct_vector[st + 473] = (int64_t)(batch.haps[r].length);
        // batch.haps_data_in_ct.push_back(encrypt_plaintext_vector_to_ciphertext(haps_ct_vector));
    }
    batch.haps_data_in_ct.push_back(encrypt_plaintext_vector_to_ciphertext(haps_ct_vector));

    batch.num_cells = total_read_length * total_hap_length;

    int_ct_vector[0] = (int64_t)(batch.num_reads);
    int_ct_vector[1] = (int64_t)(batch.num_haps);
    int_ct_vector[2] = (int64_t)(batch.id);
    int_ct_vector[3] = (int64_t)(batch.num_cells);

    batch.int_ct = encrypt_plaintext_vector_to_ciphertext(int_ct_vector);

    gettimeofday(&start_time, NULL);

    while (loop)
    {

        // #pragma omp parallel num_threads(nThreads)
        // {
        // int tid = omp_get_thread_num();
        int batch_size = batch.num_reads * batch.num_haps;
        assert(batch_size <= MAX_BATCH_SIZE);
        testcase *tc = &testcases[0];
        // For each batch, we will create testcases out of each of the different combinations of batch.reads[r] and batch.haps[h]
        // printf("tid: %d\n",tid);

        // #pragma omp for schedule(dynamic)           // The #pragma omp for directive in OpenMP is used to parallelize loops in a shared-memory parallel programming model. It allows the iterations of a loop to be executed in parallel by multiple threads.
        for (int r = 0; r < batch.num_reads; r++)
        {
            for (int h = 0; h < batch.num_haps; h++)
            {
                // printf("r: %d, h: %d\n",r,h);
                tc->rslen = batch.reads[r].length;
                tc->haplen = batch.haps[h].length;
                tc->hap = batch.haps[h].bases;
                tc->rs = batch.reads[r].bases;
                tc->q = batch.reads[r].q;
                tc->i = batch.reads[r].i;
                tc->d = batch.reads[r].d;
                tc->c = batch.reads[r].c;
                tc++;
            }
        }

        printf("---------------------------------------------------------------\n");
        printf("batch_count: %d\n",count);
        // printf("Going into computelikelihoods\n");
        // Use a separate testcase for each call
        computelikelihoodsboth(&testcases[0], batch.results, batch_size);
        // printf("Out of computeLikelihoods\n");

        // }
        loop--;
    }

    gettimeofday(&end_time, NULL);
    runtime += (end_time.tv_sec - start_time.tv_sec) * 1e6 + end_time.tv_usec - start_time.tv_usec;
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

void read_testfile(string filename, int nThreads, int loop)
{
    istream *is;
    ifstream ifs;
    if (filename == "")
    {
        printf("Reading test data from stdin");
        is = &std::cin;
    }
    else
    {
        // printf("Reading test data from file: %s\n", filename.c_str());
        ifs.open(filename.c_str());         // c_str() converts string to array of characters and terminates with \0
        if (!ifs.is_open())
        {
            printf("Cannot open file : %s", filename.c_str());
            exit(0);
        }
        is = &ifs;
    }

    // vector<Batch> batches;
    fprintf(stderr, "size of 1 batch: %lld\n", sizeof(Batch));
    int count = 0;

    while (!is->eof())
    {
        // cout<<"nThreads: "<<nThreads<<"; loop: "<<loop<<endl;
        printf("------------------\n");
        printf("count: %d\n", count);
        read_batch(*is, count, nThreads, loop);
        // batch.id = count++;
        count++;
        // batches.push_back(batch);
    }

    // return batches;
}

int main(int argc, char **argv)
{
#ifdef VTUNE_ANALYSIS
    __itt_pause();
#endif
    if (argc == 1)
    {
        std::cout << USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        // printf("c: %c, optarg: %s\n",c,optarg);
        switch (c)
        {
        case 'f':
            opt::testfile = optarg;
            break; // ../../input-datasets/phmm/small/5m.in
        case 'l':
            opt::loop = stoi(optarg);
            break;
        case 't':
            opt::nThreads = stoi(optarg);
            break; // 1
        case OPT_HELP:
            std::cout << USAGE_MESSAGE;
            exit(EXIT_SUCCESS);
        }
    }

    // disable stdout buffering
    setbuf(stdout, NULL);

    initPairHMM();

    ConvertChar::init();

    testcases = (testcase *)_mm_malloc(MAX_BATCH_SIZE * opt::nThreads * sizeof(testcase), 64);

    // vector<Batch> batches;
    read_testfile(opt::testfile, opt::nThreads, opt::loop);

    // printf("Num threads %d\n",opt::nThreads);

    // #pragma omp parallel num_threads(opt::nThreads)
    // {
    //     int tid = omp_get_thread_num();
    //     if (tid == 0) {
    //         printf("Running %d threads\n", opt::nThreads);
    //     }
    // }

    // struct timeval start_time, end_time;
    // double runtime = 0;

    // gettimeofday(&start_time, NULL);

#ifdef VTUNE_ANALYSIS
    __itt_resume();
#endif
    // cout<<"opt::loop: "<<opt::loop<<endl;

    // while (opt::loop) {
    //     #ifdef ENABLE_SORT
    //         std::sort(batches.begin(), batches.end(), SortByCells());
    //     #endif

    //     #pragma omp parallel num_threads(opt::nThreads)
    //     {
    //         int tid = omp_get_thread_num();
    //         #pragma omp for schedule(dynamic)           // The #pragma omp for directive in OpenMP is used to parallelize loops in a shared-memory parallel programming model. It allows the iterations of a loop to be executed in parallel by multiple threads.
    //             for (int i = 0; i < batches.size(); i++) {
    //                 int batch_size = batches[i].num_reads * batches[i].num_haps;
    //                 assert(batch_size <= MAX_BATCH_SIZE);
    //                 testcase *tc = &testcases[tid * MAX_BATCH_SIZE];
    //                 // For each batch, we will create testcases out of each of the different combinations of batch.reads[r] and batch.haps[h]

    //                 for (int r = 0; r < batches[i].num_reads; r++) {
    //                     for (int h = 0; h < batches[i].num_haps; h++) {
    //                         tc->rslen = batches[i].reads[r].length;
    //                         tc->haplen = batches[i].haps[h].length;
    //                         tc->hap = batches[i].haps[h].bases;
    //                         tc->rs = batches[i].reads[r].bases;
    //                         tc->q = batches[i].reads[r].q;
    //                         tc->i = batches[i].reads[r].i;
    //                         tc->d = batches[i].reads[r].d;
    //                         tc->c = batches[i].reads[r].c;
    //                         tc++;
    //                     }
    //                 }
    //                 // Use a separate testcase for each call
    //                 computelikelihoodsboth(&testcases[tid * MAX_BATCH_SIZE], batches[i].results, batch_size);
    //             }
    //     }
    //     opt::loop--;

#ifdef ENABLE_SORT
    std::sort(batches.begin(), batches.end(), SortById());
#endif
    // }
#ifdef VTUNE_ANALYSIS
    __itt_pause();
#endif

    // gettimeofday(&end_time, NULL);
    // runtime += (end_time.tv_sec - start_time.tv_sec)*1e6 + end_time.tv_usec - start_time.tv_usec;

    //     for (int i = 0; i < batches.size(); i++) {
    // #ifdef PRINT_OUTPUT
    //         int batch_size = batches[i].num_reads * batches[i].num_haps;
    //         for (int j = 0; j < batch_size; j++) {
    //             printf("%lf\n", batches[i].results[j]);
    //         }
    // #endif
    //         _mm_free(batches[i].reads);
    //         _mm_free(batches[i].haps);
    //         _mm_free(batches[i].results);
    //     }
    _mm_free(testcases);
    cout << "max_base_size: " << max_base_size << "; max_q_size: " << max_q_size << "; max_i_size: " << max_i_size << "; max_d_size: " << max_d_size << "; max_c_size: " << max_c_size << "; max_num_reads: " << max_num_reads << "; max_num_haps: " << max_num_haps << "; max_haps_bases_size: " << max_haps_bases_size << endl;
    printf("\nPairHMM completed. Kernel runtime: %.2f sec\n", runtime * 1e-6);
}
