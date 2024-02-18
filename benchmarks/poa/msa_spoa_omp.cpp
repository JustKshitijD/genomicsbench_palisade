/* To complile:
   g++ -O3 msa_spoa_omp.cpp -fopenmp -std=c++11 -I include/ -L build/lib/ -lspoa -o msa_spoa
  */
/*
Each batch can have at max 105 strings of max length 857;
We create slots of 1024 size for each string in out CT, so 1 CT has 16384/1024=16 strings;
In that case, ceil(105/16)=7 CTs needed to cover all 105 strings of batch;
If you want ith string (indexing starts from 0), access it in (i/16)th CT and 
from (1024*(i%16)) index in vector of that CT till 0 is encountered.
*/ 

#include <getopt.h>
#include <stdio.h>
#include <cstdint>
#include <iostream>
#include <cstring>
#include <fstream>
#include <string>
#include <algorithm>
#include <sys/time.h>
#include <sys/resource.h>
#include <time.h>
#include <assert.h>
#include <exception>
#include <getopt.h>
#include <omp.h>
#include "palisade_header.h"                // This palisade_header.h is in same directory
#include "spoa/spoa.hpp"
#include "spoa/graph.hpp"
#include "spoa/alignment_engine.hpp"
#include <x86intrin.h>
#include <chrono>

using namespace std::chrono;

// #define VTUNE_ANALYSIS 1

#define CLMUL 8

// #define ENABLE_SORT 1

#ifdef VTUNE_ANALYSIS
    #include <ittnotify.h>
#endif

using namespace std;
using Alignment = std::vector<std::pair<std::int32_t, std::int32_t>>;

// #define DEBUG_FILE_READ
// #define PRINT_OUTPUT

typedef struct {
    int id;
    vector<string> seqs;
    // string consensus_seq;
    vecCT consensus_seq;
    uint8_t padding[24];
    vecCT enc_seqs{vecCT(7,0)};
} Batch;

struct SortBySize
{
    bool operator()( const Batch& lx, const Batch& rx ) const {
        return lx.seqs.size() > rx.seqs.size();
    }
};

struct SortById
{
    bool operator()( const Batch& lx, const Batch& rx ) const {
        return lx.id < rx.id;
    }
};

double get_realtime()
{
	struct timeval tp;
	struct timezone tzp;
	gettimeofday(&tp, &tzp);
	return tp.tv_sec*1e6 + tp.tv_usec;
}

long peakrss(void)
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
#ifdef __linux__
	return r.ru_maxrss * 1024;
#else
	return r.ru_maxrss;
#endif
}

//     typedef struct {
//     int id;
//     vector<string> seqs;     // - max length of a sequence here is 857; max size of this vector is 105
//     string consensus_seq;
//     uint8_t padding[24];
// } Batch;

// readFile(fp_seq, batches)
void readFile(ifstream& in_file, vector<Batch>& batches) {
    string seq; vecCT enc_seq;
    if (!in_file.eof()) {
        getline(in_file, seq); // header line
        // cout<<"header seq: "<<seq<<endl;
        enc_seq.resize((int)(seq.size())+1);
        for(int j=0;j<(int)(seq.size());j++)
            enc_seq[j]=encrypt_plaintext_integer_to_ciphertext(seq[j]);
        enc_seq[(int)(seq.size())]=encrypt_plaintext_integer_to_ciphertext('\0');
    }

    int max_seq_len=-100; int max_batch_size=-100;

    int count = 0;
    while (!in_file.eof()) { // loop across batches
        // if (seq[1] == '0') {
        if(compare_enc(enc_seq[1],'0'))
        {
            Batch b;
            // cout<<"count: "<<count<<endl;
            b.id = count++; int batch_seq_count=0; vecInt v_seq(16384,0); int batch_ct_count=0;

            // batch_seq_count stores # of sequences in batch
            // batch_ct_count stores for that batch, how many CTs have yet been created bunching together 16 sequences;
            // 16 sequences are stored at gaps of 1024 in 1 CT 
            while (!in_file.eof()) { // loop across sequences in batch
                
                // 16 sequences can fit in one batch
                if(batch_seq_count%16==0 && batch_seq_count!=0){
                    // cout<<"Filling batch_ct_count: "<<batch_ct_count<<endl;
                    b.enc_seqs[batch_ct_count]=encrypt_plaintext_vector_to_ciphertext(v_seq);   // b.enc_seqs is of size 7- will store 7 CTs
                    batch_ct_count++;
                    fill(v_seq.begin(),v_seq.end(),0);
                    batch_seq_count=0;
                }

                getline(in_file, seq); // sequence line
                // cout<<"seq: "<<seq<<endl;
                // cout<<"seq.size(): "<<(int)(seq.size())<<endl;
                // if((int)(seq.size())>max_seq_len)
                // {
                //     max_seq_len=(int)(seq.size());
                //     // cout<<"Changed, max_seq_len: "<<max_seq_len<<endl;
                // }

                // cout<<"max_seq_len: "<<max_seq_len<<endl;
                b.seqs.push_back(seq);

                // cout<<"batch_seq_count: "<<batch_seq_count<<endl;
                for(int k=0;k<(int)(seq.size());k++){
                    v_seq[1024*batch_seq_count+k]=seq[k];
                }
                batch_seq_count++;

                getline(in_file, seq); // header line
                if (seq[1] == '0') {                            // fill batch with seqs till 0 is encountered
                    // cout<<"Filling batch_ct_count: "<<batch_ct_count<<endl;
                    b.enc_seqs[batch_ct_count]=encrypt_plaintext_vector_to_ciphertext(v_seq);
                    batch_ct_count++;
                    fill(v_seq.begin(),v_seq.end(),0);

                    // cout<<"Final batch_seq_count: "<<batch_seq_count<<endl;
                    // cout<<"Final batch_ct_count: "<<batch_ct_count<<endl<<endl;

                    batches.push_back(b);
                    break;
                }
            }

            max_batch_size = max(max_batch_size,(int)(b.seqs.size()));
            if (in_file.eof()) {
                // cout<<"Filling batch_ct_count: "<<batch_ct_count<<endl;
                b.enc_seqs[batch_ct_count]=encrypt_plaintext_vector_to_ciphertext(v_seq);
                batch_ct_count++;
                fill(v_seq.begin(),v_seq.end(),0);

                // cout<<"Final batch_seq_count: "<<batch_seq_count<<endl;
                // cout<<"Final batch_ct_count: "<<batch_ct_count<<endl<<endl;
                batches.push_back(b);
            }
        }
    }

    // printf("max_seq_len: %d\n",max_seq_len);
    // printf("max_batch_size: %d\n",max_batch_size);

    // for (int i = 0; i < batches.size(); i++) {
    //     cout<<"------------------------------------------------------------"<<endl;
    //     cout << "Batch " << i << endl;
    //     for (int j = 0; j < batches[i].seqs.size(); j++) {
    //         // cout << batches[i].seqs[j] << endl;
            
    //         int ct_count=j/16; int ct_index=1024*(j%16);
    //         vecInt c=decrypt_ciphertext_to_plaintext_vector(batches[i].enc_seqs[ct_count]); string s2=""; int k;
    //         for(k=ct_index;k<ct_index+1024 && c[k]!=0;k++){
    //             // printf("%c",c[k]);
    //             s2+=c[k];
    //         }
    //         // cout<<endl;
    //         // cout<<"k: "<<k<<endl;
    //         // s2+='\0';
    //         if(batches[i].seqs[j].compare(s2)!=0){
    //             printf("%s\n",batches[i].seqs[j]);
    //             printf("%s\n",s2);
    //             printf("strlen(batches[i].seqs[j]): %d\n",batches[i].seqs[j].size());
    //             printf("strlen(c2): %d\n\n",s2.size());
    //         }
    //         //cout<<endl<<endl;
    //     }
    // }

#ifdef DEBUG_FILE_READ    
    for (int i = 0; i < batches.size(); i++) {
        cout << "Batch " << i << endl;
        for (int j = 0; j < batches[i].seqs.size(); j++) {
            cout << batches[i].seqs[j] << endl;
        }
    }
#endif

}


void help() {
    std::cout <<
        "\n"
        "usage: ./poa -s input.fasta -t <num_threads> > cons.fasta\n"
        "\n"
        "    options:\n"
        "        -m <int>\n"
        "            default: 2\n"
        "            score for matching bases\n"
        "        -x <int>\n"
        "            default: 4\n"
        "            penalty for mismatching bases\n"
        "        -o <int(,int)>\n"
        "            default: gap_open1 = 4, gap_open2 = 24\n"
        "            gap opening penalty (must be non-negative)\n"
        "        -e <int(,int)>\n"
        "            default: gap_ext1 = 2, gap_ext2 = 1\n"
        "            gap extension penalty (must be non-negative)\n"
        "        -s <file>\n"
        "            default: seq.fa\n"
        "            the input sequence set\n"
        "        -n <int>\n"
        "            default: 10\n"
        "            number of sequences in each set\n"
        "        -t <int>\n"
        "            default: 1\n"
        "            number of CPU threads\n"
        "        -h \n"
        "            prints the usage\n";
}

// ../benchmarks/poa/poa -s $INPUTS_DIR/poa/small/input-1000.fasta -t 1

int main(int argc, char** argv) {
#ifdef VTUNE_ANALYSIS
    __itt_pause();
#endif
    string seq_file = "seq.fa";

    std::uint8_t algorithm = 1;
    std::int8_t m = 2; CT enc_m=encrypt_plaintext_integer_to_ciphertext(m);
    std::int8_t x = -4; CT enc_x=encrypt_plaintext_integer_to_ciphertext(x);
    std::int8_t o1 = -4; CT enc_o1=encrypt_plaintext_integer_to_ciphertext(o1);
    std::int8_t e1 = -2; CT enc_e1=encrypt_plaintext_integer_to_ciphertext(e1);
    std::int8_t o2 = -24; CT enc_o2=encrypt_plaintext_integer_to_ciphertext(o2);
    std::int8_t e2 = -1; CT enc_e2=encrypt_plaintext_integer_to_ciphertext(e2);

    if (argc == 1) {
        help();
        exit(EXIT_FAILURE);
    }
    // long int strtol(const char *str, char **endptr, int base) converts the initial part of the string in str to a long int value according to the given base, which must be between 2 and 36 inclusive, or be the special value 0
    
    // getopt(int argc, char *const argv[], const char *optstring)
    // optstring is simply  a list of characters, 
    // each representing a single character option.

    char opt, *s; int n_seqs = 0, numThreads = 1; 
    CT enc_n_seqs=encrypt_plaintext_integer_to_ciphertext(0); CT enc_numThreads = encrypt_plaintext_integer_to_ciphertext(numThreads);

    while ((opt = getopt(argc, argv, "l:m:x:o:n:e:q:c:s:t:h")) != -1) {
        switch (opt) {
            case 'm': m = atoi(optarg); enc_m=encrypt_plaintext_integer_to_ciphertext(m); break;
            case 'x': x = 0-atoi(optarg); enc_x=encrypt_plaintext_integer_to_ciphertext(x); break;
            case 'o': o1 = 0-strtol(optarg, &s, 10); enc_o1=encrypt_plaintext_integer_to_ciphertext(o1); if (*s == ',') {o2 = 0-strtol(s+1, &s, 10); enc_o2=encrypt_plaintext_integer_to_ciphertext(o2); break;}
            case 'e': e1 = 0-strtol(optarg, &s, 10); enc_e1=encrypt_plaintext_integer_to_ciphertext(e1); if (*s == ',') {e2 = 0-strtol(s+1, &s, 10); enc_e2=encrypt_plaintext_integer_to_ciphertext(e2); break; }
            case 'n': n_seqs = atoi(optarg); enc_n_seqs=encrypt_plaintext_integer_to_ciphertext(n_seqs); break;
            case 's': seq_file = optarg; break;
            case 't': numThreads = atoi(optarg); enc_numThreads=encrypt_plaintext_integer_to_ciphertext(numThreads); break;
            case 'h': help(); return 0;
            default: help(); return 1;
        }
    }

    std::int8_t oe1=o1+e1, oe2=o2+e2;
    CT enc_oe1 = cc->EvalAdd(enc_o1,enc_e1);
    CT enc_oe2 = cc->EvalAdd(enc_o2,enc_e2);

    std::unique_ptr<spoa::AlignmentEngine> alignment_engine[numThreads];
    for (int i = 0; i < numThreads; i++) {
        printf("-----------------------------------------------\n");
        printf("i: %d\n",i);

        try {
            alignment_engine[i] = spoa::createAlignmentEngine(
                    static_cast<spoa::AlignmentType>(algorithm), m, encrypt_plaintext_integer_to_ciphertext(m), 
                    x, encrypt_plaintext_integer_to_ciphertext(x), oe1, encrypt_plaintext_integer_to_ciphertext(oe1), e1, 
                    encrypt_plaintext_integer_to_ciphertext(e1), oe2, encrypt_plaintext_integer_to_ciphertext(oe2), e2, encrypt_plaintext_integer_to_ciphertext(e2));
        } catch(std::invalid_argument& exception) {
            std::cerr << exception.what() << std::endl;
            return 1;
        }
    }

    ifstream fp_seq;
    fp_seq.open(seq_file, ios::in);
    assert(fp_seq.is_open());

    vector<Batch> batches;

    struct timeval start_time, end_time, t_start, t_end;
    double runtime = 0; int seq_i; string seq;
    double realtime = 0, real_start, real_end;
    double graphCreationTime = 0, alignTime = 0, addToGraphTime = 0, generateConsensusTime = 0;

#pragma omp parallel num_threads(numThreads) 
{
    int tid = omp_get_thread_num();
    if (tid == 0) {
        fprintf(stderr, "Running with threads: %d\n", numThreads);
    }
}

    readFile(fp_seq, batches);
    fprintf(stderr, "Number of batches: %lu, Size of batch struct %d\n", batches.size(), sizeof(Batch));
    int64_t workTicks[CLMUL * numThreads];
    std::memset(workTicks, 0, CLMUL * numThreads * sizeof(int64_t));
    gettimeofday(&start_time, NULL); real_start = get_realtime();

#ifdef VTUNE_ANALYSIS
    __itt_resume();
#endif

#ifdef ENABLE_SORT
    std::sort(batches.begin(), batches.end(), SortBySize());
#endif

#pragma omp parallel num_threads(numThreads)
{
    int tid = omp_get_thread_num();
    #pragma omp for schedule(dynamic, 1)
        for (int i = 0; i < batches.size(); i++) {
            auto start = high_resolution_clock::now();
            printf("###################################################################\n");
            printf("Batch #%d\n",i);

            int64_t st1 = __rdtsc();
            // gettimeofday(&t_start, NULL);
            auto createGraph_start = high_resolution_clock::now();
            auto graph = spoa::createGraph();
            auto createGraph_end = high_resolution_clock::now();
            auto createGraph_duration = duration_cast<microseconds>(createGraph_end  - createGraph_start);
            cout<<"For batch "<<i<<"; createGraph_duration: "<<createGraph_duration.count()<<endl;
            // auto graph = spoa::createGraph();               // new graph for each batch
            // printf("After createGraph\n");

            // printf("coder_: \n");
            // printf("graph->coder_size(): %d\n",(int)(graph->coder_size()));
            // for(int k=0;k<graph->coder_size();k++){
            //     printf("%d ",(int)(graph->coder(k)));
            // }
            // printf("\n");

            // graph->print_codes_size_and_codes();
            // printf("num_codes: %d\n",graph->num_codes());           // initially 0

            // gettimeofday(&t_end, NULL);
            // graphCreationTime += (t_end.tv_sec - t_start.tv_sec)*1e6 + t_end.tv_usec - t_start.tv_usec;
            for (int j = 0; j < batches[i].seqs.size(); j++) {
                // printf("----------------------------------------------------------------------------------\n");
                
                // gettimeofday(&t_start, NULL);

                // printf("\nsequence number j: %d\n",j);
                // cout<<"Calling align - "<<endl;

                // std::unique_ptr<spoa::AlignmentEngine> alignment_engine[numThreads];
                // auto alignment = alignment_engine[tid]->align(batches[i].seqs[j], graph); 
                // printf("After align\n");

                // Initially, this too will have codes_ as all -1; this will change in add_alignment, where we'll do- 
                // codes_['A']=num_codes_; num_codes_++;
                // codes_['T']=num_codes_; num_codes_++;
                // codes_['G']=num_codes_; num_codes_++;
                // codes_['C']=num_codes_; num_codes_++;

                // sgraph->print_codes_size_and_codes();   

                // printf("coder_: \n");
                // printf("graph->coder_size(): %d\n",(int)(graph->coder_size()));
                // for(int k=0;k<graph->coder_size();k++){
                //     // printf("%d ",(int)(graph->coder(k)));
                //     graph->coder(k);
                // }
                // printf("\n");
                // printf("num_codes: %d\n",graph->num_codes());

                // cout<<"alignment.size(): "<<alignment.size()<<endl;

                // for(int k=0;k<alignment.size();k++){
                //     cout<<"k: "<<k<<endl;
                //     cout<<"("<<alignment[k].first<<","<<alignment[k].second<<")"<<endl;
                // }
                // gettimeofday(&t_end, NULL);
                // alignTime += (t_end.tv_sec - t_start.tv_sec)*1e6 + t_end.tv_usec - t_start.tv_usec;
                // gettimeofday(&t_start, NULL);

                // Retrieve batches[i].seqs[j] string and encrypt each character of string to ciphertext and push to vector
                // Indexing into this vector and decrypting it will give us each character

                auto indexing_start = high_resolution_clock::now();
                vecCT v_seq;
                int ct_count=j/16; int ct_index=1024*(j%16);
                vecInt c=decrypt_ciphertext_to_plaintext_vector(batches[i].enc_seqs[ct_count]); int k;
                for(k=ct_index;k<ct_index+1024 && c[k]!=0;k++){
                    // printf("%c",c[k]);
                    // s2+=c[k];
                    v_seq.push_back(encrypt_plaintext_integer_to_ciphertext(c[k]));
                }
                auto indexing_end = high_resolution_clock::now();
                auto indexing_duration = duration_cast<microseconds>(indexing_end  - indexing_start);
                cout<<"For sequence "<<j<<"; indexing_duration: "<<indexing_duration.count()<<endl;

                // using Alignment = std::vector<std::pair<std::int32_t, std::int32_t>>;

                auto add_alignment_start = high_resolution_clock::now();
                auto alignment = alignment_engine[tid]->align(batches[i].seqs[j], graph); 
                // cout<<"Entering add_alignment"<<endl;
                // auto add_alignment_start = high_resolution_clock::now();
                graph->add_alignment(alignment, batches[i].seqs[j], v_seq);  // default 3rd argument is uint32_t weight=1U
                auto add_alignment_end = high_resolution_clock::now();
                auto add_alignment_duration = duration_cast<microseconds>(add_alignment_end  - add_alignment_start);
                cout<<"For sequence "<<j<<"; add_alignment_duration: "<<add_alignment_duration.count()<<endl;
                // cout<<"Left add_alignment"<<endl;

                // gettimeofday(&t_end, NULL);
                // addToGraphTime += (t_end.tv_sec - t_start.tv_sec)*1e6 + t_end.tv_usec - t_start.tv_usec;
            }

            cout<<endl;
            // gettimeofday(&t_start, NULL);
            auto consensus_seq_start = high_resolution_clock::now();
            batches[i].consensus_seq = graph->generate_consensus();
            auto consensus_seq_end = high_resolution_clock::now();
            auto consensus_seq_duration = duration_cast<microseconds>(consensus_seq_end  - consensus_seq_start);
            cout<<"For batch "<<i<<"; consensus_seq_duration: "<<consensus_seq_duration.count()<<endl;

            batches[i].consensus_seq.emplace_back(encrypt_plaintext_integer_to_ciphertext('\0'));
            printf("batches[%d].consensus_seq: %s\n",i,convert_ciphertext_vector_to_plaintext_string(batches[i].consensus_seq));

            auto stop=high_resolution_clock::now();
            auto duration = duration_cast<microseconds>(stop - start);
            cout<<"Batch duration: "<<duration.count()<<endl;

            // gettimeofday(&t_end, NULL);
            // generateConsensusTime += (t_end.tv_sec - t_start.tv_sec)*1e6 + t_end.tv_usec - t_start.tv_usec;
            int64_t et1 = __rdtsc();
            workTicks[CLMUL * tid] += (et1 - st1);
        }
    printf("%d] workTicks = %ld\n", tid, workTicks[CLMUL * tid]);

}
#ifdef ENABLE_SORT
    std::sort(batches.begin(), batches.end(), SortById());
#endif

#ifdef VTUNE_ANALYSIS
    __itt_pause();
#endif
    gettimeofday(&end_time, NULL); real_end = get_realtime();
    runtime += (end_time.tv_sec - start_time.tv_sec)*1e6 + end_time.tv_usec - start_time.tv_usec;
    realtime += (real_end-real_start);

    int64_t sumTicks = 0;
    int64_t maxTicks = 0;
    for(int i = 0; i < numThreads; i++)
    {
        sumTicks += workTicks[CLMUL * i];
        if(workTicks[CLMUL * i] > maxTicks) maxTicks = workTicks[CLMUL * i];
    }
    double avgTicks = (sumTicks * 1.0) / numThreads;
    printf("avgTicks = %lf, maxTicks = %ld, load imbalance = %lf\n", avgTicks, maxTicks, maxTicks/avgTicks);
#ifdef PRINT_OUTPUT
    for (int i = 0; i < batches.size(); i++) {
        cout << ">Consensus_sequence" << endl;
        cout << batches[i].consensus_seq.c_str() << endl;
    }
#endif

    fprintf(stderr, "Runtime: %.2f, GraphCreate: %.2f, Align: %.2f, AddSeqGraph: %.2f, Consensus %.2f %.2f %.3f \n", runtime*1e-6, graphCreationTime*1e-6, alignTime*1e-6, addToGraphTime*1e-6, generateConsensusTime*1e-6, realtime*1e-6, peakrss()/1024.0/1024.0);

    fp_seq.close();
    return 0;
}
