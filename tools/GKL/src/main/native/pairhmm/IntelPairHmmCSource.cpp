#include <avx.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string.h>
#include <xmmintrin.h>
#include "IntelPairHmm.h"
#include "pairhmm_common.h"
#include "avx_impl.h"
#ifndef __APPLE__
  #include "avx512_impl.h"
#endif
#include "Context.h"
#include"debug.h"
bool g_use_double;
int g_max_threads;
bool g_use_fpga;

Context<float> g_ctxf;
Context<double> g_ctxd;

// Overall, this line declares a function pointer g_compute_full_prob_float that can point to a function returning a float and taking a testcase* parameter. 
float (*g_compute_full_prob_float)(testcase *tc);
double (*g_compute_full_prob_double)(testcase *tc);

/*
 * Method:    initNative
 */
void initPairHMM()
{
  // set function pointers
  if(is_avx512_supported())
  {
    printf("AVX512 supported!\n");
#ifndef __APPLE__
    DBG("Using CPU-supported AVX-512 instructions");
    g_compute_full_prob_float = compute_fp_avx512s;
    g_compute_full_prob_double = compute_fp_avx512d;
#else
    assert(false);
#endif
  }
  else
  {
    printf("AVX512 not supported!\n");
    g_compute_full_prob_float = compute_fp_avxs;
    g_compute_full_prob_double = compute_fp_avxd;
  }
  // init convert char table
  ConvertChar::init();
}

// typedef struct {
//   int rslen, haplen;
//   const char *q, *i, *d, *c;
//   const char *hap, *rs;
// } testcase;

// computelikelihoodsboth(&testcases[0], batch.results, batch_size);
// batch_size=batch.num_reads*batch.num_haps
void computelikelihoodsboth(testcase *testcases, double *expected_results, int batch_size)
{

#ifdef _OPENMP
  #pragma omp parallel for schedule(dynamic, 1)
#endif
  // printf("batch_size: %d\n",batch_size);

  for (int i = 0; i < batch_size; i++) {
    double result_final = 0;
    float result_float = g_compute_full_prob_float(&testcases[i]);
    // printf("result_float: %f\n",result_float);
    if (result_float < MIN_ACCEPTED) {
      double result_double = g_compute_full_prob_double(&testcases[i]);
      result_final = log10(result_double) - g_ctxd.LOG10_INITIAL_CONSTANT;
    }
    else {
      result_final = (double)(log10f(result_float) - g_ctxf.LOG10_INITIAL_CONSTANT);
    }
    printf("i: %d; result_final: %f\n",i,result_final);
    expected_results[i] = result_final;
  }
  
  return;
}

/* Computelikelihoodsfloat method takes in one vector of input and computes the pHMM result float */

void computelikelihoodsfloat(testcase *testcases, float* expected_results) {
    float result_final = 0;

    float result_float = g_compute_full_prob_float(testcases);

    result_final = (double)(log10f(result_float) - g_ctxf.LOG10_INITIAL_CONSTANT);

    (*expected_results) = (float)result_final;

    return;
}

/* Computelikelihoodsdouble method takes in one vector of input and computes the pHMM result double */

void computelikelihoodsdouble(testcase *testcases, double *expected_results)
{

    double result_final = 0;
    double result_double = g_compute_full_prob_double(testcases);
    result_final = log10(result_double) - g_ctxd.LOG10_INITIAL_CONSTANT;

    (*expected_results) = (double)result_final;

    return;

}

