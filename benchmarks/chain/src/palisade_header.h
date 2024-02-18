# include <iostream>
# include <stdio.h>
# include <stdlib.h>
# include <limits.h>
#include <signal.h>
#include "palisade.h"
#include <fstream>
#include <sstream>
#include <cmath>
#include <math.h>

#include "ciphertext-ser.h"
#include "cryptocontext-ser.h"
#include "pubkeylp-ser.h"
#include "scheme/bfvrns/bfvrns-ser.h"

// #define _CRTDBG_MAP_ALLOC
// #include <stdlib.h>
// #include <crtdbg.h>

using namespace std;
using namespace lbcrypto;

extern uint64_t p; // = 65537;
extern double sigma; //= 3.2; 
extern SecurityLevel securityLevel ; //= HEStd_128_classic;
extern uint32_t depth; // = 2;
extern CryptoContext<DCRTPoly> cc; // = CryptoContextFactory<DCRTPoly>::genCryptoContextBFVrns(p, securityLevel, sigma, 0, depth, 0, OPTIMIZED);
extern int32_t n; // = cc->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder()/2; 
extern LPKeyPair<DCRTPoly> kp; 
extern bool init_flag;
extern long sa_ms_byte_enc_counter;
extern long sa_ls_word_enc_counter;

//data types we will need
using CT = Ciphertext<DCRTPoly> ; //ciphertext
using PT = Plaintext ; //plaintext
using vecCT = vector<CT>; //vector of ciphertexts
using vecPT = vector<PT>; //vector of plaintexts
using vecInt = vector<int64_t>; // vector of ints
using vecChar = vector<char>; // vector of characters
 
void init();
vector<int64_t> decrypt_ciphertext_to_plaintext_vector(Ciphertext<DCRTPoly> ciphertext);
Ciphertext<DCRTPoly> encrypt_plaintext_integer_to_ciphertext(int64_t c);
Ciphertext<DCRTPoly> encrypt_plaintext_vector_to_ciphertext(vector<int64_t> d);
// Ciphertext<DCRTPoly> encrypt_plaintext_integer_to_ciphertext(int c);
char* convert_ciphertext_vector_to_plaintext_string(vecCT enc_v);
vector<int64_t> convert_ciphertext_vector_to_plaintext_vector(vecCT enc_v);
Plaintext encode_integer_to_plaintext(int64_t c);
Plaintext encode_integer_to_plaintext(int c);
Plaintext encode_vector_to_plaintext(vecInt c);
vector<CT> get_encrypted_bits_vector(int64_t n);
vector<CT> get_encrypted_bits_vector(int n);
CT shift_left(CT c, int n);
CT shift_encrypted_bit_vector_and_return_integer(vector<CT> encrypted_bit_vector, int64_t n);
int strlen_enc(vecCT v);
int64_t operate_and_decrypt(CT c1, string oper, CT c2);
int64_t operate_and_decrypt(CT c1, string oper, int64_t x);
int64_t compare_enc(CT c1, CT c2);
int64_t compare_enc(CT c1, int64_t x);
void assign_string_to_vecCT(vecCT &v, char* c, int length_to_be_assigned);

int64_t strcmp_enc(vecCT v, char* s);
int64_t strcmp_enc(vecCT v1, vecCT v2);
int64_t strcmp_enc(char* s, vecCT v);

void strdup_enc(vecCT s, vecCT& d);
void strcat_enc(vecCT &s,char* a, int index);
// CT deserialize_ciphertext_from_file(string ct_name);
Ciphertext<DCRTPoly> encrypt_plaintext_vector_to_ciphertext(vector<int64_t> d) ;
int strlen_string_enc(vecCT v);  
// CT serialize_ciphertext_to_file(string ct_name);

int64_t sa_ls_word_i(int64_t i);
int64_t sa_ms_byte_i(int64_t i);
CT cp_occ_cp_count_i(int64_t i, int64_t j);
CT cp_occ_one_hot_bwt_str_i(int64_t i, int64_t j);
string p_str_i(int64_t i);
int64_t compare_element_at_index_in_ct_and_other_element(CT c, int index, int64_t ele);
CT do_logical_and_of_encryted_bit_vectors(vecCT a, vecCT b);
