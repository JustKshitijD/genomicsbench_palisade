#include "palisade_header.h"

//  Among other arguments, this factory method takes the plaintext modulus, the security level, and the additive, multiplicative,
// and key-switching depths of the computation. The plaintext modulus determines an upper
// bound for the integer inputs in the BFV scheme. The security level is a parameter whose
// possible values correspond to 128-, 192-, or 256-bit security; PALISADE picks the appropriate underlying security parameters by consulting the relevant tables of the current FHE
// standard. Depending on the computation the user wants to carry out, she can set either the
// additive, multiplicative, or key-switching depth, and set all others to 0

// uint64_t p = 65537;
// uint64_t p = 12869861377; // -  used for fmi, bsw, phmm
uint64_t p=3276802621441;               // 3.2*10^12
// uint64_t p=9899511922689010000; 
double sigma = 3.2;
SecurityLevel securityLevel = HEStd_128_classic;
uint32_t depth = 3;
CryptoContext<DCRTPoly> cc = CryptoContextFactory<DCRTPoly>::genCryptoContextBFVrns(p, securityLevel, sigma, 0, depth, 0, OPTIMIZED);
int32_t n = cc->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder() / 2;
// double log2_q = std::log2(cc->GetCryptoParameters()->GetElementParams()->GetModulus().ConvertToDouble());
LPKeyPair<DCRTPoly> kp;
bool init_flag = false;
long sa_ms_byte_enc_counter = 83443;
long sa_ls_word_enc_counter = 83443;

void init()
{
    if (init_flag == false)
    {
        printf("p: %lu\n",cc->GetCryptoParameters()->GetPlaintextModulus());        
        printf("n: %d\n",n);             // 16384
        // printf("log2_q: %lf\n",log2_q);  // 240
     
        // if (!Serial::DeserializeFromFile("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/serializer_folder/cryptocontext.txt", cc,
        //                                  SerType::BINARY))
        // {
        //     std::cerr << "I cannot read serialization from "
        //               << "/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/serializer_folder/cryptocontext.txt" << std::endl;
        //     return;
        // }
        // std::cout << "The cryptocontext has been deserialized." << std::endl;

        // LPPublicKey<DCRTPoly> pk;
        // if (Serial::DeserializeFromFile("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/serializer_folder/key-public.txt", pk,
        //                                 SerType::BINARY) == false)
        // {
        //     std::cerr << "Could not read public key" << std::endl;
        //     return;
        // }
        // std::cout << "The public key has been deserialized." << std::endl;

        // std::ifstream emkeys("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/serializer_folder/key-eval-mult.txt",
        //                      std::ios::in | std::ios::binary);
        // if (!emkeys.is_open())
        // {
        //     std::cerr << "I cannot read serialization from "
        //               << "/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/serializer_folder/key-eval-mult.txt" << std::endl;
        //     return;
        // }
        // if (cc->DeserializeEvalMultKey(emkeys, SerType::BINARY) == false)
        // {
        //     std::cerr << "Could not deserialize the eval mult key file" << std::endl;
        //     return;
        // }
        // std::cout << "Deserialized the eval mult keys." << std::endl;

        // LPPrivateKey<DCRTPoly> sk;
        // if (Serial::DeserializeFromFile("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/serializer_folder/key-private.txt", sk,
        //                                 SerType::BINARY) == false)
        // {
        //     std::cerr << "Could not read secret key" << std::endl;
        //     return;
        // }
        // std::cout << "The secret key has been deserialized." << std::endl;

        cc->Enable(ENCRYPTION);
        cc->Enable(SHE);

        kp=cc->KeyGen();
        cc->EvalMultKeyGen(kp.secretKey);
        cc->EvalAtIndexKeyGen(kp.secretKey, {1, 2, -1, -2});
        // kp.publicKey = pk;
        // kp.secretKey = sk;
        init_flag = true;

        // // Serialize cryptocontext
        // if (!Serial::SerializeToFile("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/serializer_folder/cryptocontext.txt", cc,
        //                             SerType::BINARY)) {
        //     std::cerr << "Error writing serialization of the crypto context to "
        //                 "cryptocontext.txt"
        //             << std::endl;
        //     return;
        // }
        // std::cout << "The cryptocontext has been serialized." << std::endl;

        // // Serialize the public key
        // if (!Serial::SerializeToFile("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/serializer_folder/key-public.txt",
        //                             kp.publicKey, SerType::BINARY)) {
        //     std::cerr << "Error writing serialization of public key to key-public.txt"
        //             << std::endl;
        //     return;
        // }
        // std::cout << "The public key has been serialized." << std::endl;

        // // Serialize the secret key
        // if (!Serial::SerializeToFile("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/serializer_folder/key-private.txt",
        //                             kp.secretKey, SerType::BINARY)) {
        //     std::cerr << "Error writing serialization of private key to key-private.txt"
        //             << std::endl;
        //     return;
        // }
        // std::cout << "The secret key has been serialized." << std::endl;

        // // Serialize the relinearization (evaluation) key for homomorphic
        // // multiplication
        // std::ofstream emkeyfile("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/serializer_folder/key-eval-mult.txt",
        //                         std::ios::out | std::ios::binary);
        // if (emkeyfile.is_open()) {
        //     if (cc->SerializeEvalMultKey(emkeyfile, SerType::BINARY) == false) {
        //     std::cerr << "Error writing serialization of the eval mult keys to "
        //                 "key-eval-mult.txt"
        //                 << std::endl;
        //     return;
        //     }
        //     std::cout << "The eval mult keys have been serialized." << std::endl;

        //     emkeyfile.close();
        // } else {
        //     std::cerr << "Error serializing eval mult keys" << std::endl;
        //     return;
        // }
    }
}

vecInt decrypt_ciphertext_to_plaintext_vector(CT ciphertext)
{
    // printf("In decrypt_ciphertext_to_plaintext_vector\n");
    Plaintext plaintext;
    cc->Decrypt(kp.secretKey, ciphertext, &plaintext);
    vecInt v;
    // printf("Earlier v.size(): %d\n",v.size());
    // for(int i=0;i<v.size();i++){
    //     cout<<"v["<<i<<"] = "<<v[i]<<endl;
    // }

    v = plaintext->GetPackedValue();
    
    // printf("Later v.size(): %d\n",v.size());
    // for(int i=0;i<v.size();i++){
    //     cout<<"v["<<i<<"] = "<<v[i]<<endl;
    // }
    // cout<<endl;

    return v;
}

// /* Put integer as single element in a plaintext vetor, then encrypt that vector and return the ciphertext */
// Ciphertext<DCRTPoly> encrypt_plaintext_integer_to_ciphertext(int d)
// {
//     // printf("init_flag: %d\n",init_flag);
//     init();
//     int64_t c=d;
//     vector<int64_t> vec;
//     vec.clear();
//     vec.push_back(c);
//     Plaintext plaintext = cc->MakePackedPlaintext(vec);
//     auto ciphertext = cc->Encrypt(kp.publicKey, plaintext);
//     return ciphertext;
// }

/* Put integer as single element in a plaintext vetor, then encrypt that vector and return the ciphertext */
Ciphertext<DCRTPoly> encrypt_plaintext_integer_to_ciphertext(int64_t d)
{
    // printf("init_flag: %d\n",init_flag);
    init();
    // cout<<"In encrypt_plaintext_integer_to_ciphertext, d: "<<d<<endl;
    int64_t c = d;
    vector<int64_t> vec;
    vec.clear();
    vec.push_back(c);
    // cout<<"c: "<<c<<endl;
    Plaintext plaintext = cc->MakePackedPlaintext(vec);
    auto ciphertext = cc->Encrypt(kp.publicKey, plaintext);
    return ciphertext;
}

/* Put integer as single element in a plaintext vetor, then encrypt that vector and return the ciphertext */
Ciphertext<DCRTPoly> encrypt_plaintext_vector_to_ciphertext(vector<int64_t> d)
{
    // printf("init_flag: %d\n",init_flag);
    init();
    // int64_t c=d;
    // vector<int64_t> vec;
    // vec.clear();
    // vec.push_back(c);
    // // cout<<"c: "<<c<<endl;
    Plaintext plaintext = cc->MakePackedPlaintext(d);
    auto ciphertext = cc->Encrypt(kp.publicKey, plaintext);
    return ciphertext;
}

Plaintext encode_integer_to_plaintext(int64_t c)
{
    // printf("init_flag: %d\n",init_flag);
    init();
    vector<int64_t> vec;
    vec.clear();
    vec.push_back(c);
    Plaintext plaintext = cc->MakePackedPlaintext(vec);
    return plaintext;
}

Plaintext encode_vector_to_plaintext(vecInt c)
{
    // printf("init_flag: %d\n",init_flag);
    init();
    Plaintext plaintext = cc->MakePackedPlaintext(c);
    return plaintext;
}

Plaintext encode_integer_to_plaintext(int d)
{
    // printf("init_flag: %d\n",init_flag);
    init();
    int64_t c = (int64_t)(d);
    vector<int64_t> vec;
    vec.clear();
    vec.push_back(c);
    Plaintext plaintext = cc->MakePackedPlaintext(vec);
    return plaintext;
}

CT shift_left(CT c, int n) // easy shift left encrypted bit vector
{
    for (int i = 0; i < n; i++)
        c = cc->EvalAdd(c, c);
    return c;
}

// n is -ve means shift the vector with encrypted 0s and 1s to right; n is +vs means shift left
CT shift_encrypted_bit_vector_and_return_integer(vector<CT> encrypted_bit_vector, int64_t n)
{
    // cout<<"n: "<<n<<"; encrypted_bit_vector.size(): "<<encrypted_bit_vector.size()<<endl;

    // cout<<"encrypted_bit_vector: [";
    // for(int i=0;i<encrypted_bit_vector.size();i++)
    //     cout<<decrypt_ciphertext_to_plaintext_vector(encrypted_bit_vector[i])[0]<<" ";
    // cout<<"]"<<endl;

    CT ans = encrypt_plaintext_integer_to_ciphertext(0);

    if (n > 0)
    {
        int k = 0;
        for (int i = encrypted_bit_vector.size() + n - 1; i >= n; i--)
        {
            CT c = encrypt_plaintext_integer_to_ciphertext(0);
            for (int j = 0; j < i; j++)
                c = cc->EvalAdd(c, cc->EvalAdd(encrypted_bit_vector[k], encrypted_bit_vector[k]));

            ans = cc->EvalAdd(ans, c);
            k++;
        }

        return ans;
    }

    if (n < 0 && -1 * n > encrypted_bit_vector.size())
        return encrypt_plaintext_integer_to_ciphertext(0);

    CT c = encrypt_plaintext_integer_to_ciphertext(0);
    int k = 0;
    for (int i = encrypted_bit_vector.size() - (-1 * n) - 1; i >= 0; i--)
    {
        // CT c=encrypt_plaintext_integer_to_ciphertext(0);
        CT c = encrypted_bit_vector[k];

        if (decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(c, encode_integer_to_plaintext(0)))[0] == 0)
        {
            k++;
            continue;
        }

        for (int j = 0; j < i; j++)
            c = cc->EvalAdd(c, c);

        // cout<<"c: "<<decrypt_ciphertext_to_plaintext_vector(c)[0]<<endl;
        ans = cc->EvalAdd(ans, c);
        k++;

        // cout << "i: " << i << "; c: " << decrypt_ciphertext_to_plaintext_vector(c)[0] << "; ans: " << decrypt_ciphertext_to_plaintext_vector(ans)[0] << endl;
    }

    return ans;
}

vector<CT> get_encrypted_bits_vector(int64_t n)
{
    vector<CT> v;
    do
    {
        int64_t x = (int64_t)(n % 2);
        CT w = encrypt_plaintext_integer_to_ciphertext(x);
        v.insert(v.begin(), w);
        n /= 2;
    }while (n != 0);

    return v;                       // sizeof(v) is # of bits in bit-decomposition of n
}

vector<CT> get_encrypted_bits_vector(int m) // makes an encrypted vector out of the bits comprising input int m
{
    int64_t n = (int64_t)(m);
    vector<CT> v;
    
    do
    {
        int64_t x = (int64_t)(n % 2);
        CT w = encrypt_plaintext_integer_to_ciphertext(x);
        v.insert(v.begin(), w);
        n /= 2;
    }while (n != 0);

    return v;
}

CT do_logical_and_of_encryted_bit_vectors(vecCT a, vecCT b)
{
    // printf("a.size(): %d, b.size(): %d\n",(int)(a.size()),(int)(b.size()));

    vecCT v(max(a.size(),b.size()));
    vecCT v1,v2;                            // v1 has larger array; v2 has smaller array
    if(a.size()>b.size())
    {
        v1.resize(a.size()); v2.resize(b.size());
        for(int i=0;i<a.size();i++){
            v1[i]=a[i];
        }
        for(int i=0;i<b.size();i++){
            v2[i]=b[i];
        }
    }
    else{
        v1.resize(b.size()); v2.resize(a.size());
        for(int i=0;i<b.size();i++){
            v1[i]=b[i];
        }
        for(int i=0;i<a.size();i++){
            v2[i]=a[i];
        }
    }

    // for(int i=0;i<v1.size();i++){
    //     cout<<"i: "<<i<<"; v1[i]: "<<decrypt_ciphertext_to_plaintext_vector(v1[i])[0]<<endl;
    // }
    // cout<<endl;
    // for(int i=0;i<v2.size();i++){
    //     cout<<"i: "<<i<<"; v2[i]: "<<decrypt_ciphertext_to_plaintext_vector(v2[i])[0]<<endl;
    // }
    // cout<<endl;

    for(int i=0;i<v1.size()-v2.size();i++){
        v[i]=encrypt_plaintext_integer_to_ciphertext(0);
        // printf("v1[%d]: %lld, v2[%d]:%lld, v[%d]: %lld\n",i,v1[i],i,v2[i],i,v[i]);
    }

    for(int i=v1.size()-v2.size();i<v1.size();i++){
        CT c1=v1[i];
        CT c2=v2[i-(v1.size()-v2.size())];

        // if(decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(c1,c2))[0]!=0)
        //     v[i]=encrypt_plaintext_integer_to_ciphertext(0);
        // else
        // {
        //     if(cc->Eval)
        // }       

        v[i]=cc->EvalMult(c1,c2);               // Mult(c1,c2) gives c1 AND c2
        // printf("v1[%d]: %lld, v2[%d]:%lld, v[%d]: %lld\n",i,v1[i],i,v2[i],i,v[i]);
    }

    // for(int i=0;i<v.size();i++){
    //     cout<<"i: "<<i<<"; v[i]: "<<decrypt_ciphertext_to_plaintext_vector(v[i])[0]<<endl;
    // }
    // cout<<endl;

    CT ans=encrypt_plaintext_integer_to_ciphertext(0);

    for(int i=0;i<v.size();i++)
    {
        CT c = v[i];

        if (decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(c, encode_integer_to_plaintext(0)))[0] == 0)
        {
            // k++;
            continue;
        }

        // for (int j = 0; j < v.size()-1-i; j++)
        //     c = cc->EvalAdd(c, c);

        // // cout<<"c: "<<decrypt_ciphertext_to_plaintext_vector(c)[0]<<endl;
        // ans = cc->EvalAdd(ans, c);

        ans=cc->EvalAdd(ans,cc->EvalMult(v[i],encode_integer_to_plaintext((int64_t)(pow(2,v.size()-i-1)))));
    }

    return ans;
}

vector<int64_t> convert_ciphertext_vector_to_plaintext_vector(vecCT enc_v) // decrypts each element of enc_v into plaintext vector, but leaves the final NULL values
{
    vector<int64_t> v;
    v.resize(enc_v.size());
    // cout<<"enc_v.size(): "<<enc_v.size()<<endl;
    int i = 0;
    while (i < enc_v.size())
    {
        if (enc_v[i])
            v[i] = decrypt_ciphertext_to_plaintext_vector(enc_v[i])[0];
        else
            v[i] = 0;
        i++;
    }
    return v;
}

// It is assumed that enc_v already has enc('\0') as final element; so if enc_v=[enc('a'),enc('b'),enc('\0'),NULL,NULL], then, strlen_enc(enc_v)=3
char *convert_ciphertext_vector_to_plaintext_string(vecCT enc_v) // create char* array of same size as vecCT and copy the initial non-NULL characters
{
    // strlen_enc(enc_v) gives size of enc_v till \0 is hit
    char *s = (char *)malloc(strlen_enc(enc_v));
    int i = 0;
    while (i < enc_v.size() && enc_v[i] != 0)
    {
        s[i] = (char)(decrypt_ciphertext_to_plaintext_vector(enc_v[i])[0]);
        i++;
    }
    return s;
}

// count characters in vecCT till NULL(not encrypted('\0)) is hit
// But what if v has some enc('\0') in continuation and NULL only after that?
// So, we should return sz
int strlen_enc(vecCT v) // get length of vecCT till NULL character is hit
{
    int sz = 0;
    for (int i = 0; i < v.size(); i++)
    {
        if (v[i] != 0)
            sz++;
        else
            break;
    }

    return sz;
}

// count characters in vecCT till enc('\0') or NULL is hit;
// When str is converted to enc_str, as we will use a for-loop to convert each character of str to CT and put into vecCT, the final '\0' won't be there in vecCT and all trailing slots of vecCT will be enc(0)
int strlen_string_enc(vecCT v) // get length of vecCT till NULL character is hit
{
    // cout << "In strlen_string_enc" << endl;
    int sz = 0;
    for (int i = 0; i < v.size(); i++)
    {
        // cout << "decrypt_ciphertext_to_plaintext_vector(v[i])[0]: " << decrypt_ciphertext_to_plaintext_vector(v[i])[0] << endl;
        if (operate_and_decrypt(v[i], "-", '\0') != 0 && v[i]!=0)
            sz++;
        else
            break;
    }

    return sz;
}

int64_t operate_and_decrypt(CT c1, string oper, CT c2)
{
    CT res;
    if (oper == "-")
        res = cc->EvalSub(c1, c2);
    else if (oper == "+")
        res = cc->EvalAdd(c1, c2);
    else if (oper == "*")
        res = cc->EvalMult(c1, c2);

    return decrypt_ciphertext_to_plaintext_vector(res)[0];
}

int64_t operate_and_decrypt(CT c1, string oper, int64_t x)
{
    CT res;
    if (oper == "-")
        res = cc->EvalSub(c1, encode_integer_to_plaintext(x));
    else if (oper == "+")
        res = cc->EvalAdd(c1, encode_integer_to_plaintext(x));
    else if (oper == "*")
        res = cc->EvalMult(c1, encode_integer_to_plaintext(x));

    return decrypt_ciphertext_to_plaintext_vector(res)[0];
}

// If c1 and c2 are same, returns 1
int64_t compare_enc(CT c1, CT c2)
{
    // printf("Comparing %d and %d\n",decrypt_ciphertext_to_plaintext_vector(c1)[0], decrypt_ciphertext_to_plaintext_vector(c2)[0]);
    // cout<<"res: "<<operate_and_decrypt(c1,"-",c2)<<endl;
    if (operate_and_decrypt(c1, "-", c2) == 0)
        return 1;
    return 0;
}

int64_t compare_enc(CT c1, int64_t x)
{
    // printf("2_Comparing %d and %d\n",decrypt_ciphertext_to_plaintext_vector(c1)[0], x);
    // cout<<"res: "<<operate_and_decrypt(c1,"-",x)<<endl;
    if (operate_and_decrypt(c1, "-", x) == 0)
        return 1;
    return 0;
}

// Assume that c encrypts vector of size 16384
int64_t compare_element_at_index_in_ct_and_other_element(CT c1, int index, int64_t ele){
    vector<int64_t> v2(16384,0);
    v2[index]=ele;
    CT c2=encrypt_plaintext_vector_to_ciphertext(v2);
    if(decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(c1,c2))[index]==0)
        return 1;
    return 0;
}

void assign_string_to_vecCT(vecCT &v, char *c, int length_to_be_assigned) // copy string to vecCT; last char will be made \0
{
    if (length_to_be_assigned == -1)
    {
        length_to_be_assigned = strlen(c);
    }

    if(v.size()<length_to_be_assigned+1)
        v.resize(length_to_be_assigned+1);

    // if (length_to_be_assigned == strlen(c))
    //     if (v.size() < strlen(c) + 1) // only resize if v doesn't have room for char* c
    //         v.resize(strlen(c) + 1);  // v.shrink_to_fit();                     // make room to add enc('\0) at last of vector
    // else                          // assign a part of c to v, but still terminate it with \0
    //     v.resize(length_to_be_assigned + 1);

    // v.reserve(length_to_be_assigned);
    int i = 0;

    for (i = 0; i < length_to_be_assigned; i++)
    {
        v[i] = encrypt_plaintext_integer_to_ciphertext(c[i]);
    }
    v[length_to_be_assigned] = encrypt_plaintext_integer_to_ciphertext(0);
}

void strcat_enc(vecCT &s, char *a, int index)
{
    try
    {
        if (a == NULL)
            throw(0);
    }
    catch (int x)
    {
        cout << "a is NULL in strcat_enc" << endl;
    }

    int sz = strlen_enc(s);
    try
    {
        if (sz - 1 + strlen(a) > s.size())
            throw(sz - 1 + strlen(a));
    }
    catch (int x)
    {
        cout << "Total concatenated size= " << x << ", is greater that s.size()= " << s.size() << endl;
    }

    if (index == -1)
        index = strlen_enc(s) - 1; // s is most likely enrypting a string; In that case, s[s.size()-1]=enc('\n'). So, a[0] should be written to s[s.size()-1], i.e, should overwrite enc('\n')
    int i = index;
    while (i - index < strlen(a))
    {
        s[i] = encrypt_plaintext_integer_to_ciphertext(a[i - index]);
        i++;
    }

    s[i] = encrypt_plaintext_integer_to_ciphertext('\0');
}

void strdup_enc(vecCT s, vecCT &d)
{
    int sz = s.size();
    // cout<<"sz: "<<sz<<endl;

    d.resize(sz);
    // cout<<"d.size(): "<<d.size()<<endl;

    // if(d.size()<s.size())                                                                   // only resize d if d doesn't have enough space for s
    //     d.resize(sz); // d.shrink_to_fit();

    int i = 0;
    for (i = 0; i < sz; i++)
    {
        // if(s[i]!=0)
        //     cout<<"s["<<i<<"]="<<decrypt_ciphertext_to_plaintext_vector(s[i])[0]<<endl;
        // else
        //     cout<<"Empty s[i]"<<endl; 
        d[i] = s[i];
    }
    while (i < d.size())
    {
        d[i] = 0;           // NULL
        i++;
    }
}

int64_t strcmp_enc(vecCT v, char *s)
{
    vecCT v2;
    assign_string_to_vecCT(v2, s, -1);
    return strcmp_enc(v, v2);

    // int cmp_len=min(v.size(),strlen(s)); int i;
    // for(i=0;i<cmp_len;i++)
    // {
    //     int64_t x=operate_and_decrypt(v[i],"-",s[i]);
    //     if(x!=0)
    //         return x;
    // }
    // if(v.size()==strlen(s))
    //     return 0;
    // else if(v.size()>strlen(s))
    //     return decrypt_ciphertext_to_plaintext_vector(v[i])[0];
    // return -1*int64_t(s[i]);
}

int64_t strcmp_enc(vecCT v1, vecCT v2)
{
    // int cmp_len = min(v1.size(), v2.size());
    int cmp_len = min(strlen_enc(v1), strlen_enc(v2));

    int i;
    for (i = 0; i < cmp_len; i++)
    {
        int64_t x = operate_and_decrypt(v1[i], "-", v2[i]);     // x=0 if v1[i]=v2[i]
        if (x != 0)
            return x;
    }
    if (strlen_enc(v1)==strlen_enc(v2))
        return 0;
    else if (strlen_enc(v1)>strlen_enc(v2))
        return decrypt_ciphertext_to_plaintext_vector(v1[i])[0];

    return -1 * decrypt_ciphertext_to_plaintext_vector(v2[i])[0];
}

int64_t strcmp_enc(char *s, vecCT v)
{
    vecCT v2;
    assign_string_to_vecCT(v2, s, -1);
    return strcmp_enc(v2, v);
    // int cmp_len=min(v.size(),strlen(s)); int i;
    // for(i=0;i<cmp_len;i++)
    // {
    //     int64_t x=operate_and_decrypt(v[i],"-",s[i]);
    //     if(x!=0)
    //         return -1*x;
    // }
    // if(v.size()==strlen(s))
    //     return 0;
    // else if(v.size()>strlen(s))
    //     return -1*decrypt_ciphertext_to_plaintext_vector(v[i])[0];
    // return int64_t(s[i]);
}

// CT deserialize_ciphertext_from_file(string ct_name)
// {
//     CT ct;
//     if (Serial::DeserializeFromFile("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/serializer_folder/"+ct_name+".txt" , ct, SerType::BINARY) == false) {
//         std::cerr << "Could not read the ciphertext "+ct_name+".txt" << std::endl;
//         return encrypt_plaintext_integer_to_ciphertext(1) ;
//     }

//     return ct;
// }

// CT serialize_ciphertext_to_file(string ct_name)
// {
//     CT ct;
//     if (Serial::DeserializeFromFile("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/serializer_folder/"+ct_name+".txt" , ct, SerType::BINARY) == false) {
//         std::cerr << "Could not read the ciphertext "+ct_name+".txt" << std::endl;
//         return encrypt_plaintext_integer_to_ciphertext(1) ;
//     }

//     return ct;
// }

int64_t sa_ls_word_i(int64_t i)
{
    int64_t ct_count = i / 16384;
    int64_t index = i % 16384;

    CT ct;
    if (Serial::DeserializeFromFile("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/packed_cts/sa_ls_word/" + to_string(ct_count), ct, SerType::BINARY) == false)
    {
        std::cerr << "Could not read the ciphertext sa_ls_word/" + to_string(ct_count) << std::endl;
        return INT_MIN;
    }

    vector<int64_t> vv = decrypt_ciphertext_to_plaintext_vector(ct);

    return vv[index];
}

int64_t sa_ms_byte_i(int64_t i)
{
    int64_t ct_count = i / 16384;
    int64_t index = i % 16384;

    CT ct;
    if (Serial::DeserializeFromFile("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/packed_cts/sa_ms_byte/" + to_string(ct_count), ct, SerType::BINARY) == false)
    {
        std::cerr << "Could not read the ciphertext sa_ms_byte/" + to_string(ct_count) << std::endl;
        return INT_MIN;
    }

    vector<int64_t> vv = decrypt_ciphertext_to_plaintext_vector(ct);

    return vv[index];
}

void signal_handler(int sig)
{
    printf("Caught signal %d\n", sig);
    exit(1);
}

// want cp_occ_cp_count[i][j]
CT cp_occ_cp_count_i(int64_t i, int64_t j)
{
    int64_t k = 4 * i + j;
    int64_t ct_count = k / 16384;
    int64_t index = k % 16384;

    CT ct;
    if (Serial::DeserializeFromFile("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/packed_cts/cp_occ_cp_count_file/" + to_string(ct_count), ct, SerType::BINARY) == false)
    {
        std::cerr << "Could not read the ciphertext cp_occ_cp_count_file/" + to_string(ct_count) << std::endl;
        // return INT_MIN;
        signal(SIGSEGV, signal_handler);
    }

    // vector<int64_t> vv=decrypt_ciphertext_to_plaintext_vector(ct);

    // return vv[index];

    return ct;
}

// want cp_occ_one_hot_bwt_str[i][j]
CT cp_occ_one_hot_bwt_str_i(int64_t i, int64_t j)
{
    int64_t k = 4 * i + j; // kth number
    int64_t ct_count = k / 16384;
    int64_t index = k % 16384;

    CT ct;
    if (Serial::DeserializeFromFile("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/packed_cts/cp_occ_one_hot_bwt_str_file/" + to_string(ct_count), ct, SerType::BINARY) == false)
    {
        std::cerr << "Could not read the ciphertext cp_occ_cp_count_file/" + to_string(ct_count) << std::endl;
        signal(SIGSEGV, signal_handler);
    }

    // vector<int64_t> vv=decrypt_ciphertext_to_plaintext_vector(ct);

    // return vv[index];

    return ct;
}

// want cp_occ_one_hot_bwt_str[i][j]
string p_str_i(int64_t i)
{

    ifstream file("string_start_stop_p_str.txt");

    string line;
    for (int j = 1; j <= i; j++) {
        if (!getline(file, line)) {
            cout << "Error: file does not have " << i << " lines" << endl;
            return NULL;
        }
    }

    // Parse the line into 4 integers
    int start_ct_count, start_index_in_ct, fin_ct_count, fin_index_in_ct;
    stringstream ss(line);
    vector<int> nums;
    int num;
    while (ss >> num) {
        nums.push_back(num);
        if (ss.peek() == ' ') {
            ss.ignore();
        }
    }
    if (nums.size() != 4) {
        cout << "Error: line does not have 4 integers" << endl;
        return NULL;
    }
    start_ct_count = nums[0];
    start_index_in_ct = nums[1];
    fin_ct_count = nums[2];
    fin_index_in_ct = nums[3];

    string str="";

    CT c1, c2;

    if (Serial::DeserializeFromFile("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/packed_cts/p_str/" + to_string(start_ct_count), c1, SerType::BINARY) == false)
    {
        std::cerr << "Could not read the ciphertext p_str/" + to_string(start_ct_count) << std::endl;
        signal(SIGSEGV, signal_handler);
    }
    if (Serial::DeserializeFromFile("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/packed_cts/p_str/" + to_string(fin_ct_count), c2, SerType::BINARY) == false)
    {
        std::cerr << "Could not read the ciphertext p_str/" + to_string(fin_ct_count) << std::endl;
        signal(SIGSEGV, signal_handler);
    }

    vector<int64_t> v1=decrypt_ciphertext_to_plaintext_vector(c1);
    vector<int64_t> v2=decrypt_ciphertext_to_plaintext_vector(c2);

    if(start_ct_count==fin_ct_count)
    {
        for(int i=start_index_in_ct;i<=fin_index_in_ct;i++)
            str+=v1[i];
    }
    else{
        for(int i=start_index_in_ct;i<v1.size();i++)
            str+=v1[i];
        for(int i=0;i<fin_index_in_ct;i++)
            str+=v2[i];
    }
    cout<<"str: "<<str<<endl;
    // vector<int64_t> vv=decrypt_ciphertext_to_plaintext_vector(ct);

    // return vv[index];

    return str;
}

