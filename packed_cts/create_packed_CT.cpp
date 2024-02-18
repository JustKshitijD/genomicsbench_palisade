#include "palisade_header.h"
#include <fstream>
#include <string.h>
#include <iostream>
#include <vector>

// int main()
// {
//     init();
//     ifstream inputFile;

//     inputFile.open("../cp_occ_one_hot_bwt_str_file.txt");
//     int num;

//     long long k=0; vector<int64_t> vv;

//     if (inputFile.is_open()) {

//         while (inputFile >> num) {
//             if(k%16384==0)
//             {
//                 cout<<"k: "<<k<<endl;
//                 vv.clear();
//                 if(k!=0)
//                 {
//                     CT c=encrypt_plaintext_vector_to_ciphertext(vv);
//                     // Serialize cryptocontext
//                     if (!Serial::SerializeToFile("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/packed_cts/cp_occ_one_hot_bwt_str_file/"+to_string(k),c,SerType::BINARY)) {
//                         std::cerr << "Error writing serialization of the crypto context to "
//                                     "cryptocontext.txt"
//                                 << std::endl;
//                     }
//                     // std::cout << "The cryptocontext has been serialized." << std::endl;
//                 }
//             }
//             else{
//                 vv.push_back(num);
//             }

//             k++;
//         }

//         inputFile.close();
//     }
//     else {
//         cout << "Error opening file.";
//     }

//     return 0;
// } 

int main()
{
    init();
    ifstream inputFile;

    printf("p_gi_2\n");
    inputFile.open("../p_gi_2.txt");
    int num;

    long long k=0; vector<int64_t> vv;

    // if (inputFile.is_open()) {

    //     while (inputFile >> num) {

    //         cout<<"num: "<<num<<endl;
    //         cout<<"k: "<<k<<endl<<endl;

    //         if(k%16384==0)
    //         {
    //             vv.clear();
    //             if(k!=0)
    //             {
    //                 CT c=encrypt_plaintext_vector_to_ciphertext(vv);
    //                 // Serialize cryptocontext
    //                 if (!Serial::SerializeToFile("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/packed_cts/cp_occ_one_hot_bwt_str_file/"+to_string(k),c,SerType::BINARY)) {
    //                     std::cerr << "Error writing serialization of the crypto context to "
    //                                 "cryptocontext.txt"
    //                             << std::endl;
    //                 }
    //                 // std::cout << "The cryptocontext has been serialized." << std::endl;
    //             }
    //         }
    //         else{
    //             vv.push_back(num);
    //         }

    //         k++;
    //     }

    //     inputFile.close();
    // }
    // else {
    //     cout << "Error opening file.";
    // }

    // inputFile.open("../cp_occ_cp_count_file.txt");

    // if (inputFile.is_open()) {

    //     while (inputFile >> num) {

    //         cout<<"num: "<<num<<endl;
    //         cout<<"k: "<<k<<endl<<endl;

    //         if(k%16384==0)
    //         {
    //             vv.clear();
    //             if(k!=0)
    //             {
    //                 CT c=encrypt_plaintext_vector_to_ciphertext(vv);
    //                 // Serialize cryptocontext
    //                 if (!Serial::SerializeToFile("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/packed_cts/cp_occ_cp_count_file/"+to_string(k),c,SerType::BINARY)) {
    //                     std::cerr << "Error writing serialization of the crypto context to "
    //                                 "cryptocontext.txt"
    //                             << std::endl;
    //                 }
    //                 // std::cout << "The cryptocontext has been serialized." << std::endl;
    //             }
    //         }
    //         else{
    //             vv.push_back(num);
    //         }

    //         k++;
    //     }

    //     inputFile.close();
    // }
    // else {
    //     cout << "Error opening file.";
    // }

    // inputFile.open("../sa_ls_word.txt");

    // if (inputFile.is_open()) {

    //     while (inputFile >> num) {
    //         if(k%16384==0)
    //         {
    //             cout<<"k: "<<k<<endl;
    //             vv.clear();
    //             if(k!=0)
    //             {
    //                 CT c=encrypt_plaintext_vector_to_ciphertext(vv);
    //                 // Serialize cryptocontext
    //                 if (!Serial::SerializeToFile("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/packed_cts/sa_ls_word/"+to_string(k),c,SerType::BINARY)) {
    //                     std::cerr << "Error writing serialization of the crypto context to "
    //                                 "cryptocontext.txt"
    //                             << std::endl;
    //                 }
    //                 // std::cout << "The cryptocontext has been serialized." << std::endl;
    //             }
    //         }
    //         else{
    //             vv.push_back(num);
    //         }

    //         k++;
    //     }

    //     inputFile.close();
    // }
    // else {
    //     cout << "Error opening file.";
    // }

    // inputFile.open("../sa_ms_byte.txt");

    // if (inputFile.is_open()) {

    //     while (inputFile >> num) {
    //         if(k%16384==0)
    //         {
    //             cout<<"k: "<<k<<endl;
    //             vv.clear();
    //             if(k!=0)
    //             {
    //                 CT c=encrypt_plaintext_vector_to_ciphertext(vv);
    //                 // Serialize cryptocontext
    //                 if (!Serial::SerializeToFile("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/packed_cts/sa_ms_byte/"+to_string(k),c,SerType::BINARY)) {
    //                     std::cerr << "Error writing serialization of the crypto context to "
    //                                 "cryptocontext.txt"
    //                             << std::endl;
    //                 }
    //                 // std::cout << "The cryptocontext has been serialized." << std::endl;
    //             }
    //         }
    //         else{
    //             vv.push_back(num);
    //         }

    //         k++;
    //     }

    //     inputFile.close();
    // }
    // else {
    //     cout << "Error opening file.";
    // }

    if (inputFile.is_open()) {

        while (inputFile >> num) {

            cout<<"num: "<<num<<endl;
            cout<<"k: "<<k<<endl<<endl;

            if(k%16384==0)
            {
                vv.clear();
                if(k!=0)
                {
                    CT c=encrypt_plaintext_vector_to_ciphertext(vv);
                    // Serialize cryptocontext
                    if (!Serial::SerializeToFile("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/packed_cts/p_gi/"+to_string(k),c,SerType::BINARY)) {
                        std::cerr << "Error writing serialization of the crypto context to "
                                    "cryptocontext.txt"
                                << std::endl;
                    }
                    // std::cout << "The cryptocontext has been serialized." << std::endl;
                }
            }
            else{
                vv.push_back(num);
            }

            k++;
        }

        if(vv.size()!=0)
        {
            CT c=encrypt_plaintext_vector_to_ciphertext(vv);
            // Serialize cryptocontext
            if (!Serial::SerializeToFile("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/packed_cts/p_gi/"+to_string(k),c,SerType::BINARY)) {
                std::cerr << "Error writing serialization of the crypto context to "
                            "cryptocontext.txt"
                        << std::endl;
            }
        }

        inputFile.close();
    }
    else {
        cout << "Error opening file.";
    }

    k=0;

    printf("--------------------------------------\n");
    printf("p_len_2\n");
    inputFile.open("../p_len_2.txt");

    if (inputFile.is_open()) {

        while (inputFile >> num) {

            cout<<"num: "<<num<<endl;
            cout<<"k: "<<k<<endl<<endl;

            if(k%16384==0)
            {
                vv.clear();
                if(k!=0)
                {
                    CT c=encrypt_plaintext_vector_to_ciphertext(vv);
                    // Serialize cryptocontext
                    if (!Serial::SerializeToFile("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/packed_cts/p_len/"+to_string(k),c,SerType::BINARY)) {
                        std::cerr << "Error writing serialization of the crypto context to "
                                    "cryptocontext.txt"
                                << std::endl;
                    }
                    // std::cout << "The cryptocontext has been serialized." << std::endl;
                }
            }
            else{
                vv.push_back(num);
            }

            k++;
        }

        if(vv.size()!=0)
        {
            CT c=encrypt_plaintext_vector_to_ciphertext(vv);
            // Serialize cryptocontext
            if (!Serial::SerializeToFile("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/packed_cts/p_len/"+to_string(k),c,SerType::BINARY)) {
                std::cerr << "Error writing serialization of the crypto context to "
                            "cryptocontext.txt"
                        << std::endl;
            }
        }

        inputFile.close();
    }
    else {
        cout << "Error opening file.";
    }

    k=0;

    printf("--------------------------------------\n");
    printf("p_anno_2\n");
    ifstream infile2("../p_anno_2.txt");

    vector<char> buffer; int ct_count=0; buffer.clear(); string line;
    vector<pair<pair<int,int>,pair<int,int> > > string_start_stop_p_anno; pair<int,int> st; pair<int,int> en;
    // a vector for string- ((ciphertext_count_storing_first_char_of_string,index_in_the_vector_which_that_ciphertext_is_encrypting),(ciphertext_count_storing_last_char_of_string,index_in_the_vector_which_that_ciphertext_is_encrypting))
    while (getline(infile2, line)) {
        // Add each character to the buffer
       int ind=0;
       cout<<"line in anno: "<<line<<endl;

        for (char c : line) {
            buffer.push_back(c);

            if(ind==0){
                st=make_pair(ct_count,buffer.size()-1);
            }
            else if(ind==line.size()-1){
                en=make_pair(ct_count,buffer.size()-1);
                string_start_stop_p_anno.push_back(make_pair(st,en));
            }

            // If the buffer is full, create a new vector and start filling it
            if (buffer.size() == 16384) {
                vector<char> chunk(buffer.begin(), buffer.end());
                vector<int64_t> chunk_2;
                for(int k=0;k<chunk.size();k++)
                    chunk_2.push_back(int64_t(chunk[k]));

                printf("chunk_2\n");
                for(int k=0;k<chunk_2.size();k++)
                    printf("%c",chunk_2[k]);
                printf("\n");
        
                CT cc=encrypt_plaintext_vector_to_ciphertext(chunk_2);
                if (!Serial::SerializeToFile("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/packed_cts/p_anno/"+to_string(ct_count),cc,SerType::BINARY)) {
                    std::cerr << "Error writing serialization of the crypto context to "
                                "cryptocontext.txt"
                            << std::endl;
                }

                ct_count++;

                // Do something with the chunk, like store it in a vector of vectors
                buffer.clear();
            }

            ind++;
        }
    }

    // If there are any remaining characters in the buffer, create a final chunk
    if (!buffer.empty()) {
        vector<char> chunk(buffer.begin(), buffer.end());
        vector<int64_t> chunk_2;
        for(int k=0;k<chunk.size();k++)
            chunk_2.push_back(int64_t(chunk[k]));

        printf("chunk_2\n");
        for(int k=0;k<chunk_2.size();k++)
            printf("%c",chunk_2[k]);
        printf("\n");

        CT cc=encrypt_plaintext_vector_to_ciphertext(chunk_2);
        if (!Serial::SerializeToFile("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/packed_cts/p_anno/"+to_string(ct_count),cc,SerType::BINARY)) {
            std::cerr << "Error writing serialization of the crypto context to "
                        "cryptocontext.txt"
                    << std::endl;
        }

        en=make_pair(ct_count,buffer.size()-1);

        ct_count++;

        // Do something with the chunk, like store it in a vector of vectors
        buffer.clear();

        string_start_stop_p_anno.push_back(make_pair(st,en));
        // Do something with the final chunk
    }

    ofstream outFile("string_start_stop_p_anno.txt");
    printf("string_start_stop_p_anno:- \n");
    for(int i=0;i<string_start_stop_p_anno.size();i++){
        printf("i: %d:- ",i);
        pair<pair<int,int>,pair<int,int> >p=string_start_stop_p_anno[i];
        printf("((%d,%d),(%d,%d))\n",(p.first).first,(p.first).second,(p.second).first,(p.second).second);
        outFile << (p.first).first << " " << (p.first).second << " " << (p.second).first << " " << (p.second).second << endl;
    }
    printf("\n");

    // if (inputFile.is_open()) {
    //     printf("Anno is open!\n");
    //     while (inputFile >> num) {

    //         cout<<"num: "<<num<<endl;
    //         cout<<"k: "<<k<<endl<<endl;

    //         if(k%16384==0)
    //         {
    //             vv.clear();
    //             if(k!=0)
    //             {
    //                 CT c=encrypt_plaintext_vector_to_ciphertext(vv);
    //                 // Serialize cryptocontext
    //                 if (!Serial::SerializeToFile("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/packed_cts/p_anno/"+to_string(k),c,SerType::BINARY)) {
    //                     std::cerr << "Error writing serialization of the crypto context to "
    //                                 "cryptocontext.txt"
    //                             << std::endl;
    //                 }
    //                 // std::cout << "The cryptocontext has been serialized." << std::endl;
    //             }
    //         }
    //         else{
    //             vv.push_back(num);
    //         }

    //         k++;
    //     }

    //     if(vv.size()!=0)
    //     {
    //         CT c=encrypt_plaintext_vector_to_ciphertext(vv);
    //         // Serialize cryptocontext
    //         if (!Serial::SerializeToFile("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/packed_cts/p_anno/"+to_string(k),c,SerType::BINARY)) {
    //             std::cerr << "Error writing serialization of the crypto context to "
    //                         "cryptocontext.txt"
    //                     << std::endl;
    //         }
    //     }

    //     inputFile.close();
    // }
    // else {
    //     cout << "Error opening file.";
    // }

    k=0;

    printf("--------------------------------------\n");
    printf("p_len_2\n");

    inputFile.open("../p_len_2.txt");

    if (inputFile.is_open()) {

        while (inputFile >> num) {

            cout<<"num: "<<num<<endl;
            cout<<"k: "<<k<<endl<<endl;

            if(k%16384==0)
            {
                vv.clear();
                if(k!=0)
                {
                    CT c=encrypt_plaintext_vector_to_ciphertext(vv);
                    // Serialize cryptocontext
                    if (!Serial::SerializeToFile("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/packed_cts/p_len/"+to_string(k),c,SerType::BINARY)) {
                        std::cerr << "Error writing serialization of the crypto context to "
                                    "cryptocontext.txt"
                                << std::endl;
                    }
                    // std::cout << "The cryptocontext has been serialized." << std::endl;
                }
            }
            else{
                vv.push_back(num);
            }

            k++;
        }

        if(vv.size()!=0)
        {
            CT c=encrypt_plaintext_vector_to_ciphertext(vv);
            // Serialize cryptocontext
            if (!Serial::SerializeToFile("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/packed_cts/p_len/"+to_string(k),c,SerType::BINARY)) {
                std::cerr << "Error writing serialization of the crypto context to "
                            "cryptocontext.txt"
                        << std::endl;
            }
        }

        inputFile.close();
    }
    else {
        cout << "Error opening file.";
    }

    k=0;

    printf("--------------------------------------\n");
    printf("p_n_ambs_2\n");

    inputFile.open("../p_n_ambs_2.txt");

    if (inputFile.is_open()) {

        while (inputFile >> num) {

            cout<<"num: "<<num<<endl;
            cout<<"k: "<<k<<endl<<endl;

            if(k%16384==0)
            {
                vv.clear();
                if(k!=0)
                {
                    CT c=encrypt_plaintext_vector_to_ciphertext(vv);
                    // Serialize cryptocontext
                    if (!Serial::SerializeToFile("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/packed_cts/p_n_ambs/"+to_string(k),c,SerType::BINARY)) {
                        std::cerr << "Error writing serialization of the crypto context to "
                                    "cryptocontext.txt"
                                << std::endl;
                    }
                    // std::cout << "The cryptocontext has been serialized." << std::endl;
                }
            }
            else{
                vv.push_back(num);
            }

            k++;
        }

        if(vv.size()!=0)
        {
            CT c=encrypt_plaintext_vector_to_ciphertext(vv);
            // Serialize cryptocontext
            if (!Serial::SerializeToFile("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/packed_cts/p_n_ambs/"+to_string(k),c,SerType::BINARY)) {
                std::cerr << "Error writing serialization of the crypto context to "
                            "cryptocontext.txt"
                        << std::endl;
            }
        }

        inputFile.close();
    }
    else {
        cout << "Error opening file.";
    }

    k=0;

    printf("--------------------------------------\n");
    printf("p_offset_2\n");

    inputFile.open("../p_offset_2.txt");

    if (inputFile.is_open()) {
        printf("Offset is open!\n"); int64_t num2;

        while (inputFile >> num2) {

            cout<<"num2: "<<num2<<endl;
            cout<<"k: "<<k<<endl<<endl;

            if(k%16384==0)
            {
                vv.clear();
                if(k!=0)
                {
                    CT c=encrypt_plaintext_vector_to_ciphertext(vv);
                    // Serialize cryptocontext
                    if (!Serial::SerializeToFile("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/packed_cts/p_offset/"+to_string(k),c,SerType::BINARY)) {
                        std::cerr << "Error writing serialization of the crypto context to "
                                    "cryptocontext.txt"
                                << std::endl;
                    }
                    // std::cout << "The cryptocontext has been serialized." << std::endl;
                }
            }
            else{
                vv.push_back(num2);
            }

            k++;
        }

        if(vv.size()!=0)
        {
            CT c=encrypt_plaintext_vector_to_ciphertext(vv);
            // Serialize cryptocontext
            if (!Serial::SerializeToFile("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/packed_cts/p_offset/"+to_string(k),c,SerType::BINARY)) {
                std::cerr << "Error writing serialization of the crypto context to "
                            "cryptocontext.txt"
                        << std::endl;
            }
        }

        inputFile.close();
    }
    else {
        cout << "Error opening file.";
    }

    k=0;

    printf("--------------------------------------\n");
    printf("p_str_2\n");

    // inputFile.open("../p_str_2.txt");
    ifstream infile("../p_str_2.txt");

    // Read strings from the file
    // string line;
    buffer.clear(); ct_count=0;
    vector<pair<pair<int,int>,pair<int,int> > > string_start_stop_p_str;  
    // a vector for string- ((ciphertext_count_storing_first_char_of_string,index_in_the_vector_which_that_ciphertext_is_encrypting),(ciphertext_count_storing_last_char_of_string,index_in_the_vector_which_that_ciphertext_is_encrypting))
    
    // pair<int,int> st; pair<int,int> en;
    while (getline(infile, line)) {
        cout<<"line in str: "<<line<<endl;
        // Add each character to the buffer
       int ind=0; 
        for (char c : line) {
            buffer.push_back(c);

            if(ind==0){
                st=make_pair(ct_count,buffer.size()-1);
            }
            else if(ind==line.size()-1){
                en=make_pair(ct_count,buffer.size()-1);
                string_start_stop_p_str.push_back(make_pair(st,en));
            }

            // If the buffer is full, create a new vector and start filling it
            if (buffer.size() == 16384) {
                vector<char> chunk(buffer.begin(), buffer.end());
                vector<int64_t> chunk_2;
                for(int k=0;k<chunk.size();k++)
                    chunk_2.push_back(int64_t(chunk[k]));

                printf("chunk_2\n");
                for(int k=0;k<chunk_2.size();k++)
                    printf("%c",chunk_2[k]);
                printf("\n");
        
                CT cc=encrypt_plaintext_vector_to_ciphertext(chunk_2);
                if (!Serial::SerializeToFile("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/packed_cts/p_str/"+to_string(ct_count),cc,SerType::BINARY)) {
                    std::cerr << "Error writing serialization of the crypto context to "
                                "cryptocontext.txt"
                            << std::endl;
                }

                ct_count++;

                // Do something with the chunk, like store it in a vector of vectors
                buffer.clear();
            }

            ind++;
        }
    }

    // If there are any remaining characters in the buffer, create a final chunk
    if (!buffer.empty()) {
        vector<char> chunk(buffer.begin(), buffer.end());
        vector<int64_t> chunk_2;
        for(int k=0;k<chunk.size();k++)
            chunk_2.push_back(int64_t(chunk[k]));

        printf("chunk_2\n");
        for(int k=0;k<chunk_2.size();k++)
            printf("%c",chunk_2[k]);
        printf("\n");

        CT cc=encrypt_plaintext_vector_to_ciphertext(chunk_2);
        if (!Serial::SerializeToFile("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/packed_cts/p_str/"+to_string(ct_count),cc,SerType::BINARY)) {
            std::cerr << "Error writing serialization of the crypto context to "
                        "cryptocontext.txt"
                    << std::endl;
        }

        en=make_pair(ct_count,buffer.size()-1);

        ct_count++;

        // Do something with the chunk, like store it in a vector of vectors
        buffer.clear();

        string_start_stop_p_str.push_back(make_pair(st,en));
        // Do something with the final chunk
    }

    ofstream outFile2("string_start_stop_p_str.txt");
    printf("string_start_stop_p_str:- \n");
    for(int i=0;i<string_start_stop_p_str.size();i++){
        printf("i: %d:- ",i);
        pair<pair<int,int>,pair<int,int> >p=string_start_stop_p_str[i];
        printf("((%d,%d),(%d,%d))\n",(p.first).first,(p.first).second,(p.second).first,(p.second).second);
        outFile2<<(p.first).first<<" "<<(p.first).second<<" "<<(p.second).first<<" "<<(p.second).second<<endl;
    }
    printf("\n");

    // Close the file stream
    infile.close();

    // if(vv.size()!=0)
    // {
    //     CT c=encrypt_plaintext_vector_to_ciphertext(vv);
    //     // Serialize cryptocontext
    //     if (!Serial::SerializeToFile("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/packed_cts/p_str/"+to_string(k),c,SerType::BINARY)) {
    //         std::cerr << "Error writing serialization of the crypto context to "
    //                     "cryptocontext.txt"
    //                 << std::endl;
    //     }
    // }

    // if (inputFile.is_open()) {

    //     while (inputFile >> num) {

    //         cout<<"num: "<<num<<endl;
    //         cout<<"k: "<<k<<endl<<endl;

    //         if(k%16384==0)
    //         {
    //             vv.clear();
    //             if(k!=0)
    //             {
    //                 CT c=encrypt_plaintext_vector_to_ciphertext(vv);
    //                 // Serialize cryptocontext
    //                 if (!Serial::SerializeToFile("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/packed_cts/p_str/"+to_string(k),c,SerType::BINARY)) {
    //                     std::cerr << "Error writing serialization of the crypto context to "
    //                                 "cryptocontext.txt"
    //                             << std::endl;
    //                 }
    //                 // std::cout << "The cryptocontext has been serialized." << std::endl;
    //             }
    //         }
    //         else{
    //             vv.push_back(num);
    //         }

    //         k++;
    //     }

    //     if(vv.size()!=0)
    //     {
    //         CT c=encrypt_plaintext_vector_to_ciphertext(vv);
    //         // Serialize cryptocontext
    //         if (!Serial::SerializeToFile("/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/packed_cts/p_str/"+to_string(k),c,SerType::BINARY)) {
    //             std::cerr << "Error writing serialization of the crypto context to "
    //                         "cryptocontext.txt"
    //                     << std::endl;
    //         }
    //     }

    //     inputFile.close();
    // }
    // else {
    //     cout << "Error opening file.";
    // }

    return 0;
} 


