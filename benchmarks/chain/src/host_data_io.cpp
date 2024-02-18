#include "host_data_io.h"
#include "host_data.h"
#include <iostream>
#include <iomanip>

void skip_to_EOR(FILE *fp) {
    const char *loc = "EOR";
    while (*loc != '\0') {
        if (fgetc(fp) == *loc) {
            loc++;
        }
    }
}

/*
Input sequence-
In every call to read_call, follwing are read- 

long long n;
float avg_qspan;
int max_dist_x (5000 always), max_dist_y(5000 always), bw(500 always), n_segs(1 always)

Read call.n (x,y) values
*/

// typedef int64_t anchor_idx_t;

// struct anchor_t {
//     uint64_t x;
//     uint64_t y;
// };

// struct call_t {
//     anchor_idx_t n;
//     float avg_qspan;
//     int max_dist_x, max_dist_y, bw, n_segs;
//     std::vector<anchor_t> anchors;
// };

call_t read_call(FILE *fp) {
    call_t call;

    long long n;
    float avg_qspan;
    int max_dist_x, max_dist_y, bw, n_segs;

    int t = fscanf(fp, "%lld%f%d%d%d%d",
            &n, &avg_qspan, &max_dist_x, &max_dist_y, &bw, &n_segs);
    // fprintf(stderr, "read %d arguments\n", t);
    if (t != 6) {
        call.n = ANCHOR_NULL;
        call.avg_qspan = .0;
        return call;
    }

    // printf("n: %lld\n",n);
    
    call.n = n;
    call.avg_qspan = avg_qspan;
    call.max_dist_x = max_dist_x;
    call.max_dist_y = max_dist_y;
    call.bw = bw;
    call.n_segs = n_segs;
    // fprintf(stderr, "%lld\t%f\t%d\t%d\t%d\t%d\n", n, avg_qspan, max_dist_x, max_dist_y, bw, n_segs);

    call.anchors.resize(call.n);

    for (anchor_idx_t i = 0; i < call.n; i++) {
        uint64_t x, y;
        fscanf(fp, "%llu%llu", &x, &y);

        anchor_t t;
        t.x = x; t.y = y;

        call.anchors[i] = t;
    }

    skip_to_EOR(fp);
    return call;
}

// call_t read_call(FILE *fp, int &max_n) {
//     call_t call;

//     long long n;
//     float avg_qspan;
//     int max_dist_x, max_dist_y, bw, n_segs;

//     int t = fscanf(fp, "%lld%f%d%d%d%d",
//             &n, &avg_qspan, &max_dist_x, &max_dist_y, &bw, &n_segs);
//     // fprintf(stderr, "read %d arguments\n", t);
//     if (t != 6) {
//         call.n = ANCHOR_NULL;
//         call.avg_qspan = .0;
//         return call;
//     }

//     call.n = n; call.n_ct=encrypt_plaintext_integer_to_ciphertext(n);
//     call.avg_qspan = avg_qspan; call.avg_qspan_ct=encrypt_plaintext_integer_to_ciphertext(avg_qspan);
//     call.max_dist_x = max_dist_x;   call.max_dist_x_ct=encrypt_plaintext_integer_to_ciphertext(max_dist_x);
//     call.max_dist_y = max_dist_y;   call.max_dist_y_ct=encrypt_plaintext_integer_to_ciphertext(max_dist_y);
//     call.bw = bw;   call.bw_ct=encrypt_plaintext_integer_to_ciphertext(bw);
//     call.n_segs = n_segs;   call.n_segs_ct=encrypt_plaintext_integer_to_ciphertext(n_segs);

//     call.index_ct_ends=-1;
//     // fprintf(stderr, "%lld\t%f\t%d\t%d\t%d\t%d\n", n, avg_qspan, max_dist_x, max_dist_y, bw, n_segs);

//     // std::cout<<"avg_qspan: "<<avg_qspan<<std::endl;
//     // std::cout << std::fixed << std::setprecision(15) << "avg_qspan: "<< avg_qspan << std::endl;

//     string x_str=to_string(avg_qspan);                      // can be 31.251556396484375
//     int ind=x_str.find('.');
    
//     string int_x=x_str.substr(0,ind);                       // 31   
//     string float_x=x_str.substr(ind+1,x_str.size()-ind-1);  // 251556
//     string new_num_x=int_x+float_x;
//     size_t sz;
//     int64_t new_x_val = std::stoll(new_num_x,&sz,0);          // so, new_x_val is 31251556; we can recover original float number from this because int_x.size() is always 2. So, old float number is 31.251556
//     // printf("new_x_val: %lld\n\n",new_x_val);

//     // size_t sz;
//     // int64_t int_x_val = std::stoll (int_x,&sz,0);
//     // printf("%lld\n",int_x_val);
    
//     // int64_t float_x_val = std::stoll (float_x,&sz,0);
//     // printf("%lld\n\n",float_x_val);

//     if(n>max_n)
//         max_n=n;

//     // if((int)(int_x.length())>max_n)
//     //     max_n=(int)(int_x.length());

//     vecInt other_vars_v; 
//     other_vars_v.push_back(n); other_vars_v.push_back(new_x_val); other_vars_v.push_back((int64_t)(max_dist_x)); other_vars_v.push_back((int64_t)(max_dist_y));
//     other_vars_v.push_back((int64_t)(bw)); other_vars_v.push_back((int64_t)(n_segs));
//     call.other_vars=encrypt_plaintext_vector_to_ciphertext(other_vars_v);

//     call.anchors.resize(call.n);

//     for (anchor_idx_t i = 0; i < call.n; i++) {
//         uint64_t x, y;
//         fscanf(fp, "%llu%llu", &x, &y);         // x values increase cumulatively; y fluctuate randomly

//         // printf("original i: %d, x: %llu, y: %llu\n",i,x,y);

//         if((x>p/2 || y>p/2) && call.index_ct_ends==-1){
//             call.index_ct_ends=i;                       // At index_ct_ends, the first x>p value lies
//         }

//         anchor_t t;
//         t.x = x; t.y = y;

//         call.anchors[i] = t;
//     }

//     cout<<"call.index_ct_ends: "<<call.index_ct_ends<<endl;         // can be 0!
//     int anchor_ct_count;  int anchor_ct_index;

//     // if(call.index_ct_ends==0)
//     // {
//     //     anchor_ct_count=1; anchor_ct_index=0;
//     // }
//     if(call.index_ct_ends!=0){
//         anchor_ct_count=ceil(call.index_ct_ends*1.0/16384); anchor_ct_index=(call.index_ct_ends)%16384;
//     }

//     if(call.index_ct_ends!=0)
//     {
//         // printf("anchor_ct_count: %d\n",anchor_ct_count);
//         call.anchors_x.resize(anchor_ct_count); call.anchors_y.resize(anchor_ct_count);
//     }

//     vecInt vv_x(16384,0); vecInt vv_y(16384,0); int ct_x=0, ct_y=0; int ct_count=0;

//     for (anchor_idx_t i =0; i <call.index_ct_ends ; i++) {
//         // cout<<"encrypting i: "<<i<<endl;

//         if(i!=0 && i%16384==0){
//             printf("encrypting i= %d\n",i);

//             call.anchors_x[ct_count]=encrypt_plaintext_vector_to_ciphertext(vv_x);
//             call.anchors_y[ct_count]=encrypt_plaintext_vector_to_ciphertext(vv_y);

//             fill(vv_x.begin(),vv_x.end(),0);
//             fill(vv_y.begin(),vv_y.end(),0);
//             ct_count++;
//         }

//         uint64_t x = call.anchors[i].x;
//         uint64_t y = call.anchors[i].y;
//         // cout<<"i: "<<i<<"; x: "<<x<<"; y: "<<y<<endl;

//         vv_x[i-16384*ct_count]=x;
//         vv_y[i-16384*ct_count]=y;

//     }
    
//     if(call.index_ct_ends!=0)
//     {
//         call.anchors_x[ct_count]=encrypt_plaintext_vector_to_ciphertext(vv_x);
//         call.anchors_y[ct_count]=encrypt_plaintext_vector_to_ciphertext(vv_y);
//         fill(vv_x.begin(),vv_x.end(),0);
//         fill(vv_y.begin(),vv_y.end(),0);
//     }

//     int ct_number=0; int ct_index=0;

//     for(int i=0;i<call.n;i++){
//         // cout<<"----------------------------------"<<endl;
//         // printf("i: %d\n",i);
//         if(i<call.index_ct_ends){
//             // cout<<"----------------------------------"<<endl;
//             // printf("i: %d\n",i);

//             if(i%16384==0){
//                 ct_number=i/16384;
//                 vv_x=decrypt_ciphertext_to_plaintext_vector(call.anchors_x[ct_number]);
//                 vv_y=decrypt_ciphertext_to_plaintext_vector(call.anchors_y[ct_number]);
//             }

//             // printf("call.anchors[i].x: %llu, vv_x[%d]: %llu\n",call.anchors[i].x,i%16384,vv_x[i%16384]);
//             // printf("call.anchors[i].y: %llu, vv_y[%d]: %llu\n",call.anchors[i].y,i%16384,vv_y[i%16384]);
//             // printf("\n");
//         }
//         // else{
//         //     if (i==call.index_ct_ends)
//         //     {
//         //         printf("###################################\n");
//         //         cout<<"index ct ends!!! "<<i<<endl<<endl;
//         //     }
//         //     // printf("call.anchors[i].x: %llu, call.anchors[i].y: %llu\n",call.anchors[i].x,call.anchors[i].y);
//         // }
//     }

//     skip_to_EOR(fp);
//     return call;
// }

void print_return(FILE *fp, const return_t &data)
{
    fprintf(fp, "%lld\n", (long long)data.n);
    for (anchor_idx_t i = 0; i < data.n; i++) {
        fprintf(fp, "%d\t%d\n", (int)data.scores[i], (int)data.parents[i]);
    }
    fprintf(fp, "EOR\n");
}
