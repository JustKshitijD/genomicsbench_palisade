#include <vector>
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include "omp.h"
#include "host_kernel.h"
#include "common.h"
#include "minimap.h"
#include "mmpriv.h"
#include "kalloc.h"

#include <chrono>
using namespace std::chrono;

static const char LogTable256[256] = {
#define LT(n) n, n, n, n, n, n, n, n, n, n, n, n, n, n, n, n
	-1, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
	LT(4), LT(5), LT(5), LT(6), LT(6), LT(6), LT(6),
	LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7)
};

static inline int ilog2_32(uint32_t v)
{
	uint32_t t, tt;
	if ((tt = v>>16)) return (t = tt>>8) ? 24 + LogTable256[t] : 16 + LogTable256[tt];
	return (t = v>>8) ? 8 + LogTable256[t] : LogTable256[v];
}

static inline int ilog2_32_ct(CT v)
{
	uint32_t t, tt;
	int64_t v_plain=decrypt_ciphertext_to_plaintext_vector(v)[0];
	tt=decrypt_ciphertext_to_plaintext_vector(shift_encrypted_bit_vector_and_return_integer(get_encrypted_bits_vector(v_plain),-16))[0];

	if(tt!=0)
	{
		CT t_ct=shift_encrypted_bit_vector_and_return_integer(get_encrypted_bits_vector((int64_t)(tt)),-8);
		if(decrypt_ciphertext_to_plaintext_vector(t_ct)[0]!=0)
			return 24+LogTable256[decrypt_ciphertext_to_plaintext_vector(t_ct)[0]];
		return 16+LogTable256[tt];
	}

	t=decrypt_ciphertext_to_plaintext_vector(shift_encrypted_bit_vector_and_return_integer(get_encrypted_bits_vector(v_plain),-8))[0];
	if(t!=0)
	{
		return 8+LogTable256[t];
	}
	return LogTable256[v_plain];
	
	// if ((tt = v>>16)) return (t = tt>>8) ? 24 + LogTable256[t] : 16 + LogTable256[tt];
	// return (t = v>>8) ? 8 + LogTable256[t] : LogTable256[v];
}

const int BACKSEARCH = 65;
#define MM_SEED_SEG_SHIFT  48
#define MM_SEED_SEG_MASK   (0xffULL<<(MM_SEED_SEG_SHIFT))

void chain_dp(call_t* a, return_t* ret)
{

	// TODO: make sure this works when n has more than 32 bits
	int64_t i, j, st = 0;
	int is_cdna = 0;
    const float gap_scale = 1.0f;
    const int max_iter = 5000;
    const int max_skip = 25;
    int max_dist_x = a->max_dist_x, max_dist_y = a->max_dist_y, bw = a->bw;
	CT max_dist_x_ct=encrypt_plaintext_integer_to_ciphertext(a->max_dist_x); 
	CT max_dist_y_ct=encrypt_plaintext_integer_to_ciphertext(a->max_dist_y); 
	CT bw_ct=encrypt_plaintext_integer_to_ciphertext(a->bw);

	int scaling_factor_in_avg_qspan=10;
    float avg_qspan = a->avg_qspan;
	// printf("avg_qspan: %f\n",avg_qspan);
	CT avg_qspan_ct=encrypt_plaintext_integer_to_ciphertext((int64_t)(a->avg_qspan*pow(2,scaling_factor_in_avg_qspan)));

    int n_segs = a->n_segs; CT n_segs_ct = encrypt_plaintext_integer_to_ciphertext(a->n_segs);
    int64_t n = a->n;	CT n_ct = encrypt_plaintext_integer_to_ciphertext(a->n);
	// printf("n: %lld\n",n);

	ret->n = n;
	ret->scores.resize(n);	ret->scores_ct.resize(ceil(n*1.0/16384));
	ret->parents.resize(n);	ret->parents_ct.resize(ceil(n*1.0/16384));
    ret->targets.resize(n);	ret->targets_ct.resize(ceil(n*1.0/16384));
    ret->peak_scores.resize(n);	ret->peak_scores_ct.resize(ceil(n*1.0/16384));
	// printf("ceil(n*1.0/16384): %d\n",(int)(ceil(n*1.0/16384)));
	// printf("ret->scores_ct.size(): %d\n",(int)(ret->scores_ct.size()));

	for(int k=0;k<(int)(ceil(n*1.0/16384));k++){
		// printf("k: %d\n",k);
		vecInt v(16384,0);
		ret->scores_ct[k]=encrypt_plaintext_vector_to_ciphertext(v);
		ret->parents_ct[k]=encrypt_plaintext_vector_to_ciphertext(v);
		ret->targets_ct[k]=encrypt_plaintext_vector_to_ciphertext(v);
		ret->peak_scores_ct[k]=encrypt_plaintext_vector_to_ciphertext(v);
		// printf("Initialized\n");
	}

	// fill the score and backtrack arrays
	for (i = 0; i < n; ++i) {
		auto i_start = high_resolution_clock::now();

		// if(i<a->index_ct_ends)												// Can work with CTs here
		if(a->anchors[i].x<p/2 && a->anchors[i].y<p/2)
		{
			// printf("\ni1: %d\n",i);
			
			//  ri = a->anchors[i].x;									// call_t.anchors[i].x is in range 9.2*10^18
			CT ri_ct=encrypt_plaintext_integer_to_ciphertext(a->anchors[i].x);

			int64_t max_j = -1;
			// int32_t qi = (int32_t)a->anchors[i].y, q_span = a->anchors[i].y>>32&0xff; // NB: only 8 bits of span is used!!!
			// uint64_t y_shifted_plain=a->anchors[i].y>>32;
			// printf("y_shifted_plain: %llu\n",y_shifted_plain);
			// uint64_t y_shift_and=y_shifted&0xff;

			CT qi_ct=encrypt_plaintext_integer_to_ciphertext((int32_t)(a->anchors[i].y)); 
			// printf("qi_ct: %lld\n",decrypt_ciphertext_to_plaintext_vector(qi_ct)[0]);
			
			vecCT y_vec=get_encrypted_bits_vector((int64_t)(a->anchors[i].y));
			CT y_shifted=shift_encrypted_bit_vector_and_return_integer(y_vec,-1*32);
			// printf("y_shifted: %llu\n",decrypt_ciphertext_to_plaintext_vector(y_shifted)[0]);
			CT q_span_2_ct=do_logical_and_of_encryted_bit_vectors(get_encrypted_bits_vector(decrypt_ciphertext_to_plaintext_vector(y_shifted)[0]),get_encrypted_bits_vector(255));
			int q_span_2=decrypt_ciphertext_to_plaintext_vector(q_span_2_ct)[0];
			// int q_span_2=decrypt_ciphertext_to_plaintext_vector(y_shifted)[0]&0xff;
			// printf("q_span_2: %lld\n",q_span_2);

			// int32_t max_f = q_span, n_skip = 0, min_d;
			int32_t max_f = q_span_2, n_skip = 0, min_d;
			// CT max_f_ct=encrypt_plaintext_integer_to_ciphertext(max_f);
			CT max_f_ct=q_span_2_ct;

			// int32_t sidi = (a->anchors[i].y & MM_SEED_SEG_MASK) >> MM_SEED_SEG_SHIFT;
			// int sidi_1=(a->anchors[i].y & MM_SEED_SEG_MASK) ;
			// printf("sidi_1: %d\n",sidi_1);
			CT sidi_1_ct=do_logical_and_of_encryted_bit_vectors(get_encrypted_bits_vector((int64_t)(a->anchors[i].y)),get_encrypted_bits_vector((int64_t)(MM_SEED_SEG_MASK)));
			// printf("sidi_1_ct: %lld\n",decrypt_ciphertext_to_plaintext_vector(sidi_1_ct)[0]);
			CT sidi_ct=shift_encrypted_bit_vector_and_return_integer(get_encrypted_bits_vector(decrypt_ciphertext_to_plaintext_vector(sidi_1_ct)[0]),-1*MM_SEED_SEG_SHIFT);
			// printf("sidi_ct: %lld\n",decrypt_ciphertext_to_plaintext_vector(sidi_ct)[0]);

			// while (st < i && ri > a->anchors[st].x + max_dist_x) ++st;					// Adding, comparing for ri
			while (st < i && decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(ri_ct,cc->EvalAdd(encrypt_plaintext_integer_to_ciphertext(a->anchors[st].x),max_dist_x_ct)))[0]>0) ++st;					// Adding, comparing for ri
			// printf("st: %lld\n",st);

			if (i - st > max_iter) st = i - max_iter;

			for (j = i - 1; j >= st; --j) {

				// printf("ri_ct: %lld\n",decrypt_ciphertext_to_plaintext_vector(ri_ct)[0]);
				// printf("(int64_t)(a->anchors[j].x): %lld, (int32_t)(a->anchors[j].y): %lld\n",(int64_t)(a->anchors[j].x),(int32_t)(a->anchors[j].y));

				// printf("j: %d\n",j);
				// int64_t dr = ri - a->anchors[j].x;
				CT dr_ct = cc->EvalSub(ri_ct,encrypt_plaintext_integer_to_ciphertext((int64_t)(a->anchors[j].x)));

				// int32_t dq = qi - (int32_t)a->anchors[j].y, 
				// int dd, sc, log_dd, gap_cost;
				CT dq_ct=cc->EvalSub(qi_ct,encrypt_plaintext_integer_to_ciphertext((int32_t)(a->anchors[j].y)));

				// printf("dr_ct: %lld, dq_ct: %lld\n",decrypt_ciphertext_to_plaintext_vector(dr_ct)[0],decrypt_ciphertext_to_plaintext_vector(dq_ct)[0]);

				// int32_t sidj = (a->anchors[j].y & MM_SEED_SEG_MASK) >> MM_SEED_SEG_SHIFT;
				CT sidj_ct=shift_encrypted_bit_vector_and_return_integer(get_encrypted_bits_vector(decrypt_ciphertext_to_plaintext_vector(do_logical_and_of_encryted_bit_vectors(get_encrypted_bits_vector((int64_t)(a->anchors[j].y)),get_encrypted_bits_vector((int64_t)(MM_SEED_SEG_MASK))))[0]),-1*MM_SEED_SEG_SHIFT);
				// printf("sidj: %d, sidj_ct: %lld\n",sidj,decrypt_ciphertext_to_plaintext_vector(sidj_ct)[0]);

				// printf("sidi: %lld, sidj: %lld, dr: %d, dq: %lld, max_dist_x: %lld, max_dist_y: %lld\n",decrypt_ciphertext_to_plaintext_vector(sidi_ct)[0],decrypt_ciphertext_to_plaintext_vector(sidj_ct)[0],decrypt_ciphertext_to_plaintext_vector(dr_ct)[0],decrypt_ciphertext_to_plaintext_vector(dq_ct)[0],max_dist_x,max_dist_y);
				// printf("decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(sidi_ct,sidj_ct))[0]: %d\n",decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(sidi_ct,sidj_ct))[0]);

				// if ((sidi == sidj && dr == 0) || dq <= 0) continue; // don't skip if an anchor is used by multiple segments; see below
				if (((int)(decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(sidi_ct,sidj_ct))[0])==0 && (int)(decrypt_ciphertext_to_plaintext_vector(dr_ct)[0]) == 0) || (int)(decrypt_ciphertext_to_plaintext_vector(dq_ct)[0]) <= 0) continue; // don't skip if an anchor is used by multiple segments; see below
				// printf("1\n");
				// printf("decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(dq_ct,encode_integer_to_plaintext(max_dist_y)))[0]: %lld\n",decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(dq_ct,encode_integer_to_plaintext(max_dist_y)))[0]);
				// printf("decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(dq_ct,encode_integer_to_plaintext(max_dist_x)))[0]: %lld\n",decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(dq_ct,encode_integer_to_plaintext(max_dist_x)))[0]);
				int y2=decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(dq_ct,encode_integer_to_plaintext(max_dist_y)))[0];
				int y3=decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(dq_ct,encode_integer_to_plaintext(max_dist_x)))[0];
				// printf("y2: %lld, y3: %lld\n",y2,y3);
				
				bool b1= ((int)(decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(sidi_ct,sidj_ct))[0])==0);
				bool b2= (y2 > 0);
				bool b3= (y3 > 0);

				// printf("b1: %d, b2: %d, b3: %d\n",b1,b2,b3);

				// if ((sidi == sidj && dq > max_dist_y) || dq > max_dist_x) continue;
				if ((b1 && b2) || b3) continue; // don't skip if an anchor is used by multiple segments; see below
				// printf("2\n");

				// dd = dr > dq? dr - dq : dq - dr;
				CT dd_ct = (int)(decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(dr_ct,dq_ct))[0])>0?cc->EvalSub(dr_ct,dq_ct):cc->EvalSub(dq_ct,dr_ct);
				// printf("dd: %d, dd_ct: %lld\n",dd,decrypt_ciphertext_to_plaintext_vector(dd_ct)[0]);
				// printf("dd: %lld, bw: %lld, bw_ct: %lld\n",decrypt_ciphertext_to_plaintext_vector(dd_ct)[0],bw,decrypt_ciphertext_to_plaintext_vector(bw_ct)[0]);
				// printf("decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(dd_ct,bw_ct))[0]: %d\n",decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(dd_ct,bw_ct))[0]);
				int x2=(int)(decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(dd_ct,bw_ct))[0]);
				// printf("x2: %d\n",x2);

				b1=((int)(decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(sidi_ct,sidj_ct))[0])==0);
				b2=(x2>0);
				// printf("bool (decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(sidi_ct,sidj_ct))[0]==0): %d, x2>0: %d\n",b1,b2);

				// if (sidi == sidj && dd > bw) continue;
				if(b1 && b2) continue;
				// rintf("3\n");

				b1=((int)(decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(n_segs_ct,encode_integer_to_plaintext(1))))[0]>0);
				b2=((int)(decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(sidi_ct,sidj_ct))[0]==0));
				b3=((int)(decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(dr_ct,encode_integer_to_plaintext(max_dist_y)))[0])>0);
				// printf("bool (decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(n_segs_ct,encode_integer_to_plaintext(1)))[0]>0): %d, (decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(sidi_ct,sidj_ct))[0]==0): %d, (decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(dr_ct,encode_integer_to_plaintext(max_dist_y)))[0]>0): %d\n",b1,b2,b3);

				// if (n_segs > 1 && !is_cdna && sidi == sidj && dr > max_dist_y) continue;
				if(b1 && !is_cdna && b2 && b3) continue;

				// min_d = dq < dr? dq : dr;
				b1=(int)(decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(dq_ct,dr_ct))[0])<0;
				CT min_d_ct = b1? dq_ct:dr_ct;

				// printf("min_d_ct: %lld\n",decrypt_ciphertext_to_plaintext_vector(min_d_ct)[0]);

				// sc = min_d > q_span? q_span : dq < dr? dq : dr;
				CT sc_ct = (int)(decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(min_d_ct,q_span_2_ct))[0]>0)?q_span_2_ct:(int)(decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(dq_ct,dr_ct))[0])<0?dq_ct:dr_ct;
				// printf("sc_ct1: %lld\n",decrypt_ciphertext_to_plaintext_vector(sc_ct)[0]);

				// log_dd = dd? ilog2_32(dd) : 0; 
				// printf("dd: %lld\n",(decrypt_ciphertext_to_plaintext_vector(dd_ct)[0]));
				CT log_dd_ct; b1=((int)(decrypt_ciphertext_to_plaintext_vector(dd_ct)[0])!=0); 
				// printf("bool decrypt_ciphertext_to_plaintext_vector(dd_ct)[0]!=0: %d\n",b1);
				// printf("dd: %lld\n",decrypt_ciphertext_to_plaintext_vector(dd_ct)[0]);
				if(b1)
				{
					log_dd_ct=encrypt_plaintext_integer_to_ciphertext(ilog2_32_ct(dd_ct));
				}
				else
				{
					// printf("Inside dd: %lld\n",decrypt_ciphertext_to_plaintext_vector(dd_ct)[0]);
					log_dd_ct=encrypt_plaintext_integer_to_ciphertext(0);
				}
				// printf("dd: %lld\n",decrypt_ciphertext_to_plaintext_vector(dd_ct)[0]);
				// printf("log_dd_ct: %lld\n",decrypt_ciphertext_to_plaintext_vector(log_dd_ct)[0]);

				// gap_cost = 0; 
				CT gap_cost_ct=encrypt_plaintext_integer_to_ciphertext(0);
				// printf("is_cdna: %d\n",is_cdna);

				b1=((int)(decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(sidi_ct,sidj_ct))[0])!=0);
				// printf("bool decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(sidi_ct,sidj_ct))[0]!=0: %d\n",b1);
				// printf("dd: %lld\n",decrypt_ciphertext_to_plaintext_vector(dd_ct)[0]);

				if (is_cdna || b1) {
					int c_log, c_lin; CT c_log_ct, c_lin_ct;
					// c_lin = (int)(dd * .01 * avg_qspan);
					CT g1=cc->EvalMult(dd_ct,encode_integer_to_plaintext((int64_t)(0.01*pow(2,scaling_factor_in_avg_qspan))));
					c_lin_ct=cc->EvalMult(g1,avg_qspan_ct);	// This is scaled up by 2^20
					// printf("Scaled c_lin_ct: %lld\n",decrypt_ciphertext_to_plaintext_vector(c_lin_ct)[0]);

					// c_log = log_dd;
					c_log_ct=log_dd_ct;
					CT c_lin_ct_2=shift_encrypted_bit_vector_and_return_integer(get_encrypted_bits_vector(decrypt_ciphertext_to_plaintext_vector(c_lin_ct)[0]),-20);
					
					// printf("Descaled c_lin_ct: %lld\n",decrypt_ciphertext_to_plaintext_vector(c_lin_ct)[0]);
					
					// If decrypt(c_lin_ct_2)==(int)(dd*0.01*avg_qspan), then, our scaling and descaling was successful

					// printf("c_lin: %d, c_lin_ct_2: %lld, c_log: %d, c_log_ct: %lld\n",c_lin,decrypt_ciphertext_to_plaintext_vector(c_lin_ct_2)[0],c_log,decrypt_ciphertext_to_plaintext_vector(c_log_ct)[0]);
					b1=((int)(decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(sidi_ct,sidj_ct))[0])!=0);
					b2=((int)(decrypt_ciphertext_to_plaintext_vector(dr_ct)[0])==0);
					// printf("(decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(sidi_ct,sidj_ct))[0]!=0): %d\n",b1);
					// printf("decrypt_ciphertext_to_plaintext_vector(dr_ct)[0]: %d\n",b2);
					b3=((int)(decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(dr_ct,dq_ct))[0])>0);
					bool b4=((int)(decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(sidi_ct,sidj_ct))[0])!=0);
					// printf("else if (decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(dr_ct,dq_ct))[0]>0): %d\n",b3);
					// printf("else if (decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(sidi_ct,sidj_ct))[0]!=0): %d\n",b4);

					if(b1 && b2)
					{
						//sc++;
						sc_ct=cc->EvalAdd(sc_ct,encode_integer_to_plaintext(1));
						// printf("Inner sc2: %lld\n",decrypt_ciphertext_to_plaintext_vector(sc_ct)[0]);
						
						//printf("sc: %d, sc_ct: %lld\n",sc,decrypt_ciphertext_to_plaintext_vector(sc_ct)[0]);
					}
					else if(b3 || b4)
					{
						if(decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(c_lin_ct_2,c_log_ct))[0]<0)
							gap_cost_ct=c_lin_ct_2;
						else
							gap_cost_ct=c_log_ct;

						// printf("gap_cost2: %d\n",decrypt_ciphertext_to_plaintext_vector(gap_cost_ct)[0]);
						// gap_cost = c_lin < c_log? c_lin : c_log;
						// printf("else if: gap_cost: %d, gap_cost_ct: %lld\n",gap_cost,decrypt_ciphertext_to_plaintext_vector(gap_cost_ct)[0]);
					}
					else
					{
						gap_cost_ct=cc->EvalAdd(c_lin_ct_2,shift_encrypted_bit_vector_and_return_integer(get_encrypted_bits_vector(decrypt_ciphertext_to_plaintext_vector(c_log_ct)[0]),-1));
						// printf("gap_cost3: %d\n",decrypt_ciphertext_to_plaintext_vector(gap_cost_ct)[0]);
						// gap_cost = c_lin + (c_log>>1);
						// printf("else1: gap_cost: %d, gap_cost_ct: %lld\n",gap_cost,decrypt_ciphertext_to_plaintext_vector(gap_cost_ct)[0]);
					}
					// c_lin = (int)(dd * .01 * avg_qspan);
					// c_log = log_dd;
					// if (sidi != sidj && dr == 0) ++sc; // possibly due to overlapping paired ends; give a minor bonus
					// else if (dr > dq || sidi != sidj) gap_cost = c_lin < c_log? c_lin : c_log;
					// else gap_cost = c_lin + (c_log>>1);
				}
				else{
					// printf("dd: %lld\n",decrypt_ciphertext_to_plaintext_vector(dd_ct)[0]);
					CT pow_scaled_ct=encrypt_plaintext_integer_to_ciphertext((int64_t)(0.01*pow(2,scaling_factor_in_avg_qspan)));
					// printf("pow_scaled_ct: %lld\n",decrypt_ciphertext_to_plaintext_vector(pow_scaled_ct)[0]);
					CT g1=cc->EvalMult(dd_ct,pow_scaled_ct);
					// printf("avg_qspan_ct: %lld\n",decrypt_ciphertext_to_plaintext_vector(avg_qspan_ct)[0]);
					// printf("g1: %lld\n",decrypt_ciphertext_to_plaintext_vector(g1)[0]);
					CT g2=cc->EvalMult(g1,avg_qspan_ct);	// This is scaled up by 2^20
					// printf("g2: %lld\n",decrypt_ciphertext_to_plaintext_vector(g2)[0]);
					CT g3=shift_encrypted_bit_vector_and_return_integer(get_encrypted_bits_vector(decrypt_ciphertext_to_plaintext_vector(g2)[0]),-20);
					// printf("g3: %lld\n",decrypt_ciphertext_to_plaintext_vector(g3)[0]);
					// If decrypt(c_lin_ct_2)==(int)(dd*0.01*avg_qspan), then, our scaling and descaling was successful

					CT shift_log_dd=shift_encrypted_bit_vector_and_return_integer(get_encrypted_bits_vector(decrypt_ciphertext_to_plaintext_vector(log_dd_ct)[0]),-1);
					// printf("shift_log_dd: %lld\n",decrypt_ciphertext_to_plaintext_vector(shift_log_dd)[0]);
					gap_cost_ct=cc->EvalAdd(g3,shift_log_dd);
					// gap_cost = (int)(dd * .01 * avg_qspan) + (log_dd>>1);
					// printf("gap_cost4: %d\n",decrypt_ciphertext_to_plaintext_vector(gap_cost_ct)[0]);

					// printf("else 2:- (int)(dd * .01 * avg_qspan): %d, g3: %lld\n",(int)(dd * .01 * avg_qspan), decrypt_ciphertext_to_plaintext_vector(g3)[0]);
					// printf("gap_cost: %d, gap_cost_ct: %lld\n",gap_cost,decrypt_ciphertext_to_plaintext_vector(gap_cost_ct)[0]);
				}

				// sc -= (int)((double)gap_cost * gap_scale + .499);
				// int k1=(int)((double)gap_cost * gap_scale + .499);
				CT k2=cc->EvalMult(gap_cost_ct,encode_integer_to_plaintext((int64_t)((int64_t)(gap_scale)*pow(2,10))));
				CT k1_ct=cc->EvalAdd(k2,encode_integer_to_plaintext((int64_t)(0.499*pow(2,scaling_factor_in_avg_qspan))));
				CT k1_ct_shifted=shift_encrypted_bit_vector_and_return_integer(get_encrypted_bits_vector(decrypt_ciphertext_to_plaintext_vector(k1_ct)[0]),-20);
				// printf("k1: %d, k1_ct: %lld\n",k1,decrypt_ciphertext_to_plaintext_vector(k1_ct)[0]);
				int k1_shifted=(int)(decrypt_ciphertext_to_plaintext_vector(k1_ct_shifted)[0]);
				// printf("sub: %lld\n",k1_shifted);

				sc_ct=cc->EvalSub(sc_ct,encode_integer_to_plaintext(k1_shifted));
				// printf("sc2: %lld\n",decrypt_ciphertext_to_plaintext_vector(sc_ct)[0]);
				//printf("sc: %d, sc_ct: %lld\n",sc,decrypt_ciphertext_to_plaintext_vector(sc_ct)[0]);

				// sc += ret->scores[j];
				sc_ct=cc->EvalAdd(sc_ct,encode_integer_to_plaintext(decrypt_ciphertext_to_plaintext_vector(ret->scores_ct[j/16384])[j%16384]));
				// printf("sc3: %lld\n",decrypt_ciphertext_to_plaintext_vector(sc_ct)[0]);
				// printf("sc: %d, sc_ct: %lld\n",sc,decrypt_ciphertext_to_plaintext_vector(sc_ct)[0]);

				if(decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(sc_ct,max_f_ct))[0]>0){
					max_f_ct=sc_ct; max_j=j; max_f=decrypt_ciphertext_to_plaintext_vector(sc_ct)[0];
					if (n_skip > 0) --n_skip;
				}else if(decrypt_ciphertext_to_plaintext_vector(ret->targets_ct[j/16384])[j%16384]==i){
					if (++n_skip > max_skip) {
						break;
					}
				}

				// printf("ret->parents[%d]: %d\n",j,decrypt_ciphertext_to_plaintext_vector(ret->parents_ct[j/16384])[j%16384]);

				if(decrypt_ciphertext_to_plaintext_vector(ret->parents_ct[j/16384])[j%16384]>=0)
				{
					// ret->targets_ct[decrypt_ciphertext_to_plaintext_vector(ret->parents_ct[j])[0]]=encrypt_plaintext_integer_to_ciphertext(i);
					ret->targets[ret->parents[j]] = i;

					int p_index=decrypt_ciphertext_to_plaintext_vector(ret->parents_ct[j/16384])[j%16384];
					vecInt vv=decrypt_ciphertext_to_plaintext_vector(ret->targets_ct[p_index/16384]);
					vv[p_index%16384]=i;
					ret->targets_ct[p_index/16384]=encrypt_plaintext_vector_to_ciphertext(vv);

					// printf("ret->targets[%d]: %d\n",p_index,decrypt_ciphertext_to_plaintext_vector(ret->targets_ct[p_index/16384])[p_index%16384]); // ret->parents[j]=p_index
				}
				// if (sc > max_f) {
				// 	max_f = sc, max_j = j;
				// 	if (n_skip > 0) --n_skip;
				// } else if (ret->targets[j] == i) {
				// 	if (++n_skip > max_skip) {
				// 		break;
				// 	}
				// }
				// if (ret->parents[j] >= 0) ret->targets[ret->parents[j]] = i;
			}
			ret->scores[i] = max_f, ret->parents[i] = max_j;
			ret->peak_scores[i] = max_j >= 0 && ret->peak_scores[max_j] > max_f ? ret->peak_scores[max_j] : max_f;
			
			vecInt vv=decrypt_ciphertext_to_plaintext_vector(ret->scores_ct[i/16384]);	
			vv[i%16384]=max_f;
			ret->scores_ct[i/16384]=encrypt_plaintext_vector_to_ciphertext(vv);

			vv=decrypt_ciphertext_to_plaintext_vector(ret->parents_ct[i/16384]);	
			vv[i%16384]=max_j;
			ret->parents_ct[i/16384]=encrypt_plaintext_vector_to_ciphertext(vv);

			vv=decrypt_ciphertext_to_plaintext_vector(ret->peak_scores_ct[i/16384]);
			vv[i%16384]=max_j >= 0 && decrypt_ciphertext_to_plaintext_vector(ret->peak_scores_ct[max_j/16384])[max_j%16384] > max_f ? decrypt_ciphertext_to_plaintext_vector(ret->peak_scores_ct[max_j/16384])[max_j%16384]: max_f;
			ret->peak_scores_ct[i/16384]=encrypt_plaintext_vector_to_ciphertext(vv);

			printf("ret->scores[%d]: %d, %lld\n",i,ret->scores[i],decrypt_ciphertext_to_plaintext_vector(ret->scores_ct[i/16384])[i%16384]);
			printf("ret->parents[%d]: %d, %lld\n",i,ret->parents[i],decrypt_ciphertext_to_plaintext_vector(ret->parents_ct[i/16384])[i%16384]);
			printf("ret->peak_scores[%d]: %d, %lld\n",i,ret->peak_scores[i],decrypt_ciphertext_to_plaintext_vector(ret->peak_scores_ct[i/16384])[i%16384]);

			// ret->scores_ct[i] = max_f_ct, ret->parents_ct[i] = encrypt_plaintext_integer_to_ciphertext(max_j);
			// ret->peak_scores_ct[i] = max_j >= 0 && decrypt_ciphertext_to_plaintext_vector(cc->EvalSub(ret->peak_scores_ct[max_j], max_f_ct))[0]>0 ? ret->peak_scores_ct[max_j] : max_f_ct;
		
			// printf("ret->scores[%d]: %d, ret->scores_ct[%d]: %lld\n",i,ret->scores[i],i,decrypt_ciphertext_to_plaintext_vector(ret->scores_ct[i])[0]);
			// printf("ret->parents[%d]: %d, ret->parents_ct[%d]: %lld\n",i,ret->parents[i],i,decrypt_ciphertext_to_plaintext_vector(ret->parents_ct[i])[0]);
			// printf("ret->peak_scores[%d]: %d, ret->peak_scores_ct[%d]: %lld\n",i,ret->peak_scores[i],i,decrypt_ciphertext_to_plaintext_vector(ret->peak_scores_ct[i])[0]);
		
		}
		else{
			printf("\ni2: %d\n",i);
			uint64_t ri = a->anchors[i].x;									// call_t.anchors[i].x is in range 9.2*10^18
			int64_t max_j = -1;
			int32_t qi = (int32_t)a->anchors[i].y, q_span = a->anchors[i].y>>32&0xff; // NB: only 8 bits of span is used!!!
			int32_t max_f = q_span, n_skip = 0, min_d;
			int32_t sidi = (a->anchors[i].y & MM_SEED_SEG_MASK) >> MM_SEED_SEG_SHIFT;
			while (st < i && ri > a->anchors[st].x + max_dist_x) ++st;					// Adding, comparing for ri
			if (i - st > max_iter) st = i - max_iter;
			for (j = i - 1; j >= st; --j) {
				int64_t dr = ri - a->anchors[j].x;
				int32_t dq = qi - (int32_t)a->anchors[j].y, dd, sc, log_dd, gap_cost;
				int32_t sidj = (a->anchors[j].y & MM_SEED_SEG_MASK) >> MM_SEED_SEG_SHIFT;
				if ((sidi == sidj && dr == 0) || dq <= 0) continue; // don't skip if an anchor is used by multiple segments; see below
				if ((sidi == sidj && dq > max_dist_y) || dq > max_dist_x) continue;
				dd = dr > dq? dr - dq : dq - dr;
				if (sidi == sidj && dd > bw) continue;
				if (n_segs > 1 && !is_cdna && sidi == sidj && dr > max_dist_y) continue;
				min_d = dq < dr? dq : dr;
				sc = min_d > q_span? q_span : dq < dr? dq : dr;
				log_dd = dd? ilog2_32(dd) : 0;
				gap_cost = 0;
				if (is_cdna || sidi != sidj) {
					int c_log, c_lin;
					c_lin = (int)(dd * .01 * avg_qspan);
					c_log = log_dd;
					if (sidi != sidj && dr == 0) ++sc; // possibly due to overlapping paired ends; give a minor bonus
					else if (dr > dq || sidi != sidj) gap_cost = c_lin < c_log? c_lin : c_log;
					else gap_cost = c_lin + (c_log>>1);
				} else gap_cost = (int)(dd * .01 * avg_qspan) + (log_dd>>1);
				sc -= (int)((double)gap_cost * gap_scale + .499);
				sc += ret->scores[j];
				if (sc > max_f) {
					max_f = sc, max_j = j;
					if (n_skip > 0) --n_skip;
				} else if (ret->targets[j] == i) {
					if (++n_skip > max_skip) {
						break;
					}
				}
				if (ret->parents[j] >= 0) {
					ret->targets[ret->parents[j]] = i;

					int p_index=decrypt_ciphertext_to_plaintext_vector(ret->parents_ct[j/16384])[j%16384];
					vecInt vv=decrypt_ciphertext_to_plaintext_vector(ret->targets_ct[p_index/16384]);
					vv[p_index%16384]=i;
					ret->targets_ct[p_index/16384]=encrypt_plaintext_vector_to_ciphertext(vv);
				}
			}
			ret->scores[i] = max_f, ret->parents[i] = max_j;
			ret->peak_scores[i] = max_j >= 0 && ret->peak_scores[max_j] > max_f ? ret->peak_scores[max_j] : max_f;

			vecInt vv=decrypt_ciphertext_to_plaintext_vector(ret->scores_ct[i/16384]);	
			vv[i%16384]=max_f;
			ret->scores_ct[i/16384]=encrypt_plaintext_vector_to_ciphertext(vv);

			vv=decrypt_ciphertext_to_plaintext_vector(ret->parents_ct[i/16384]);	
			vv[i%16384]=max_j;
			ret->parents_ct[i/16384]=encrypt_plaintext_vector_to_ciphertext(vv);

			vv=decrypt_ciphertext_to_plaintext_vector(ret->peak_scores_ct[i/16384]);
			vv[i%16384]=max_j >= 0 && decrypt_ciphertext_to_plaintext_vector(ret->peak_scores_ct[max_j/16384])[max_j%16384] > max_f ? decrypt_ciphertext_to_plaintext_vector(ret->peak_scores_ct[max_j/16384])[max_j%16384]: max_f;
			ret->peak_scores_ct[i/16384]=encrypt_plaintext_vector_to_ciphertext(vv);

			printf("ret->scores[%d]: %d, %lld\n",i,ret->scores[i],decrypt_ciphertext_to_plaintext_vector(ret->scores_ct[i/16384])[i%16384]);
			printf("ret->parents[%d]: %d, %lld\n",i,ret->parents[i],decrypt_ciphertext_to_plaintext_vector(ret->parents_ct[i/16384])[i%16384]);
			printf("ret->peak_scores[%d]: %d, %lld\n",i,ret->peak_scores[i],decrypt_ciphertext_to_plaintext_vector(ret->peak_scores_ct[i/16384])[i%16384]);
		}

		auto i_end = high_resolution_clock::now();

		auto duration = duration_cast<microseconds>(i_end - i_start);
		printf("Time taken by i: %d is %lld microseconds\n",i,duration.count());
	}
}

void host_chain_kernel(std::vector<call_t> &args, std::vector<return_t> &rets, int numThreads)
{
	printf("args.size(): %d\n",(int)(args.size()));
    #pragma omp parallel num_threads(numThreads)
    {
        #pragma omp for schedule(dynamic)
            for (size_t batch = 0; batch < args.size(); batch++) {
				printf("-----------------------------------------------batch: %d\n",batch);
                call_t* arg = &args[batch];
                return_t* ret = &rets[batch];
                // fprintf(stderr, "%lld\t%f\t%d\t%d\t%d\t%d\n", arg->n, arg->avg_qspan, arg->max_dist_x, arg->max_dist_y, arg->bw, arg->n_segs);

				auto batch_start = high_resolution_clock::now();
				chain_dp(arg, ret);
				auto batch_end = high_resolution_clock::now();

				auto duration = duration_cast<microseconds>(batch_end - batch_start);
				printf("Time taken by batch %d is: %lld microseconds\n",(int)(batch),duration.count());
            }
    }
}
