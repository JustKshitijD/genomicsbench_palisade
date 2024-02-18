#ifndef HOST_INPUT_H
#define HOST_INPUT_H

#include <vector>
#include <cstdint>
#include "../../../palisade_header.h"

typedef int64_t anchor_idx_t;
typedef uint32_t tag_t;
typedef int32_t loc_t;
typedef int32_t loc_dist_t;
typedef int32_t score_t;
typedef int32_t parent_t;
typedef int32_t target_t;
typedef int32_t peak_score_t;

#define ANCHOR_NULL (anchor_idx_t)(-1)

struct anchor_t {
    uint64_t x;
    uint64_t y;
};

struct call_t {
    anchor_idx_t n; CT n_ct;
    float avg_qspan; CT avg_qspan_ct; int scaling_factor_in_avg_qspan=3;
    int max_dist_x, max_dist_y, bw, n_segs;
    CT max_dist_x_ct, max_dist_y_ct, bw_ct, n_segs_ct;

    std::vector<anchor_t> anchors;
    int index_ct_ends;    // Of the n seeds, the seeds before this index are represented in anchors_x, anchors_y; the seeds after this index
                            // are represented by x,y in vector<anchor_t> anchors, because these x,y values are very big

    CT other_vars;
    vecCT anchors_x; vecCT anchors_y;
};

struct return_t {
    anchor_idx_t n;
    std::vector<score_t> scores;
    std::vector<parent_t> parents;
    std::vector<target_t> targets;
    std::vector<peak_score_t> peak_scores;

    vecCT scores_ct, parents_ct, targets_ct,peak_scores_ct;
};

#endif // HOST_INPUT_H
