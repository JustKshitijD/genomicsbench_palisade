/*!
 * @file alignment_engine.cpp
 *
 * @brief AlignmentEngine class source file
 */

#include <limits>
#include <algorithm>
#include <exception>
#include <stdexcept>

#include <iostream>

#include "sisd_alignment_engine.hpp"
#include "simd_alignment_engine.hpp"
#include "spoa/alignment_engine.hpp"
#include "palisade_header.h"

namespace spoa {

std::unique_ptr<AlignmentEngine> createAlignmentEngine(AlignmentType type,
    std::int8_t m, std::int8_t n, std::int8_t g) {

    return createAlignmentEngine(type, m, n, g, g);
}

std::unique_ptr<AlignmentEngine> createAlignmentEngine(AlignmentType type,
    std::int8_t m, std::int8_t n, std::int8_t g, std::int8_t e) {

    return createAlignmentEngine(type, m, encrypt_plaintext_integer_to_ciphertext(m), n, encrypt_plaintext_integer_to_ciphertext(n), g, encrypt_plaintext_integer_to_ciphertext(g), e, encrypt_plaintext_integer_to_ciphertext(e), g, encrypt_plaintext_integer_to_ciphertext(g), e, encrypt_plaintext_integer_to_ciphertext(e));
}

std::unique_ptr<AlignmentEngine> createAlignmentEngine(AlignmentType type,
    std::int8_t m, CT enc_m, std::int8_t n, CT enc_n, std::int8_t g, CT enc_g, std::int8_t e, CT enc_e, 
    std::int8_t q, CT enc_q, std::int8_t c, CT enc_c) {

        printf("In createAlignmentEngine, m: %d, n: %d, g: %d, e: %d, q: %d, c: %d\n",m,n,g,e,q,c);

    if (type != AlignmentType::kSW &&
        type != AlignmentType::kNW &&
        type != AlignmentType::kOV) {

        throw std::invalid_argument("[spoa::createAlignmentEngine] error: "
            "invalid alignment type!");
    }
    if (operate_and_decrypt(enc_g,"-",0) > 0 || operate_and_decrypt(enc_q,"-",0) > 0 ) {
        throw std::invalid_argument("[spoa::createAlignmentEngine] error: "
            "gap opening penalty must be non-positive!");
    }
    if (operate_and_decrypt(enc_e,"-",0) > 0  || operate_and_decrypt(enc_c,"-",0) > 0 ) {
        throw std::invalid_argument("[spoa::createAlignmentEngine] error: "
            "gap extension penalty must be non-positive!");
    }

    AlignmentSubtype subtype = operate_and_decrypt(enc_g,"-",e) < 0  ?
       ( operate_and_decrypt(enc_g,"-",q) <= 0 || operate_and_decrypt(enc_e,"-",c) >= 0 ?
        AlignmentSubtype::kAffine : AlignmentSubtype::kConvex) : AlignmentSubtype::kLinear;

    if (subtype == AlignmentSubtype::kLinear) {
        e = g;
        enc_e=enc_g;
    } else if (subtype == AlignmentSubtype::kAffine) {
        q = g;
        enc_q=enc_g;

        c = e;
        enc_c=enc_e;
    }

    auto alignment_engine = createSimdAlignmentEngine(type, subtype, m, n, g, e,
        q, c);

    if (alignment_engine == nullptr) {
        return createSisdAlignmentEngine(type, subtype, m, n, g, e, q, c);
    }

    return alignment_engine;
}

AlignmentEngine::AlignmentEngine(AlignmentType type, AlignmentSubtype subtype,
    std::int8_t m, std::int8_t n, std::int8_t g, std::int8_t e,
    std::int8_t q, std::int8_t c)
        : type_(type), subtype_(subtype), m_(m), n_(n), g_(g), e_(e), q_(q), c_(c) {
}

// using Alignment = std::vector<std::pair<std::int32_t, std::int32_t>>;
Alignment AlignmentEngine::align(const std::string& sequence,
    const std::unique_ptr<Graph>& graph) {
    // std::cout<<"In align, with sequence: "<<sequence<<std::endl;
    // std::cout<<"sequence.size(): "<<sequence.size()<<std::endl;

    // Goes to- 
    // Breakpoint 1, 0x0000555555664450 in spoa::SimdAlignmentEngine<(spoa::Arch)3>::align(char const*, unsigned int, std::unique_ptr<spoa::Graph, std::default_delete<spoa::Graph> > const&) ()
    return align(sequence.c_str(), sequence.size(), graph);
}

}


