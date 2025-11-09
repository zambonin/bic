#ifndef MATH_H
#define MATH_H

#include "common.h"

uint32_t min(const uint32_t a, const uint32_t b);

int32_t max(const int32_t a, const int32_t b);

long double lg(const uintx u);

double asqrt(double x);

// §6.1 of 10.1007/978-3-642-14764-7_6
uintx compute_bin(const uint32_t n, const uint32_t k, const bic_ctx_t ctx);

uintx compute_bic_no_cache(const uint16_t n, const uint16_t k,
                           const uint16_t d);

// remark 60 of 10.24033/asens.136
uintx compute_bic_with_sums(const uint16_t n, const uint16_t k,
                            const uint16_t d, intx *partial_sums,
                            const bic_ctx_t ctx);

uintx compute_bic(const uint16_t n, const uint16_t k, const uint16_t d,
                  const bic_ctx_t ctx);

uint16_t compute_acc(uintx *rop, const uint16_t n, const uint16_t k,
                     const uint16_t d, const bic_ctx_t ctx);

// proposition 3 of 10.1007/s13389-021-00264-9
uintx compute_dir(const uint16_t n, const uint16_t k, const uint16_t d,
                  const uint16_t l, const bic_ctx_t ctx);

double exp_part_sum(const uint16_t n, const uint16_t k, const uint16_t d,
                    const uint16_t j, const uint16_t l, const bic_ctx_t ctx);

double stddev_part_sum(const uint16_t n, const uint16_t k, const uint16_t d,
                       const uint16_t j, const uint16_t l, const bic_ctx_t ctx);

#endif
