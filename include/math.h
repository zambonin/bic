#ifndef MATH_H
#define MATH_H

#include "common.h"

uint32_t min(const uint32_t a, const uint32_t b);

int32_t max(const int32_t a, const int32_t b);

long double lg(const uintx u);

double asqrt(double x);

uintx compute_bin(const uint32_t n, const uint32_t k, const bic_ctx_t ctx);

uintx compute_bic_no_cache(const uint32_t n, const uint32_t k,
                           const uint32_t d);

uintx compute_bic_with_sums(const uint32_t n, const uint32_t k,
                            const uint32_t d, intx *partial_sums,
                            const bic_ctx_t ctx);

uintx compute_bic(const uint32_t n, const uint32_t k, const uint32_t d,
                  const bic_ctx_t ctx);

uint16_t compute_acc(uintx *rop, const uint32_t n, const uint32_t k,
                     const uint32_t d, const bic_ctx_t ctx);

uintx compute_dir(const uint32_t n, const uint32_t k, const uint32_t d,
                  const uint32_t l, const bic_ctx_t ctx);

double exp_part_sum(const uint32_t n, const uint32_t k, const uint32_t d,
                    const uint32_t j, const uint32_t l, const bic_ctx_t ctx);

double stddev_part_sum(const uint32_t n, const uint32_t k, const uint32_t d,
                       const uint32_t j, const uint32_t l, const bic_ctx_t ctx);

void calc_hset(uint16_t k, uint16_t L, uint16_t *H, size_t *h_size);

#endif
