#ifndef COLEX_H
#define COLEX_H

#include "common.h"

void colex_unrank(const bic_ctx_t *ctx, uint32_t *rop, const uint16_t n,
                  const uint16_t k, const uint16_t d, const uintx r);

void colex_unrank_part_sums(const bic_ctx_t *ctx, uint32_t *rop,
                            const uint16_t n, const uint16_t k,
                            const uint16_t d, uintx r);

void colex_unrank_acc_linear(const bic_ctx_t *ctx, uint32_t *rop,
                             const uint16_t n, const uint16_t k,
                             const uint16_t d, const uintx r);

void colex_unrank_acc_bisect(const bic_ctx_t *ctx, uint32_t *rop,
                             const uint16_t n, const uint16_t k,
                             const uint16_t d, const uintx r);

// algorithm 3 of 10.1007/s13389-021-00264-9
void colex_unrank_acc_direct(const bic_ctx_t *ctx, uint32_t *rop,
                             const uint16_t n, const uint16_t k,
                             const uint16_t d, const uintx r);

uintx colex_rank(const bic_ctx_t *ctx, const uint16_t n, const uint16_t k,
                 const uint16_t d, const uint32_t *comb);

#endif
