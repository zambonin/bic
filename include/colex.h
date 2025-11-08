#ifndef COLEX_H
#define COLEX_H

#include "common.h"

void colex_unrank(uint32_t *rop, const uint16_t n, const uint16_t k,
                  const uint16_t d, const uintx r, const bic_ctx_t ctx);

void colex_unrank_part_sums(uint32_t *rop, const uint16_t n, const uint16_t k,
                            const uint16_t d, const uintx r,
                            const bic_ctx_t ctx);

void colex_unrank_acc_linear(uint32_t *rop, const uint16_t n, const uint16_t k,
                             const uint16_t d, const uintx r,
                             const bic_ctx_t ctx);

void colex_unrank_acc_bisect(uint32_t *rop, const uint16_t n, const uint16_t k,
                             const uint16_t d, const uintx r,
                             const bic_ctx_t ctx);

// algorithm 3 of 10.1007/s13389-021-00264-9
void colex_unrank_acc_direct(uint32_t *rop, const uint16_t n, const uint16_t k,
                             const uint16_t d, const uintx r,
                             const bic_ctx_t ctx);

uintx colex_rank(const uint16_t n, const uint16_t k, const uint16_t d,
                 const uint32_t *comb, const bic_ctx_t ctx);

#endif
