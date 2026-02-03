#ifndef COLEX_H
#define COLEX_H

#include "common.h"
#include "generic.h"

void colex_unrank(uint32_t *rop, const uint32_t n, const uint32_t k,
                  const uint32_t d, const uintx r, const bic_ctx_t ctx);

void colex_unrank_part_sums(uint32_t *rop, const uint32_t n, const uint32_t k,
                            const uint32_t d, const uintx r,
                            const bic_ctx_t ctx);

void colex_unrank_acc_linear(uint32_t *rop, const uint32_t n, const uint32_t k,
                             const uint32_t d, const uintx r,
                             const bic_ctx_t ctx);

void colex_unrank_acc_bisect(uint32_t *rop, const uint32_t n, const uint32_t k,
                             const uint32_t d, const uintx r,
                             const bic_ctx_t ctx);

void colex_unrank_acc_direct(uint32_t *rop, const uint32_t n, const uint32_t k,
                             const uint32_t d, const uintx r,
                             const bic_ctx_t ctx);

uintx colex_rank(const uint32_t n, const uint32_t k, const uint32_t d,
                 const uint32_t *comb, const bic_ctx_t ctx);

uintx colex_rank_acc_direct(const uint32_t n, const uint32_t k,
                            const uint32_t d, const uint32_t *comb,
                            const bic_ctx_t ctx);

void colex_generic_unrank(uint32_t *rop, const uint32_t n, const uint32_t k,
                          const uint32_t d, const uintx r, const bic_ctx_t ctx);

uintx colex_generic_rank(const uint32_t n, const uint32_t k, const uint32_t d,
                         const uint32_t *comb, const bic_ctx_t ctx);

const bic_iter_strategy_t *bic_strat_colex(void);

#endif
