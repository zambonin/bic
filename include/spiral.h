#ifndef SPIRAL_H
#define SPIRAL_H

#include "common.h"
#include "generic.h"

void spiral_unrank(uint32_t *rop, const uint32_t n, const uint32_t k,
                   const uint32_t d, const uintx r, const bic_ctx_t ctx);

uintx spiral_rank(const uint32_t n, const uint32_t k, const uint32_t d,
                  const uint32_t *comb, const bic_ctx_t ctx);

void spiral_generic_unrank(uint32_t *rop, const uint32_t n, const uint32_t k,
                           const uint32_t d, const uintx r,
                           const bic_ctx_t ctx);

uintx spiral_generic_rank(const uint32_t n, const uint32_t k, const uint32_t d,
                          const uint32_t *comb, const bic_ctx_t ctx);

const bic_iter_strategy_t *bic_strat_spiral(void);

#endif
