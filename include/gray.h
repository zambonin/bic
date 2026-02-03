#ifndef GRAY_H
#define GRAY_H

#include "common.h"
#include "generic.h"

void gray_unrank(uint32_t *rop, const uint32_t n, const uint32_t k,
                 const uint32_t d, const uintx r, const bic_ctx_t ctx);

uintx gray_rank(const uint32_t n, const uint32_t k, const uint32_t d,
                const uint32_t *comb, const bic_ctx_t ctx);

void gray_generic_unrank(uint32_t *rop, const uint32_t n, const uint32_t k,
                         const uint32_t d, const uintx r, const bic_ctx_t ctx);

uintx gray_generic_rank(const uint32_t n, const uint32_t k, const uint32_t d,
                        const uint32_t *comb, const bic_ctx_t ctx);

const bic_iter_strategy_t *bic_strat_gray(void);

#endif
