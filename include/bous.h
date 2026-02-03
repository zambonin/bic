#ifndef BOUS_H
#define BOUS_H

#include "common.h"
#include "generic.h"

void bous_unrank(uint32_t *rop, const uint32_t n, const uint32_t k,
                 const uint32_t d, const uintx r, const bic_ctx_t ctx);

uintx bous_rank(const uint32_t n, const uint32_t k, const uint32_t d,
                const uint32_t *comb, const bic_ctx_t ctx);

void bous_generic_unrank(uint32_t *rop, const uint32_t n, const uint32_t k,
                         const uint32_t d, const uintx r, const bic_ctx_t ctx);

uintx bous_generic_rank(const uint32_t n, const uint32_t k, const uint32_t d,
                        const uint32_t *comb, const bic_ctx_t ctx);

const bic_iter_strategy_t *bic_strat_bous(void);

#endif
