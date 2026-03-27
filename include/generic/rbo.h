#ifndef GENERIC_RBO_H
#define GENERIC_RBO_H

#include "generic.h"

void rbo_generic_unrank(uint32_t *rop, const uint32_t n, const uint32_t k,
                        const uint32_t d, const uintx r, const bic_ctx_t ctx);

uintx rbo_generic_rank(const uint32_t n, const uint32_t k, const uint32_t d,
                       const uint32_t *comb, const bic_ctx_t ctx);

const bic_path_strategy_t *bic_strat_rbo(void);

#endif
