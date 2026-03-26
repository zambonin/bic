#ifndef GENERIC_COLEX_H
#define GENERIC_COLEX_H

#include "generic.h"

void colex_generic_unrank(uint32_t *rop, const uint32_t n, const uint32_t k,
                          const uint32_t d, const uintx r, const bic_ctx_t ctx);

uintx colex_generic_rank(const uint32_t n, const uint32_t k, const uint32_t d,
                         const uint32_t *comb, const bic_ctx_t ctx);

const bic_path_strategy_t *bic_strat_colex(void);

#endif
