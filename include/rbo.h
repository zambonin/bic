#ifndef RBO_H
#define RBO_H

#include "common.h"
#include "generic.h"

void rbo_unrank(uint32_t *rop, const uint32_t n, const uint32_t k,
                const uint32_t d, const uintx r, const bic_ctx_t ctx);

uintx rbo_rank(const uint32_t n, const uint32_t k, const uint32_t d,
               const uint32_t *comb, const bic_ctx_t ctx);

void rbo_generic_unrank(uint32_t *rop, const uint32_t n, const uint32_t k,
                        const uint32_t d, const uintx r, const bic_ctx_t ctx);

uintx rbo_generic_rank(const uint32_t n, const uint32_t k, const uint32_t d,
                       const uint32_t *comb, const bic_ctx_t ctx);

const bic_iter_strategy_t *bic_strat_rbo(void);

void rbo_mo_unrank(uint32_t *rop, const uint32_t n, const uint32_t k,
                   const uint32_t d, const uintx r, const bic_ctx_t ctx);

uintx rbo_mo_rank(const uint32_t n, const uint32_t k, const uint32_t d,
                  const uint32_t *comb, const bic_ctx_t ctx);

void rbo_mo_generic_unrank(uint32_t *rop, const uint32_t n, const uint32_t k,
                           const uint32_t d, const uintx r,
                           const bic_ctx_t ctx);

uintx rbo_mo_generic_rank(const uint32_t n, const uint32_t k, const uint32_t d,
                          const uint32_t *comb, const bic_ctx_t ctx);

const bic_iter_strategy_t *bic_strat_rbo_mo(void);

#endif
