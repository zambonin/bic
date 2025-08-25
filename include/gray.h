#ifndef GRAY_H
#define GRAY_H

#include "common.h"

void gray_unrank(const bic_ctx_t *ctx, uint32_t *rop, const uint16_t n,
                 const uint16_t k, const uint16_t d, const uintx r);

uintx gray_rank(const bic_ctx_t *ctx, const uint16_t n, const uint16_t k,
                const uint16_t d, const uint32_t *comb);

#endif
