#ifndef RB_SPIRAL_H
#define RB_SPIRAL_H

#include "common.h"

void inner_rb_spiral_unrank(uint32_t *rop, const uint32_t n, const uint32_t k,
                            const uint32_t d, const uintx r,
                            const uint32_t start, const bic_ctx_t ctx);

void rb_spiral_unrank(uint32_t *rop, const uint32_t n, const uint32_t k,
                      const uint32_t d, const uintx r, const bic_ctx_t ctx);

uintx rb_spiral_rank(const uint32_t n, const uint32_t k, const uint32_t d,
                     const uint32_t *comb, const bic_ctx_t ctx);

#endif
