#ifndef RBO_H
#define RBO_H

#include "common.h"

void inner_rbo_unrank(uint32_t *rop, const uint16_t n, const uint16_t k,
                      const uint16_t d, const uintx r, const uint16_t start,
                      const bic_ctx_t ctx);

// §4.5 of 10.1007/978-3-031-22969-5_1
void rbo_unrank(uint32_t *rop, const uint16_t n, const uint16_t k,
                const uint16_t d, const uintx r, const bic_ctx_t ctx);

// §4.2 of 10.1007/978-3-031-22969-5_1
uintx rbo_rank(const uint16_t n, const uint16_t k, const uint16_t d,
               const uint32_t *comb, const bic_ctx_t ctx);

#endif
