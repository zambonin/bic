#ifndef CACHE_COMB_H
#define CACHE_COMB_H

#include "types.h"

uint8_t comb_build_cache(const uint32_t n, const uint32_t k, const uint32_t d,
                         bic_ctx_t ctx);
void comb_free_cache(bic_ctx_t ctx);
uintx comb_cache_get(const uint32_t n, const uint32_t k, const uint32_t d,
                     const bic_ctx_t ctx);

uintx *comb_cache_get_ptr(const uint32_t n, const uint32_t k,
                          const bic_ctx_t ctx);

size_t comb_cache_length(const uint32_t rows, const uint32_t cols);

#endif
