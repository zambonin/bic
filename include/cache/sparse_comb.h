#ifndef CACHE_SPARSE_COMB_H
#define CACHE_SPARSE_COMB_H

#include "types.h"

uint8_t scomb_build_cache(const uint32_t n, const uint32_t k, const uint32_t d,
                          bic_ctx_t ctx);
void scomb_free_cache(bic_ctx_t ctx);
uintx scomb_cache_get(const uint32_t n, const uint32_t k, const uint32_t d,
                      const bic_ctx_t ctx);

uintx *scomb_cache_get_ptr(const uint32_t n, const uint32_t k,
                           const bic_ctx_t ctx);

size_t scomb_cache_length(const uint32_t rows, const uint32_t cols,
                          const uint32_t d, const uint8_t level,
                          bic_meta_t *meta, const bic_ctx_t ctx);

#endif
