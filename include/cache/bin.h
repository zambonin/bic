#ifndef CACHE_BIN_H
#define CACHE_BIN_H

#include "types.h"

uint8_t bin_build_cache(const uint32_t n, const uint32_t k, bic_ctx_t ctx);
void bin_free_cache(bic_ctx_t ctx);
uintx bin_cache_get(const uint32_t n, const uint32_t k, const bic_ctx_t ctx);

uintx *bin_cache_get_ptr(const uint32_t n, const uint32_t k,
                         const bic_ctx_t ctx);

size_t bin_cache_length(const uint32_t rows, const uint32_t cols,
                        bic_meta_t *const meta);

#endif
