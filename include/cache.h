#ifndef CACHE_H
#define CACHE_H

#include "types.h"

uintx *bin_cache_get_ptr(const uint16_t n, const uint16_t k,
                         const bic_ctx_t ctx);

uintx *comb_cache_get_ptr(const uint16_t n, const uint16_t k,
                          const bic_ctx_t ctx);

uintx *scomb_cache_get_ptr(const uint16_t n, const uint16_t k,
                           const bic_ctx_t ctx);

uintx *acc_cache_get_ptr(const uint16_t n, const uint16_t k, uint16_t *len,
                         const bic_ctx_t ctx);

uint8_t small_acc_get_bounds(const uint16_t j, const uint16_t l,
                             uint32_t *lower, uint32_t *upper,
                             const bic_ctx_t ctx);

size_t bin_cache_length(const uint16_t rows, const uint16_t cols,
                        bic_meta_t *const meta);

size_t comb_cache_length(const uint16_t rows, const uint16_t cols);

size_t scomb_cache_length(const uint16_t rows, const uint16_t cols,
                          const uint16_t d, const uint8_t level,
                          bic_meta_t *meta, const bic_ctx_t ctx);

size_t acc_cache_length(const uint16_t rows, const uint16_t cols,
                        const uint16_t d, bic_meta_t *const meta);

size_t small_acc_cache_length(const uint16_t rows, const uint16_t cols,
                              const uint16_t n, const uint16_t k,
                              const uint16_t d, const uint8_t level,
                              bic_meta_t *const meta, const bic_ctx_t ctx);

uint8_t build_cache_rec(const bic_cache_t type, const uint16_t n,
                        const uint16_t k, const uint16_t d, bic_ctx_t ctx);

void free_all_caches(bic_ctx_t ctx);

#endif
