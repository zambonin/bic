#ifndef CACHE_H
#define CACHE_H

#include "types.h"

typedef struct {
  uint16_t *columns;
  size_t count;
  int16_t *col_map;
  size_t col_map_size;
  bool active;
} sparse_cols_t;

sparse_cols_t sparse_cols_init(uint32_t k, const bic_ctx_t ctx);

void sparse_cols_filter(sparse_cols_t *sc, uint32_t min_col, uint32_t max_col);

void sparse_cols_transfer(sparse_cols_t *sc, cache_t *c);

void sparse_cols_free(sparse_cols_t *sc);

uintx *bin_cache_get_ptr(const uint32_t n, const uint32_t k,
                         const bic_ctx_t ctx);

uintx *comb_cache_get_ptr(const uint32_t n, const uint32_t k,
                          const bic_ctx_t ctx);

uintx *scomb_cache_get_ptr(const uint32_t n, const uint32_t k,
                           const bic_ctx_t ctx);

uintx *acc_cache_get_ptr(const uint32_t n, const uint32_t k, uint32_t *len,
                         const bic_ctx_t ctx);

uint8_t small_acc_get_bounds(const uint16_t j, const uint16_t l,
                             uint32_t *lower, uint32_t *upper,
                             const bic_ctx_t ctx);

size_t bin_cache_length(const uint32_t rows, const uint32_t cols,
                        bic_meta_t *const meta);

size_t comb_cache_length(const uint32_t rows, const uint32_t cols);

size_t scomb_cache_length(const uint32_t rows, const uint32_t cols,
                          const uint32_t k, const uint32_t d,
                          const uint8_t level, bic_meta_t *meta,
                          const uint16_t *k_values, const bic_ctx_t ctx);

size_t acc_cache_length(const uint32_t rows, const uint32_t cols,
                        const uint32_t d, bic_meta_t *const meta);

size_t small_acc_cache_length(const uint32_t rows, const uint32_t cols,
                              const uint32_t n, const uint32_t k,
                              const uint32_t d, const uint8_t level,
                              bic_meta_t *const meta, const uint16_t *k_values,
                              const bic_ctx_t ctx);

uint8_t build_cache_rec(const bic_cache_t type, const uint32_t n,
                        const uint32_t k, const uint32_t d, bic_ctx_t ctx);

void free_all_caches(bic_ctx_t ctx);

#endif
