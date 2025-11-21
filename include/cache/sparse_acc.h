#ifndef CACHE_SPARSE_ACC_H
#define CACHE_SPARSE_ACC_H

#include "types.h"

uint8_t small_acc_build_cache(const uint32_t n, const uint32_t k,
                              const uint32_t d, const bic_ctx_t ctx);
void small_acc_free_cache(bic_ctx_t ctx);
uintx small_acc_dir_get(const uint32_t n, const uint32_t k, const uint32_t d,
                        const uint32_t l, const bic_ctx_t ctx);

uint8_t small_acc_get_bounds(const uint32_t j, const uint32_t l,
                             uint32_t *lower, uint32_t *upper,
                             const bic_ctx_t ctx);

size_t small_acc_cache_length(const uint32_t rows, const uint32_t cols,
                              const uint32_t n, const uint32_t k,
                              const uint32_t d, const uint8_t level,
                              bic_meta_t *const meta, const bic_ctx_t ctx);

#endif
