#ifndef CACHE_ACC_H
#define CACHE_ACC_H

#include "types.h"

uint8_t acc_build_cache(const uint32_t n, const uint32_t k, const uint32_t d,
                        bic_ctx_t ctx);
void acc_free_cache(bic_ctx_t ctx);
uint32_t acc_cache_get(uintx *rop, const uint32_t n, const uint32_t k,
                       const uint32_t d, const bic_ctx_t ctx);
uintx dir_cache_get(const uint32_t n, const uint32_t k, const uint32_t d,
                    const uint32_t l, const bic_ctx_t ctx);

uintx *acc_cache_get_ptr(const uint32_t n, const uint32_t k, uint32_t *len,
                         const bic_ctx_t ctx);

size_t acc_cache_length(const uint32_t rows, const uint32_t cols,
                        const uint32_t d, bic_meta_t *const meta);

#endif
