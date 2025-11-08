#ifndef CACHE_H
#define CACHE_H

#include "types.h"

void *cache_get_element(const cache_t *cache, const uint32_t row,
                        const uint32_t col);

uint8_t build_cache_rec(const bic_cache_t type, const uint16_t n,
                        const uint16_t k, const uint16_t d, bic_ctx_t ctx);

void free_all_caches(bic_ctx_t ctx);

#endif
