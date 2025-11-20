#ifndef CACHE_COMMON_H
#define CACHE_COMMON_H

#include "cache/acc.h"
#include "cache/bin.h"
#include "cache/sparse_comb.h"
#include "math.h"
#include "types.h"

static inline bool cache_allocation_fits(const size_t count,
                                         const size_t element_size) {
  return element_size != 0 && count <= SIZE_MAX / element_size;
}

static inline bool cache_product_fits(const size_t a, const size_t b,
                                      const size_t element_size) {
  return b == 0 || (cache_allocation_fits(a, b) &&
                    cache_allocation_fits(a * b, element_size));
}

void *cache_get_element(const cache_t *cache, const uint32_t row,
                        const uint32_t col);

uint8_t build_cache_rec(const bic_cache_t type, const uint32_t n,
                        const uint32_t k, const uint32_t d, bic_ctx_t ctx);

void free_all_caches(bic_ctx_t ctx);

#endif
