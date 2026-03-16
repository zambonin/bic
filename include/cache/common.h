#ifndef CACHE_COMMON_H
#define CACHE_COMMON_H

#include "cache/acc.h"
#include "cache/bin.h"
#include "cache/comb.h"
#include "cache/sparse_acc.h"
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

static inline uint32_t get_mapped_col(const cache_t *c, uint32_t k) {
  if (c->col_map) {
    if (k >= c->col_map_size || c->col_map[k] == -1) {
      return UINT32_MAX;
    }
    return (uint32_t)c->col_map[k];
  }
  return k;
}
uint8_t init_sparse_structures(const uint32_t k, const bic_ctx_t ctx,
                               uint32_t **H_out, size_t *h_size_out,
                               cache_t *c);

uint8_t build_cache_rec(const bic_cache_t type, const uint32_t n,
                        const uint32_t k, const uint32_t d, bic_ctx_t ctx);

void free_all_caches(bic_ctx_t ctx);

#endif
