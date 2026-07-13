#ifndef SEARCH_H
#define SEARCH_H

#include "common.h"

void linear_search(uint32_t *rop, const uint32_t lo, const uint32_t hi,
                   bool (*p)(const uint32_t val, const void *ctx),
                   const void *ctx);

void binary_search_first_true(uint32_t *rop, const uint32_t lo,
                              const uint32_t hi,
                              bool (*p)(const uint32_t val, const void *ctx),
                              const void *ctx);

void exponential_search(uint32_t *rop, const uint32_t lo, const uint32_t hi,
                        bool (*p)(const uint32_t val, const void *ctx),
                        const void *ctx);

size_t uintx_upper_bound(const uintx *values, size_t count, const uintx *key);

#endif
