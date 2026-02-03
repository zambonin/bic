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

size_t bsearch_insertion(const void *key, const void *base, size_t nel,
                         size_t width);

#endif
