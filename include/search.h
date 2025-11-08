#ifndef SEARCH_H
#define SEARCH_H

#include "common.h"

void linear_search(uint16_t *rop, const uint16_t lo, const uint16_t hi,
                   bool (*p)(const uint16_t val, const void *ctx),
                   const void *ctx);

size_t bsearch_insertion(const void *key, const void *base, size_t nel,
                         size_t width);

#endif
