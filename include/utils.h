#ifndef UTILS_H
#define UTILS_H

#include "common.h"

void linear_search(uint16_t *rop, const uint16_t lo, const uint16_t hi,
                   bool (*p)(const uint16_t, const void *), const void *ctx);

size_t bsearch_insertion(const void *key, const void *base, size_t nel,
                         size_t width);

uint16_t bits_fit_bic(const uint16_t n, const uint16_t k, const uint16_t d);

uintx random_rank(const uint16_t n, const uint16_t k, const uint16_t d);

#endif
