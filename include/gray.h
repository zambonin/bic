#ifndef GRAY_H
#define GRAY_H

#include "common.h"

void gray_unrank(uint32_t *rop, const uint16_t n, const uint16_t k,
                 const uint16_t d, const uintx r);

uintx gray_rank(const uint16_t n, const uint16_t k, const uint16_t d,
                const uint32_t *comb);

static const order gray = {.unrank = gray_unrank, .rank = gray_rank};

#endif
