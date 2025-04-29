#ifndef GRAY_H
#define GRAY_H

#include "common.h"
#include "math.h"

void gray(uint32_t *rop, uint16_t n, uint16_t k, uint16_t d, uintx r);

uintx gray_rank(const uint16_t n, const uint16_t k, const uint16_t d,
                const uint32_t *comb);

#endif
