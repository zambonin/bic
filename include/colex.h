#ifndef COLEX_H
#define COLEX_H

#include "common.h"
#include "math.h"

void colex(uint32_t *rop, uint16_t n, uint16_t k, uint16_t d, uintx r);

void colex_bs(uint32_t *rop, uint16_t n, uint16_t k, uint16_t d, uintx r);

void colex_part(uint32_t *rop, uint16_t n, uint16_t k, uint16_t d, uintx r);

void colex_dbcs(uint32_t *rop, uint16_t n, uint16_t k, uint16_t d, uintx r);

uintx colex_rank(const uint16_t n, const uint16_t k, const uint16_t d,
                 const uint32_t *comb);

#endif
