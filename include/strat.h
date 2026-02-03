#ifndef STRAT_H
#define STRAT_H

#include "common.h"

void mingen(const uint32_t m, uint32_t *n, const uint32_t k, uint32_t *d);

void minver(const uint32_t m, uint32_t *n, const uint32_t k, uint32_t *d);

void gen_params_random(const uint32_t m, uint32_t *n, const uint32_t k,
                       uint32_t *d);

#endif
