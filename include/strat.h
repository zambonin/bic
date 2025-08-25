#ifndef STRAT_H
#define STRAT_H

#include "common.h"

// §4 of 10.1007/s13389-021-00264-9
void mingen(const uint32_t m, uint32_t *n, const uint32_t k, uint32_t *d);

// §2.3 of 10.1007/s13389-021-00264-9
void minver(const uint32_t m, uint32_t *n, const uint32_t k, uint32_t *d);

void gen_params_random(const uint32_t m, uint32_t *n, const uint32_t k,
                       uint32_t *d);

#endif // STRAT_H
