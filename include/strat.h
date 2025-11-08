#ifndef STRAT_H
#define STRAT_H

#include "common.h"

// §4 of 10.1007/s13389-021-00264-9
void mingen(const uint16_t m, uint16_t *n, const uint16_t k, uint16_t *d);

// §2.3 of 10.1007/s13389-021-00264-9
void minver(const uint16_t m, uint16_t *n, const uint16_t k, uint16_t *d);

void gen_params_random(const uint16_t m, uint16_t *n, const uint16_t k,
                       uint16_t *d);

#endif // STRAT_H
