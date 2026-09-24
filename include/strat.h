#ifndef STRAT_H
#define STRAT_H

#include "common.h"
#include <stdbool.h>

typedef enum {
  BIC_COUNT_UNSUPPORTED = -1,
  BIC_COUNT_BELOW = 0,
  BIC_COUNT_AT_LEAST = 1,
} bic_count_comparison_t;

bic_count_comparison_t bic_compare_2_pow_m(const uint32_t m, const uint32_t n,
                                           const uint32_t k, const uint32_t d);

bool bic_geq_2_pow_m(const uint32_t m, const uint32_t n, const uint32_t k,
                     const uint32_t d);

// §4 of 10.1007/s13389-021-00264-9
void mingen(const uint32_t m, uint32_t *n, const uint32_t k, uint32_t *d);

// §2.3 of 10.1007/s13389-021-00264-9
void minver(const uint32_t m, uint32_t *n, const uint32_t k, uint32_t *d);

void gen_params_random(const uint32_t m, uint32_t *n, const uint32_t k,
                       uint32_t *d);

#endif // STRAT_H
