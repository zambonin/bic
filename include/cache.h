#ifndef CACHE_H
#define CACHE_H

#include <assert.h>
#include <stdlib.h>

#include "common.h"

#define BUILD_CACHE(rows, cols, _int, var, type, logic)                        \
  uint32_t n_rows = rows;                                                      \
  uint32_t n_cols = cols;                                                      \
  var##_cols = n_cols;                                                         \
                                                                               \
  var = (_int *)calloc(n_rows * n_cols, sizeof(_int));                         \
  assert(var != NULL);                                                         \
                                                                               \
  long double total_time = 0;                                                  \
  long double total_cycles = 0;                                                \
                                                                               \
  PERF(total_time, total_cycles, logic, cache);

#define GET_CACHE(var, i, j) *(var + ((i) * var##_cols) + (j))

#define GET_CACHE_OR_CALC(type, var, math)                                     \
  if (cache_type >= type) {                                                    \
    if (cache_type == type && print_type == PRINT_ACCESS) {                    \
      (GET_CACHE(access_pattern, n, k))++;                                     \
    }                                                                          \
    return GET_CACHE(var, n, k);                                               \
  }                                                                            \
  return math(n, k, d);

enum {
  NO_CACHE = 0,
  BIN_CACHE = 1,
  COMB_CACHE = 2,
  ACC_COMB_CACHE = 3,
};

extern int cache_type;
extern uintx *bin_cache;
extern uint16_t bin_cache_cols;
extern uintx *comb_cache;
extern uint16_t comb_cache_cols;
extern uintx **acc_cache;
extern uint16_t acc_cache_rows;
extern uint16_t acc_cache_cols;

void build_cache(const uint16_t n, const uint16_t k, const uint16_t d);

void free_cache();

#endif
