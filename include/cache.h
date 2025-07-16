#ifndef CACHE_H
#define CACHE_H

#include <assert.h>
#include <stdlib.h>

#include "common.h"

#define GET_CACHE_BIN(row, col)                                                \
  (*(uintx *)cache_get_element(&bin_cache_t, row, col))

#define GET_CACHE_COMB(row, col)                                               \
  (*(uintx *)cache_get_element(&comb_cache_t, row, col))

#define GET_CACHE_SCOMB(row, col)                                              \
  (*(uintx **)cache_get_element(&scomb_cache_t, row, col))

#define GET_CACHE_ACC(row, col)                                                \
  (*(uintx **)cache_get_element(&acc_cache_t, row, col))

#define GET_CACHE_OR_CALC(type, logic, math)                                   \
  if (cache_type >= type) {                                                    \
    return logic;                                                              \
  }                                                                            \
  return math(n, k, d);

enum {
  NO_CACHE = 0,
  BIN_CACHE = 1,
  COMB_CACHE = 2,
  SMALL_COMB_CACHE = 3,
  ACC_COMB_CACHE = 4,
  SENTINEL_LENGTH = 5,
};

extern int cache_type;

typedef struct {
  void *data;
  uint32_t rows;
  uint32_t cols;
  size_t elem_size;
  char *name;
  uint8_t type;
  size_t total_size;
} cache_t;

extern cache_t bin_cache_t;
extern cache_t comb_cache_t;
extern cache_t scomb_cache_t;
extern cache_t acc_cache_t;

void generic_setup_cache(cache_t *cache, const uint32_t rows,
                         const uint32_t cols, const size_t elem_size,
                         char *name, uint8_t type);

void *cache_get_element(const cache_t *cache, const uint32_t row,
                        const uint32_t col);

typedef void (*build_cache_funcptr_t)(const uint16_t n, const uint16_t k,
                                      const uint16_t d);

void build_caches(const uint16_t n, const uint16_t k, const uint16_t d);
void bin_build_cache(const uint16_t n, const uint16_t k, const uint16_t d);
void comb_build_cache(const uint16_t n, const uint16_t k, const uint16_t d);
void scomb_build_cache(const uint16_t n, const uint16_t k, const uint16_t d);
void acc_build_cache(const uint16_t n, const uint16_t k, const uint16_t d);

static const build_cache_funcptr_t cache_builders[SENTINEL_LENGTH] = {
    build_caches,      bin_build_cache, comb_build_cache,
    scomb_build_cache, acc_build_cache,
};

typedef void (*free_cache_funcptr_t)(void);

void free_caches();
void bin_free_cache();
void comb_free_cache();
void scomb_free_cache();
void acc_free_cache();

static const free_cache_funcptr_t cache_demolishers[SENTINEL_LENGTH] = {
    free_caches,      bin_free_cache, comb_free_cache,
    scomb_free_cache, acc_free_cache,
};

#endif
