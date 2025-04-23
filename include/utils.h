#ifndef UTILS_H
#define UTILS_H

#include <assert.h>
#include <time.h>

#include "common.h"

#define PERF(total_time, total_cycles, logic)                                  \
  struct timespec time_start;                                                  \
  struct timespec time_stop;                                                   \
                                                                               \
  clock_gettime(CLOCK_MONOTONIC_RAW, &time_start);                             \
  uint64_t cycle_start = cycles();                                             \
                                                                               \
  logic;                                                                       \
                                                                               \
  uint64_t cycle_stop = cycles();                                              \
  clock_gettime(CLOCK_MONOTONIC_RAW, &time_stop);                              \
                                                                               \
  total_time +=                                                                \
      (long double)(time_stop.tv_sec - time_start.tv_sec) * NS_TO_SEC +        \
      (long double)(time_stop.tv_nsec - time_start.tv_nsec);                   \
  total_cycles += cycle_stop - cycle_start;

// from https://stackoverflow.com/a/40245287
typedef struct {
  const void *key;
  void *last_visited;
} bsearch_insertion_state;

uint64_t cycles(void);

uint32_t min(const uint32_t a, const uint32_t b);

int32_t max(const int32_t a, const int32_t b);

void check_valid_bounded_composition(const uint32_t *c, const uint16_t n,
                                     const uint16_t k, const uint16_t d);

size_t bsearch_insertion(const void *key, const void *base, size_t nel,
                         size_t width);

#endif
