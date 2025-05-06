#include "cache.h"
#include "io.h"
#include "math.h"
#include "utils.h"

cache_t bin_cache_t;
cache_t comb_cache_t;
cache_t acc_cache_t;
cache_t scomb_cache_t;
cache_t access_pattern_cache_t;

int cache_type = NO_CACHE;

void *cache_get_element(const cache_t *cache, const uint32_t row,
                        const uint32_t col) {
  uint32_t index = row * cache->cols + col;
  size_t offset = index * cache->elem_size;
  return (char *)cache->data + offset;
}

void generic_setup_cache(cache_t *cache, const uint32_t rows,
                         const uint32_t cols, const size_t elem_size,
                         char *name, uint8_t type) {
  cache->rows = rows;
  cache->cols = cols;
  cache->elem_size = elem_size;
  cache->total_size = cache->rows * cache->cols * cache->elem_size;
  cache->name = name;
  cache->type = type;

  cache->data = calloc(cache->rows * cache->cols, cache->elem_size);
  assert(cache->data != NULL);
}

void print_cache_data(cache_t *cache, long double total_time,
                      long double total_cycles) {
  if (cache_type == cache->type && print_type == PRINT_ACCESS) {
    generic_setup_cache(&access_pattern_cache_t, cache->rows, cache->cols,
                        sizeof(uint32_t), (char *)"access", 0);
  } else if (print_type == PRINT_BUILD) {
    printf(
        "c = %6s, size = %12zu, avg time = %14.2Lf ns, avg cycles = %14.2Lf\n",
        cache->name, cache->total_size, total_time, total_cycles);
  }
}

void bin_build_cache(const uint16_t n, const uint16_t k, const uint16_t d) {
  (void)d;
  generic_setup_cache(&bin_cache_t, n + k + 1, k, sizeof(uintx), (char *)"bin", BIN_CACHE);

  PERF_CACHE(bin_cache_t, {
    GET_CACHE_BIN(0, 0) = 1;
    for (uint32_t row = 1; row < bin_cache_t.rows; ++row) {
      GET_CACHE_BIN(row, 0) = 1;
      for (uint16_t col = 1; col <= min(row, bin_cache_t.cols - 1); ++col) {
        GET_CACHE_BIN(row, col) = (uintx)GET_CACHE_BIN(row - 1, col - 1) +
                                  (uintx)GET_CACHE_BIN(row - 1, col);
      }
    }
  });
}

void comb_build_cache(const uint16_t n, const uint16_t k, const uint16_t d) {
  generic_setup_cache(&comb_cache_t, n + 1, k + 1, sizeof(uintx),
                      (char *)"comb", COMB_CACHE);

  PERF_CACHE(comb_cache_t, {
    GET_CACHE_COMB(0, 0) = 1;
    for (uint16_t row = 0; row < comb_cache_t.rows; ++row) {
      for (uint16_t col = 1; col < comb_cache_t.cols; ++col) {
        GET_CACHE_COMB(row, col) = inner_bic(row, col, d);
      }
    }
  });
}

void scomb_build_cache(const uint16_t n, const uint16_t k, const uint16_t d) {
  generic_setup_cache(&scomb_cache_t, 1, k - 1, sizeof(uintx *),
                      (char *)"scomb", SMALL_COMB_CACHE);

  PERF_CACHE(scomb_cache_t, {
    double variance = d * (d + 2) / 12;
    uint8_t level = 4;

    for (uint16_t col = 0; col < scomb_cache_t.cols; ++col) {
      uint16_t j = col + 1;
      double mean = j * n / k;
      double stddev = sqrt(j * variance * (k - j) / k);
      uint16_t left = max(mean - level * stddev, 0);
      uint16_t right = min(mean + level * stddev, n);

      uint16_t length = 2 + (right - left + 1);
      uintx *part = (uintx *)calloc(length, sizeof(uintx));

      part[0] = left;
      part[1] = right;
      for (uint16_t i = 2; i < length; ++i) {
        part[i] = inner_bic(left + i - 2, j, d);
      }

      GET_CACHE_SCOMB(0, col) = part;
      scomb_cache_t.total_size += length * sizeof(uintx);
    }
  });
}

void acc_build_cache(const uint16_t n, const uint16_t k, const uint16_t d) {
  generic_setup_cache(&acc_cache_t, n + 1, k, sizeof(uintx *), (char *)"acc", ACC_COMB_CACHE);

  PERF_CACHE(acc_cache_t, {
    for (uint16_t row = 0; row < acc_cache_t.rows; ++row) {
      for (uint16_t col = 0; col < acc_cache_t.cols; ++col) {
        GET_CACHE_ACC(row, col) = inner_acc(row, col, d);
      }
    }
  });
}

void build_caches(const uint16_t n, const uint16_t k, const uint16_t d) {
  for (uint8_t i = 1; i <= cache_type; ++i) {
    cache_builders[i](n, k, d);
  }
}

void bin_free_cache() { free(bin_cache_t.data); }

void comb_free_cache() { free(comb_cache_t.data); }

void scomb_free_cache() {
  for (uint16_t j = 0; j < scomb_cache_t.cols; ++j) {
    free(GET_CACHE_SCOMB(0, j));
  }
  free(scomb_cache_t.data);
}

void acc_free_cache() {
  for (uint16_t i = 0; i < acc_cache_t.rows; ++i) {
    for (uint16_t j = 0; j < acc_cache_t.cols; ++j) {
      free(GET_CACHE_ACC(i, j));
    }
  }
  free(acc_cache_t.data);
}

void free_caches() {
  for (uint8_t i = 1; i <= cache_type; ++i) {
    cache_demolishers[i]();
  }
}
