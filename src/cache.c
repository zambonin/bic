#include "cache.h"
#include "math.h"
#include "utils.h"

#if defined(DHAT)
#include <valgrind/dhat.h>
#endif

cache_t bin_cache_t;
cache_t rbo_cache_t;
cache_t comb_cache_t;
cache_t acc_cache_t;
cache_t scomb_cache_t;

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

void after_cache_build(cache_t *cache) {
  (void)cache;
#if defined(DHAT)
  DHAT_HISTOGRAM_MEMORY(cache->data);
#endif
}

void bin_build_cache(const uint16_t n, const uint16_t k, const uint16_t d) {
  (void)d;
  generic_setup_cache(&bin_cache_t, n + k + 1, k, sizeof(uintx), (char *)"bin",
                      BIN_CACHE);

  GET_CACHE_BIN(0, 0) = 1;
  for (uint32_t row = 1; row < bin_cache_t.rows; ++row) {
    GET_CACHE_BIN(row, 0) = 1;
    for (uint16_t col = 1; col <= min(row, bin_cache_t.cols - 1); ++col) {
      GET_CACHE_BIN(row, col) = (uintx)GET_CACHE_BIN(row - 1, col - 1) +
                                (uintx)GET_CACHE_BIN(row - 1, col);
    }
  }

  after_cache_build(&bin_cache_t);
}

void rbo_build_cache(const uint16_t n, const uint16_t k, const uint16_t d) {
  if (k == 0) return; // k can't be 0 to calculate the halving set

  if (k == 1) { // halving set is null when k = 1 (cache is just one column)
    generic_setup_cache(&rbo_cache_t, n + 1, k + 1, sizeof(uintx),
                      (char *)"rbo", RBO_COMB_CACHE);
    rbo_cache_t.col_of_original = NULL;
    GET_CACHE_RBO(0, 0) = 1;
    for (uint16_t row = 0; row < rbo_cache_t.rows; ++row) {
        GET_CACHE_RBO(row, 1) = inner_bic(row, 1, d);
    }
    return;
  }

  uint16_t L = lg(k);
  uint16_t *H = (uint16_t *)malloc(2 * L * sizeof(uint16_t)); // 2 * L is an upper bound on the size of the halving set
  size_t h_size = 0; // h_size will be the actual size of the halving set
  calc_hset(k, L, H, &h_size); 

  generic_setup_cache(&rbo_cache_t, n + 1, (uint16_t)h_size, sizeof(uintx),
                      (char *)"rbo", RBO_COMB_CACHE);

  int16_t *col_of_original = (int16_t *)malloc((k + 1) * sizeof(int16_t));
  for (int i = 0; i < k + 1; i++) col_of_original[i] = -1; // -1 means that the col doesn´t exists in H
  for (size_t idx = 0; idx < h_size; idx++) col_of_original[H[idx]] = idx; // maps the columns
  rbo_cache_t.col_of_original = col_of_original;

  for (uint32_t row = 0; row < rbo_cache_t.rows; ++row) {
    for (size_t col = 0; col < h_size; ++col) {
      GET_CACHE_RBO(row, col) = inner_bic(row, H[col], d);
    }
  }

  free(H);
  after_cache_build(&rbo_cache_t);
}

void comb_build_cache(const uint16_t n, const uint16_t k, const uint16_t d) {
  generic_setup_cache(&comb_cache_t, n + 1, k + 1, sizeof(uintx),
                      (char *)"comb", COMB_CACHE);

  GET_CACHE_COMB(0, 0) = 1;
  for (uint16_t row = 0; row < comb_cache_t.rows; ++row) {
    for (uint16_t col = 1; col < comb_cache_t.cols; ++col) {
      GET_CACHE_COMB(row, col) = inner_bic(row, col, d);
    }
  }

  after_cache_build(&comb_cache_t);
}

void scomb_build_cache(const uint16_t n, const uint16_t k, const uint16_t d) {
  generic_setup_cache(&scomb_cache_t, 1, k - 1, sizeof(uintx *),
                      (char *)"scomb", SMALL_COMB_CACHE);

  double variance = d * (d + 2) / 12;
  uint8_t level = 4;

  for (uint16_t col = 0; col < scomb_cache_t.cols; ++col) {
    uint16_t j = col + 1;
    double mean = j * n / k;
    double stddev = asqrt(j * variance * (k - j) / k);
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

  after_cache_build(&scomb_cache_t);
}

void acc_build_cache(const uint16_t n, const uint16_t k, const uint16_t d) {
  generic_setup_cache(&acc_cache_t, n + 1, k, sizeof(uintx *), (char *)"acc",
                      ACC_COMB_CACHE);

  for (uint16_t row = 0; row < acc_cache_t.rows; ++row) {
    for (uint16_t col = 0; col < acc_cache_t.cols; ++col) {
      GET_CACHE_ACC(row, col) = inner_acc(row, col, d);
    }
  }

  after_cache_build(&acc_cache_t);
}

void build_caches(const uint16_t n, const uint16_t k, const uint16_t d) {
  for (uint8_t i = 1; i <= cache_type; ++i) {
    cache_builders[i](n, k, d);
  }
}

void bin_free_cache() { free(bin_cache_t.data); }

void rbo_free_cache() {
  free(rbo_cache_t.data);
  if (rbo_cache_t.col_of_original != NULL) free(rbo_cache_t.col_of_original);
}

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
