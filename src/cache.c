#include "cache.h"
#include "io.h"
#include "math.h"
#include "utils.h"

int cache_type = NO_CACHE;
uintx *bin_cache;
uint16_t bin_cache_cols = 0;
uint64_t bin_cache_size = 0;
uintx *comb_cache;
uint16_t comb_cache_cols = 0;
uint64_t comb_cache_size = 0;
uintx **acc_cache;
uint16_t acc_cache_rows = 0;
uint16_t acc_cache_cols = 0;
uint64_t acc_cache_size = 0;
uintx **small_comb_cache;
uint16_t small_comb_cache_cols = 0;
uint64_t small_comb_cache_size = 0;

void build_bin_cache(const uint16_t n, const uint16_t k, const uint16_t d) {
  (void)d;
  BUILD_CACHE(n + k + 1, k, uintx, bin_cache, BIN_CACHE, {
    GET_CACHE(bin_cache, 0, 0) = 1;
    for (uint32_t row = 1; row < n_rows; ++row) {
      GET_CACHE(bin_cache, row, 0) = 1;
      for (uint16_t col = 1; col <= min(row, n_cols - 1); ++col) {
        GET_CACHE(bin_cache, row, col) =
            (uintx)GET_CACHE(bin_cache, row - 1, col - 1) +
            (uintx)GET_CACHE(bin_cache, row - 1, col);
      }
    }
  });

  PRINT_CACHE_STATS(bin_cache, BIN_CACHE, "bin");
}

void build_small_comb_cache(const uint16_t n, const uint16_t k,
                            const uint16_t d) {
  BUILD_CACHE(1, k - 1, uintx*, small_comb_cache, SMALL_COMB_CACHE, {
    double variance = d * (d + 2) / 12;
    uint8_t level = 4;

    for (uint16_t col = 0; col < n_cols; ++col) {
      uint16_t j = col + 1;
      double mean = j * n / k;
      double stddev = sqrt(j * variance * (k - j) / k);
      uint16_t left = max(mean - level * stddev, 0);
      uint16_t right = min(mean + level * stddev, n);

      uint16_t length = 2 + (right - left + 1);
      uintx* part = (uintx *)calloc(length, sizeof(uintx));
      small_comb_cache[col] = part;

      part[0] = left;
      part[1] = right;
      for (uint16_t i = 2; i < length; ++i) {
        part[i] = inner_bic(left + i - 2, j, d);
      }

      small_comb_cache_size += length * sizeof(uintx);
    }
  });

  (void)total_time;
  (void)total_cycles;
}

void build_comb_cache(const uint16_t n, const uint16_t k, const uint16_t d) {
  BUILD_CACHE(n + 1, k + 1, uintx, comb_cache, COMB_CACHE, {
    GET_CACHE(comb_cache, 0, 0) = 1;
    for (uint16_t row = 0; row < n_rows; ++row) {
      for (uint16_t col = 1; col < n_cols; ++col) {
        GET_CACHE(comb_cache, row, col) = inner_bic(row, col, d);
      }
    }
  });

  PRINT_CACHE_STATS(comb_cache, COMB_CACHE, "comb");
}

void build_acc_cache(const uint16_t n, const uint16_t k, const uint16_t d) {
  BUILD_CACHE(n + 1, k, uintx *, acc_cache, ACC_COMB_CACHE, {
    for (uint16_t row = 0; row < n_rows; ++row) {
      for (uint16_t col = 0; col < n_cols; ++col) {
        uintx* part = inner_acc(row, col, d);
        GET_CACHE(acc_cache, row, col) = part;
        free(part);
      }
    }
  });

  (void)total_time;
  (void)total_cycles;
}

void build_cache(const uint16_t n, const uint16_t k, const uint16_t d) {
  if (cache_type >= BIN_CACHE) {
    build_bin_cache(n, k, d);
  }
  if (cache_type >= COMB_CACHE) {
    build_comb_cache(n, k, d);
  }
  if (cache_type >= SMALL_COMB_CACHE) {
    build_small_comb_cache(n, k, d);
  }
  if (cache_type >= ACC_COMB_CACHE) {
    build_acc_cache(n, k, d);
  }
}

void free_cache() {
  if (cache_type >= BIN_CACHE) {
    free(bin_cache);
  }
  if (cache_type >= COMB_CACHE) {
    free(comb_cache);
  }
  if (cache_type >= SMALL_COMB_CACHE) {
    for (uint16_t i = 0; i < small_comb_cache_cols; ++i) {
      uintx* part = small_comb_cache[i];
      free(part);
    }
    free(small_comb_cache);
  }
  if (cache_type >= ACC_COMB_CACHE) {
    for (uint16_t i = 0; i < acc_cache_rows; ++i) {
      for (uint16_t j = 0; j < acc_cache_cols; ++j) {
        free(GET_CACHE(acc_cache, i, j));
      }
    }
    free(acc_cache);
  }
  free(access_pattern);
}
