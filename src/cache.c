#include "cache.h"

int cache_type = NO_CACHE;
uintx *bin_cache;
uint16_t bin_cache_cols = 0;
uintx *comb_cache;
uint16_t comb_cache_cols = 0;
uintx **acc_cache;
uint16_t acc_cache_rows = 0;
uint16_t acc_cache_cols = 0;

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
        GET_CACHE(acc_cache, row, col) = inner_acc(row, col, d);
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
  if (cache_type >= ACC_COMB_CACHE) {
    build_acc_cache(n, k, d);
  }
}
