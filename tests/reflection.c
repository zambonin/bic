#include "api.h"
#include <assert.h>
#include <stdio.h>

int main(void) {
  bic_ctx_t ctx = bic_ctx_init();
  unsigned checked = 0, rejected = 0;
  for (uint32_t d = 0; d <= 6; ++d) {
    uint64_t row[61] = {1};
    for (uint32_t k = 1; k <= 10; ++k) {
      uint64_t next[61] = {0};
      for (uint32_t n = 0; n <= k * d; ++n) {
        for (uint32_t x = 0; x <= d && x <= n; ++x) {
          next[n] += row[n - x];
        }
      }
      for (uint32_t n = 0; n <= k * d; ++n) {
        uintx count;
        if (bic_count_checked(&count, n, k, d)) {
          assert(count == next[n]);
          assert(compute_bic(n, k, d, ctx) == count);
          assert(compute_bic_no_cache(k * d - n, k, d) == count);
          ++checked;
        } else {
          ++rejected;
        }
      }
      for (uint32_t n = 0; n < 61; ++n) {
        row[n] = next[n];
      }
    }
  }
  intx terms[2] = {0};
  assert(compute_bic_with_sums(5, 3, 2, terms, ctx) == 3);
  assert(terms[0] == 21 && terms[1] == -18);
  assert(compute_bic(5, 3, 2, ctx) == 3);
  assert(compute_bic_no_cache(0, 0, 2) == 1);
  assert(compute_bic_no_cache(1, 0, 2) == 0);
  assert(compute_bic_no_cache(100, 4, 2) == 0);
#if !defined(BITINT) || BITINT >= 64

  assert(compute_bic_no_cache(1, 3, INT32_MAX) == 3);
  assert(compute_bic_no_cache(INT32_MAX, 2, INT32_MAX) == (uintx)INT32_MAX + 1);
#endif
  bic_ctx_destroy(ctx);
  printf("reflection: %u exact counts checked; %u unsupported counts\n",
         checked, rejected);
  return 0;
}
