#include <assert.h>
#include <errno.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "api.h"
#include "strat.h"

#if defined(BITINT) && BITINT == 16
static void narrow_oracle(void) {
  unsigned checks = 0, classes = 0;
  for (uint32_t k = 2; k <= 80; ++k) {
    uint32_t lo = 0, hi = INT32_MAX;
    while (lo < hi) {
      const uint32_t mid = lo + (hi - lo + 1) / 2;
      if (bic_params_supported(mid, k, mid, BIC_CACHE_NONE)) {
        lo = mid;
      } else {
        hi = mid - 1;
      }
    }
    const uint32_t cap = lo;
    assert(cap <= 254);
    uint32_t gn[16], gd[16], vn[16], vd[16];
    for (unsigned m = 0; m < 16; ++m) {
      gn[m] = gd[m] = vn[m] = vd[m] = UINT32_MAX;
    }
    gn[0] = gd[0] = vn[0] = vd[0] = 0;
    for (uint32_t d = 0; d <= cap; ++d) {
      uint64_t row[255] = {1}, next[255];
      for (uint32_t part = 0; part < k; ++part) {
        for (uint32_t n = 0; n <= cap; ++n) {
          next[n] = 0;
          for (uint32_t x = 0; x <= d && x <= n; ++x) {
            assert(UINT64_MAX - next[n] >= row[n - x]);
            next[n] += row[n - x];
          }
        }
        memcpy(row, next, (cap + 1) * sizeof(*row));
      }
      for (uint32_t n = 0; n <= cap; ++n) {
        ++classes;
        if (!bic_params_supported(n, k, d, BIC_CACHE_NONE)) {
          continue;
        }
        for (unsigned m = 1; m < 16; ++m) {
          if (row[n] < (UINT64_C(1) << m)) {
            continue;
          }
          if (d < gd[m] || (d == gd[m] && n < gn[m])) {
            gd[m] = d;
            gn[m] = n;
          }
          if (n < vn[m] || (n == vn[m] && d < vd[m])) {
            vn[m] = n;
            vd[m] = d;
          }
        }
      }
    }
    for (unsigned m = 0; m < 16; ++m) {
      uint32_t n, d;
      errno = EDOM;
      mingen(m, &n, k, &d);
      assert(n == gn[m] && d == gd[m]);
      assert(errno == (n == UINT32_MAX ? ERANGE : 0));
      errno = EDOM;
      minver(m, &n, k, &d);
      assert(n == vn[m] && d == vd[m]);
      assert(errno == (n == UINT32_MAX ? ERANGE : 0));
      checks += 2;
    }
  }
  assert(checks == 2528 && classes == 71203);
  printf("Parameter oracle: %u optima/failures across %u classes passed\n",
         checks, classes);
}
#endif

#if !defined(BITINT) || BITINT >= 512

static uintx choose_exact(uint32_t n, uint32_t k) {
  if (k > n) {
    return 0;
  }
  if (k > n - k) {
    k = n - k;
  }
  uintx value = 1;
  for (uint32_t i = 1; i <= k; ++i) {
    value = value * (n - k + i) / i;
  }
  return value;
}

static uintx count_exact(uint32_t n, uint32_t k, uint32_t d) {
  if ((uint64_t)n > (uint64_t)k * d) {
    return 0;
  }
  uint32_t limit = n / (d + 1);
  if (limit > k) {
    limit = k;
  }
  intx sum = 0;
  for (uint32_t i = 0; i <= limit; ++i) {
    const intx term = (intx)(choose_exact(k, i) *
                             choose_exact(n - i * (d + 1) + k - 1, k - 1));
    sum += (i % 2) ? -term : term;
  }
  assert(sum >= 0);
  return (uintx)sum;
}

static void wide_oracle(void) {
  const uint32_t fixtures[][6] = {
      {260, 10, INT32_MAX, 686334074, 2061453506, 1877716767},
      {448, 17, INT32_MAX, 298129809, 1825486398, 1366537962},
  };
  for (unsigned i = 0; i < 2; ++i) {
    const uint32_t *f = fixtures[i];
    const uint32_t m = f[0], k = f[1], gn = f[2], gd = f[3], vn = f[4],
                   vd = f[5];
    const uintx target = (uintx)1 << m;
    uint64_t peak = ((uint64_t)k * (gd - 1) + 1) / 2;
    if (peak > INT32_MAX) {
      peak = INT32_MAX;
    }
    assert(count_exact(gn, k, gd) >= target);
    assert(count_exact(gn - 1, k, gd) < target);
    assert(count_exact((uint32_t)peak, k, gd - 1) < target);
    assert(count_exact(vn, k, vd) >= target);
    assert(count_exact(vn, k, vd - 1) < target);
    assert(choose_exact(vn + k - 2, k - 1) < target);
    uint32_t n, d;
    mingen(m, &n, k, &d);
    assert(n == gn && d == gd && errno == 0);
    minver(m, &n, k, &d);
    assert(n == vn && d == vd && errno == 0);
  }
  puts("Wide parameter fixtures and failed neighbors passed");
}
#endif

int main(void) {
#if defined(BITINT) && BITINT == 16
  narrow_oracle();
#elif !defined(BITINT) || BITINT >= 512
  wide_oracle();
#else
  puts("Parameter oracle requires BitInt16 or a backend of at least 512 bits");
#endif
  return 0;
}
