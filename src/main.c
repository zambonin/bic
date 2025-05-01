#include "cache.h"
#include "colex.h"
#include "io.h"
#include "math.h"
#include "utils.h"

int32_t main(int32_t argc, char **argv) {
  uint16_t n = 0;
  uint16_t k = 0;
  uint16_t d = 0;
  uint16_t m = 0;
  uint32_t iterations = 1;
  order ord = colex;
  strategy_func strategy = mingen;

  if (parse_args(argc, argv, &n, &k, &d, &iterations, &ord, &m, &strategy) >
      0) {
    return 1;
  }

  srandom(time(NULL));

  build_cache(n, k, d);

  long double utime = 0;
  long double ucycles = 0;

  long double rtime = 0;
  long double rcycles = 0;

  uint32_t *comp = (uint32_t *)malloc(k * sizeof(uint32_t));

  for (uint32_t it = 0; it < iterations; ++it) {
    memset(comp, 0, k * sizeof(uint32_t));

    const uintx r = random_rank(n, k, d);

    PERF(utime, ucycles, (*ord.unrank)(comp, n, k, d, r), unrank);

    PERF(rtime, rcycles, const uintx rr = (*ord.rank)(n, k, d, comp), rank);

    assert(r == rr);

    check_valid_bounded_composition(comp, n, k, d);
  }

  pprint(n, k, d, iterations, utime, ucycles, rtime, rcycles);

  free_cache();
  free(comp);

  return 0;
}
