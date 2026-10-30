#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "extra.h"

int main(void) {
  assert(compute_bic_no_cache(0, 3, 4) == 1);
  assert(random_rank(0, 3, 4) == 0);

  const uintx count = compute_bic_no_cache(4, 9, 4);
  assert(count == 495);

  srandom(0);
  bool sampled_upper_ranks = false;
  for (uint32_t i = 0; i < 10000; ++i) {
    const uintx rank = random_rank(4, 9, 4);
    assert(rank < count);
    sampled_upper_ranks |= (rank > 255);
  }
  assert(sampled_upper_ranks);

  puts("sampling regressions passed");
  return 0;
}
