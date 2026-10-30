#include <assert.h>
#include <stdio.h>

#include "math.h"

int main(void) {
  assert(compute_bin(10007, 1, NULL) == 10007);
  assert(compute_bin(10007, 10006, NULL) == 10007);
  assert(compute_bin(10007, 2, NULL) == 50065021);
  assert(compute_bic_no_cache(1, 10007, 1) == 10007);
  puts("binomial regressions passed");
  return 0;
}
