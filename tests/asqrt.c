#include <assert.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>

#include "math.h"

static void check(const double x) {
  const double y = asqrt(x);
  const double error = y - x / y;
  assert(y > 0.0 && y <= DBL_MAX);
  assert(error <= 4.0 * DBL_EPSILON * y);
  assert(error >= -4.0 * DBL_EPSILON * y);
}

int main(void) {
  assert(asqrt(-1.0) == -1.0);
  assert(asqrt(0.0) == 0.0);
  assert(asqrt(0.25) == 0.5);
  assert(asqrt(1.0) == 1.0);
  assert(asqrt(4.0) == 2.0);
  check(1e-20);
  check(2.0);
  check(10.0);
  check(1e20);

  double x = 1.0, expected = 1.0;
  while (x <= DBL_MAX / 4.0) {
    assert(asqrt(x) == expected);
    x *= 4.0;
    expected *= 2.0;
  }
  assert(asqrt(x) == expected);
  x = expected = 1.0;
  while (x >= DBL_MIN * 4.0) {
    assert(asqrt(x) == expected);
    x *= 0.25;
    expected *= 0.5;
  }
  assert(asqrt(x) == expected);
  check(DBL_TRUE_MIN);
  check(DBL_MIN);
  check(DBL_MAX);

  const double infinity = strtod("inf", NULL);
  const double nan = strtod("nan", NULL);
  assert(asqrt(infinity) == infinity);
  assert(asqrt(-infinity) == -1.0);
  assert(asqrt(nan) != asqrt(nan));

  char text[32];
  snprintf(text, sizeof(text), "%a", asqrt(-0.0));
  assert(text[0] == '-');
  return 0;
}
