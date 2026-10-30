#include <assert.h>
#include <stdio.h>

#include "extra.h"

int main(void) {
  uint32_t wide = 0;
  uint16_t narrow = 0;
  assert(parse_uint32(&wide, "70000") == 0 && wide == 70000);
  assert(parse_uint32(&wide, "4294967295") == 0 && wide == UINT32_MAX);
  assert(parse_uint32(&wide, "4294967296") != 0 && wide == UINT32_MAX);
  assert(parse_uint32(&wide, "12x") != 0 && wide == UINT32_MAX);
  assert(parse_uint32(&wide, "-0") != 0 && wide == UINT32_MAX);
  assert(parse_uint16(&narrow, "65535") == 0 && narrow == UINT16_MAX);
  assert(parse_uint16(&narrow, "65536") != 0 && narrow == UINT16_MAX);
  puts("numeric parser regressions passed");
  return 0;
}
