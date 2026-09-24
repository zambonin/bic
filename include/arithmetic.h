#ifndef BIC_ARITHMETIC_H
#define BIC_ARITHMETIC_H

#include "common.h"
#include <stdbool.h>

static inline void bic_reduce_ratio(uint32_t *numerator,
                                    uint32_t *denominator) {
  uint32_t a = *numerator, b = *denominator;
  while (b != 0) {
    const uint32_t remainder = a % b;
    a = b;
    b = remainder;
  }
  *numerator /= a;
  *denominator /= a;
}

static inline uintx bic_muldiv_exact(uintx value, uint32_t numerator,
                                     uint32_t denominator) {
  bic_reduce_ratio(&numerator, &denominator);
  return (value / denominator) * numerator;
}

static inline intx bic_signed_muldiv_exact(const intx value,
                                           const uint32_t numerator,
                                           const uint32_t denominator) {
  const bool negative = value < 0;
  const uintx magnitude = negative ? (uintx)(-value) : (uintx)value;
  const uintx result = bic_muldiv_exact(magnitude, numerator, denominator);
  return negative ? -(intx)result : (intx)result;
}

#endif
