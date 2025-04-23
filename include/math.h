#ifndef MATH_H
#define MATH_H

#include <math.h>

#include "common.h"
#include "cache.h"

extern long double logl(long double);

extern double log(double);

long double logl2(const uintx u);

long double lg_bic(const uint16_t n, const uint16_t k, const uint16_t d);

void mingen(const uint16_t m, uint16_t *n, const uint16_t k, uint16_t *d);

void minver(const uint16_t m, uint16_t *n, const uint16_t k, uint16_t *d);

uintx bin_uiui(const uint16_t n, const uint16_t k, const uint16_t d);

uintx bin(const uint16_t n, const uint16_t k, const uint16_t d);

uintx inner_bic_with_sums(const uint16_t n, const uint16_t k, const uint16_t d,
                          intx *partial_sums, math_func bin_impl);

uintx inner_bic(const uint16_t n, const uint16_t k, const uint16_t d);

uintx bic(const uint16_t n, const uint16_t k, const uint16_t d);

uintx *inner_acc(const uint16_t n, const uint16_t k, const uint16_t d);

uintx *acc(const uint16_t n, const uint16_t k, const uint16_t d);

uintx bic_acc(const uint16_t n, const uint16_t k, const uint16_t d,
              const uint16_t l);

uintx random_rank(const uint16_t n, const uint16_t k, const uint16_t d);

#endif
