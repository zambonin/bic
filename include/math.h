#ifndef MATH_H
#define MATH_H

#include "common.h"

extern long double logl(long double);

extern double sqrt(double);

long double lg(const uintx u);

long double lg_bic(const uint16_t n, const uint16_t k, const uint16_t d);

// §4 of 10.1007/s13389-021-00264-9
void mingen(const uint16_t m, uint16_t *n, const uint16_t k, uint16_t *d);

// §2.3 of 10.1007/s13389-021-00264-9
void minver(const uint16_t m, uint16_t *n, const uint16_t k, uint16_t *d);

// §6.1 of 10.1007/978-3-642-14764-7_6
uintx inner_bin(const uint16_t n, const uint16_t k, const uint16_t d);

uintx bin(const uint16_t n, const uint16_t k, const uint16_t d);

// remark 60 of 10.24033/asens.136
uintx inner_bic_with_sums(const uint16_t n, const uint16_t k, const uint16_t d,
                          intx *partial_sums, math_func bin_impl);

uintx inner_bic(const uint16_t n, const uint16_t k, const uint16_t d);

uintx bic(const uint16_t n, const uint16_t k, const uint16_t d);

uintx *inner_acc(const uint16_t n, const uint16_t k, const uint16_t d);

uintx *acc(const uint16_t n, const uint16_t k, const uint16_t d);

// proposition 3 of 10.1007/s13389-021-00264-9
uintx bic_acc(const uint16_t n, const uint16_t k, const uint16_t d,
              const uint16_t l);

uintx random_rank(const uint16_t n, const uint16_t k, const uint16_t d);

#endif
