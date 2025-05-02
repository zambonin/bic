#!/usr/bin/env gnuplot --persist

set terminal png size 2560, 1440

n = ARG1
k = ARG2
d = ARG3

set xrange[1:k - 1]
set yrange[0:n]

set linetype 1 lw 2 lc rgb "red"
set linetype 2 lw 2 lc rgb "sienna1"
set linetype 3 lw 2 lc rgb "orange"
set linetype 4 lw 2 lc rgb "goldenrod"
set linetype 5 lw 2 lc rgb "yellow"
set linetype 6 lw 2 lc rgb "goldenrod"
set linetype 7 lw 2 lc rgb "orange"
set linetype 8 lw 2 lc rgb "sienna1"
set linetype 9 lw 2 lc rgb "red"

variance_factor = d * (d + 2) / 12.0
mean(j) = j * n / k
stddev(j) = sqrt(j * variance_factor * (k - j) / k)
sigma(j, c) = mean(j) + c * stddev(j)

level = 3

plot '/dev/stdin' matrix with image notitle, \
    for [i=-level:level] sigma(x, i) notitle with points
