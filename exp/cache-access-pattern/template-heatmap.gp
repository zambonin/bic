set terminal pngcairo size 1280, 720

unset key

set logscale cb

set autoscale fix

set xlabel "n"
set ylabel "k" rotate by 0
set cblabel "Cache access rate" offset 2,0

set cbtics format "10^{%T}"

plot "/dev/stdin" matrix using 2:1:($3 == -1 ? NaN : $3) with image
