set terminal cairolatex png standalone size 18cm, 12cm
set output ARG2

unset key

set logscale cb

set autoscale fix

set xlabel "$n$"
set ylabel "$k$" rotate by 0
set cblabel "Cache access rate" offset 2,0

set cbtics format "$10^{%T}$"

plot ARG1 matrix using 2:1:($3 == -1 ? NaN : $3) with image
