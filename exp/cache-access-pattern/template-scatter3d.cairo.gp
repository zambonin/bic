set terminal cairolatex png standalone size 18cm, 12cm
set output ARG2

unset key

set logscale cb
set colorbox horizontal user origin screen 0.09, 0.95 size 0.86, 0.02

set xlabel "$n$"
set ylabel "$d$"
set zlabel "$k$"
set cblabel "Cache access rate" offset 30,0

set xtics offset 0, -0.5
set ytics offset 0, -0.5
set ztics offset -0.5, 0
set cbtics format "$10^{%T}$"

set ticslevel 0

set view 60, 240, 1.15
set origin 0.02, 0

splot "< sort -n -k4 " . ARG1 using 1:3:2:4 with points \
    palette pointsize 0.5 pointtype 15
