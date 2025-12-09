set terminal cairolatex png standalone size 20cm, 25cm
set output ARG2

set datafile separator ","
set key off

set palette defined ( \
  0 '#440154', 0.25 '#3b528b', 0.5 '#21918c', 0.75 '#5ec962', 1 '#fde725' \
)

set xrange [0:64]
set yrange [0:64]
set grid front linetype 1 linecolor rgb "white" linewidth 0.5

tt = "Cache size reductions (ratio) for Gaussian/packing optimizations ($n = \\frac{kd}{2}$)"

set multiplot \
    layout 4, 3 \
    margins 0.08, 0.93, 0.06, 0.97 \
    spacing 0.01, 0.01 \
    title tt

array metrics[3] = [5, 6, 7]
array names[3] = [ \
    "Binomial cache", "Comp. count cache", "Accumulated sums cache", \
]

do for [stddev = 1:4] {
  do for [j = 1:3] {
    if (stddev == 4) {
      set format x "%g"
      set xlabel sprintf ("$k$\n(%s)", names[j])
    } else {
      set format x ""
      unset xlabel
    }

    if (j == 1) {
      set format y "%g"
      set ylabel sprintf("$d$\n(%d$\\sigma$)", stddev) rotate by 0 offset -0.5,0
    } else {
      set format y ""
      unset ylabel
    }

    if (j < 3) {
      unset colorbox
    } else {
      set colorbox
    }

    plot ARG1 \
        using ($4 == stddev ? $2 : 1/0): \
             ($3): \
             ($2-0.5):($2+0.5): \
             ($3-0.5):($3+0.5): \
             metrics[j] \
        with boxxyerror fs solid 1.0 linecolor palette
  }
}

unset multiplot
