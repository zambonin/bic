#!/usr/bin/env sh

IT="${1:-$(make -pqrR | awk '/^IT/ { print $NF }')}"

SECURITY="128 256"
BACKENDS="bitint boost-fix boost-arb"
ALGORITHMS="colex colexpart colexbs colexdbcs gray rbo"
CACHE_STRAT="bin comb acc"

for LEVEL in $SECURITY ; do
  ALL_DATA_PATH="$TMPDIR/c-cycles-all-m-$LEVEL-it-$IT.dat"
  rm -f $TMPDIR/c-*-cycles-test-m-$LEVEL-*.dat "$ALL_DATA_PATH"

  for IMPL in $BACKENDS ; do
    for ALG in $ALGORITHMS ; do
      for CACHE in $CACHE_STRAT ; do
        RAW_DATA_PATH="$TMPDIR/c-$IMPL-cycles-test-m-$LEVEL-it-$IT-a-$ALG-c-$CACHE.dat"
        echo "$IMPL-$ALG-$CACHE" | tee "$RAW_DATA_PATH"
        make --silent IT="$IT" ALG="$ALG" CACHE="$CACHE" \
          clean "$IMPL" "test-$LEVEL" >> "$RAW_DATA_PATH"
      done
    done
  done

  paste -d" " $TMPDIR/c-*-cycles-test-m-$LEVEL-*.dat \
    | awk '
        NR == 1 { print "k", $0 }
        NR != 1 {
          printf "%6d", int($6)
          for (i = 24; i <= NF; i += 24) {
            printf "%14.2Lf", $i
          }
          printf "\n"
        }' \
    | column -t \
    > "$ALL_DATA_PATH"

  gnuplot -e "
    stats '/dev/stdin' skip 1 nooutput;
    max_col = STATS_columns;
    set linetype cycle max_col;
    set key autotitle columnhead;
    set logscale y;
    set terminal png size 2560, 1440;
    plot for [i = 2:max_col] '/dev/stdin' using 1:i
      with lines linewidth 3 smooth mcsplines
  " < "$ALL_DATA_PATH" > "${ALL_DATA_PATH/dat/png}"
done
