#!/usr/bin/env sh

IT="${1:-$(make -pqrR | awk '/^IT/ { print $NF }')}"

for LEVEL in 256 ; do
  RAW_DATA_PATH="$TMPDIR/bitint-time-all-versions-it-$IT-m-$LEVEL.dat"
  rm -f $TMPDIR/bitint-test-*.dat "$RAW_DATA_PATH"

  for ALG in colex enup emk ; do
    for CACHE in none bin comb ; do
    (
      echo "$ALG-$CACHE"
      make --silent IT=$IT ALG=$ALG CACHE=$CACHE clean test-$LEVEL
    ) >| "$TMPDIR/bitint-test-$LEVEL-$ALG-$CACHE.dat"
    done
  done

  paste -d" " $TMPDIR/bitint-test-$LEVEL-*.dat \
    | awk '
        NR == 1 { print "k", $0 }
        NR > 1 {
          printf "%3d ", $6
          for (i = 19; i <= NF; i += 24) {
            printf "%14.2lf", $i
          }
          printf "\n"
        }' \
    | column -t \
    > "$RAW_DATA_PATH"
  gnuplot -e "
    set linetype cycle 30;
    set key autotitle columnhead;
    set logscale y;
    set terminal png size 2560, 1440;
    plot for [i = 2:10] '/dev/stdin' using 1:i
      with lines linewidth 3 smooth mcsplines
  " < "$RAW_DATA_PATH" > "${RAW_DATA_PATH/dat/png}"
done
