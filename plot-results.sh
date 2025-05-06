#!/usr/bin/env sh

IT="${1:-$(make -pqrR | awk '/^IT/ { print $NF }')}"

SECURITY="128 192 256"
BACKENDS="bitint boost-fix"
ORDERS="colex gray rbo"
ALGORITHMS="default ps al ab ad"
CACHE_STRAT="bin comb scomb acc"

for LEVEL in $SECURITY ; do
  RANK_DATA_PATH="$TMPDIR/c-cycles-rank-all-m-$LEVEL-it-$IT.dat"
  UNRANK_DATA_PATH="$TMPDIR/c-cycles-unrank-all-m-$LEVEL-it-$IT.dat"
  rm -f $TMPDIR/c-*-cycles-test-m-$LEVEL-*.dat "$RANK_DATA_PATH" "$UNRANK_DATA_PATH"

  for IMPL in $BACKENDS ; do
    for ORD in $ORDERS ; do
      [ $ORD != "colex" ] && ALG_FIX="default" || ALG_FIX="$ALGORITHMS"
      for ALG in $ALG_FIX ; do
        for CACHE in $CACHE_STRAT ; do
          RAW_DATA_PATH="$TMPDIR/c-$IMPL-cycles-test-m-$LEVEL-it-$IT-o-$ORD-a-$ALG-c-$CACHE.dat"
          echo "$IMPL-$ORD-$ALG-$CACHE" | tee "$RAW_DATA_PATH"
          make --silent IT="$IT" ORDER="$ORD" ALG="$ALG" CACHE="$CACHE" \
            clean "$IMPL" "test-$LEVEL" >> "$RAW_DATA_PATH"
        done
      done
    done
  done

  paste -d" " $TMPDIR/c-*-cycles-test-m-$LEVEL-*.dat \
    | awk '
        NR == 1 { print "k", $0 }
        NR != 1 {
          printf "%6d", int($6)
          for (i = 21; i <= NF; i += 29) {
            printf "%14.2Lf", $i
          }
          printf "\n"
        }' \
    | column -t \
    > "$UNRANK_DATA_PATH"

  paste -d" " $TMPDIR/c-*-cycles-test-m-$LEVEL-*.dat \
    | awk '
        NR == 1 { print "k", $0 }
        NR != 1 {
          printf "%6d", int($6)
          for (i = 28; i <= NF; i += 29) {
            printf "%14.2Lf", $i
          }
          printf "\n"
        }' \
    | column -t \
    > "$RANK_DATA_PATH"

  for VAR in $UNRANK_DATA_PATH $RANK_DATA_PATH ; do
    gnuplot -e "
      stats '/dev/stdin' skip 1 nooutput;
      max_col = STATS_columns;
      set linetype cycle max_col;
      set key autotitle columnhead;
      set logscale y;
      set terminal png size 2560, 1440;
      plot for [i = 2:max_col] '/dev/stdin' using 1:i
        with lines linewidth 3 smooth mcsplines
    " < "$VAR" > "${VAR/dat/png}"
  done
done
