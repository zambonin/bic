#!/usr/bin/env sh

pprint_table() {
  # shellcheck disable=SC1010,SC2016
  mlr \
    --opprint \
    clean-whitespace \
    then put '
      $alg = $o . "-" . $a;
      if ($stddev != 0) {
        $cache = $c . $stddev
      } else {
        $cache = $c
      }
      $size = $size / 1000;
      $unrank_avg_cyc = $unrank_avg_cyc / 1000' \
    then cut -o -f n,k,d,cache,size,alg,unrank_avg_cyc,hit_rate \
    then rename hit_rate,unrank_hit_rate \
    then reshape -r 'unrank_avg_cyc|unrank_hit_rate' -o metric,value \
    then put '$col = $alg . "_" . $metric' \
    then cut -x -f alg,metric \
    then reshape -s col,value \
    then unsparsify \
    then rename -r '^(.*)_unrank_avg_cyc$,\1' \
    then rename -r '^(.*)_unrank_hit_rate$,\1_hit' \
    then rename -r '^(.*)-default_hit$,\1_hit' \
    then rename -r '^(.*)-default$,\1' \
    then reorder -f n,k,d,cache,size
}

measure() {
  timing=$(./bin/cli -r"$SEED" -i"$ITERATIONS" \
      -n"$1" -k"$2" -d"$3" \
      -o"$4" -a"$5" -c"$6" -v"$7") || return
  accesses=$(./bin/cache-access -r"$SEED" -i"$ITERATIONS" \
      -q1 -n"$1" -k"$2" -d"$3" \
      -o"$4" -a"$5" -c"$6" -v"$7" --operation=unrank) || return
  printf '%s,%s\n' "$timing" "$accesses"
}

measure_order() {
  measure "$1" "$2" "$3" "$4" "$5" none 0 || return
  measure "$1" "$2" "$3" "$4" "$5" bin 0 || return
  measure "$1" "$2" "$3" "$4" "$5" scomb 2 || return
  measure "$1" "$2" "$3" "$4" "$5" scomb 4 || return
  measure "$1" "$2" "$3" "$4" "$5" comb 0 || return
  measure "$1" "$2" "$3" "$4" "$5" sacc 2 || return
  measure "$1" "$2" "$3" "$4" "$5" sacc 4 || return
  measure "$1" "$2" "$3" "$4" "$5" acc 0 || return
}

measure_case() {
  measure_order "$1" "$2" "$3" colex default || return
  measure_order "$1" "$2" "$3" colex ps || return
  measure_order "$1" "$2" "$3" colex ab || return
  measure_order "$1" "$2" "$3" rbo default || return
  measure_order "$1" "$2" "$3" spiral default || return
}

SEED=
ITERATIONS=
OPTION=

while getopts s:i:-: OPTION ; do
  if [ "$OPTION" = "-" ] ; then
    OPTION="${OPTARG%%=*}"
    OPTARG="${OPTARG#"$OPTION"}"
    OPTARG="${OPTARG#=}"
  fi
  case "$OPTION" in
    s | seed ) SEED="$OPTARG" ;;
    i | iterations ) ITERATIONS="$OPTARG" ;;
    ??* ) exit 1 ;;
    \? ) exit 2 ;;
  esac
done
shift $((OPTIND - 1))

if [ -z "$SEED" ] || [ -z "$ITERATIONS" ] || [ "$#" -ne 0 ] ; then
  cat <<-EOF
Usage: sh ${0##*/} [options]
  -s,  --seed=SEED
  -i,  --iterations=ITERATIONS
EOF
  exit 1
fi

make --silent clean || exit 1
make --silent CC="${CC:-clang}" INTWIDTH="${INTWIDTH:-256}" \
  BACKEND="${BACKEND:-bitint}" bin/cli bin/cache-access || exit 1

rows=$(
  measure_case 384 25 86 &&
  measure_case 66 67 9 &&
  measure_case 39 136 4
) || exit 1
printf '%s\n' "$rows" | pprint_table
