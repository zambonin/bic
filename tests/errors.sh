#!/bin/sh
set -eu

scratch=$(mktemp -d)
trap 'rm -rf "$scratch"' EXIT HUP INT TERM

check_status() {
  expected=$1
  shift
  status=0
  "$@" > "$scratch/out" 2> "$scratch/err" || status=$?
  test "$status" -eq "$expected"
  test ! -s "$scratch/err"
}

for program in cli enum cache-access cache-hits; do
  check_status 1 "./bin/$program" --unknown
  check_status 1 "./bin/$program" -nwrong -n0
  check_status 1 "./bin/$program" -n0 -k2 -d1 -ocolex extra
  if [ "$program" = enum ]; then
    check_status 2 "./bin/$program" -n4294967295 -k2 -d4294967295 -ocolex
  else
    check_status 2 "./bin/$program" -n4294967295 -k2 -d4294967295 -ocolex -cnone
  fi
done
check_status 1 ./bin/cache-access -n0 -k2 -d1 -cnone --operation=invalid
check_status 0 ./bin/cache-access -n0 -k2 -d1 -cnone --operation=unrank -q1
check_status 0 ./bin/cli -n0 -k2 -d1 -i1
check_status 1 ./bin/cli -n0 -k2 -d1 -i0
check_status 1 ./bin/cli -n0 -k2 -d1 -v256
