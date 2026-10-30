#!/bin/sh
set -eu
for option in o a c R; do
    if ./bin/perf -m8 --k0 4 --k1 4 -i1 "-$option"32768 >/dev/null 2>&1; then
        echo "invalid benchmark option accepted: $option" >&2
        exit 1
    fi
done
./bin/perf -m8 --k0 4 --k1 4 -i1 -r0 -o4 -a3 -R1 -c1 |
    awk 'NR == 1 { if (NF != 3 || $2 != "rbo-unrank:default-none-ucycles" || $3 != "rbo-rank:default-none-rcycles") exit 1 }
         END { if (NR != 2 || NF != 3) exit 1 }'
echo 'Benchmark argument checks passed'
