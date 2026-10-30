#!/bin/sh
set -eu
root=$(pwd)
scratch=$(mktemp -d)
trap 'rm -rf "$scratch"' EXIT HUP INT TERM
for selection in colex:default colex:ps colex:ab rbo:default spiral:default; do
    order=${selection%:*}
    algorithm=${selection#*:}
    for cache in none bin comb scomb acc sacc; do
        : > "$scratch/partition"
        for operation in rank unrank round-trip; do
            ./bin/cache-access -n6 -k6 -d2 -r42 -i32 -q1 "-o$order" "-a$algorithm" "-c$cache" "--operation=$operation" |
                awk -v expected="$operation" '
                    { n = split($0, fields, /, */)
                      for (i = 1; i <= n; i++) {
                        split(fields[i], pair, / *= */); gsub(/^ +| +$/, "", pair[1]); values[pair[1]] = pair[2]
                      }
                      if (values["operation"] != expected) exit 1
                      print values["lookups"], values["cache_hits"]
                    }' >> "$scratch/partition"
        done
        awk 'NR == 1 {a=$1;b=$2} NR == 2 {a+=$1;b+=$2}
             NR == 3 {if ($1 != a || $2 != b) exit 1} END {if (NR != 3) exit 1}' "$scratch/partition"
    done
done
mkdir "$scratch/bin"
cat > "$scratch/bin/cli" <<'EOF'
#!/bin/sh
set -eu
program=${0##*/}
if [ "${BIC_DRIVER_FAIL:-}" = "$program" ]; then exit 3; fi
printf '%s' "$program" >> "$BIC_DRIVER_LOG"
for arg do printf '\t%s' "$arg" >> "$BIC_DRIVER_LOG"; done
printf '\n' >> "$BIC_DRIVER_LOG"
for arg do
    case "$arg" in
        -n*) n=${arg#-n} ;; -k*) k=${arg#-k} ;; -d*) d=${arg#-d} ;;
        -o*) o=${arg#-o} ;; -a*) a=${arg#-a} ;; -c*) c=${arg#-c} ;;
        -v*) v=${arg#-v} ;;
    esac
done
if [ "$program" = cli ]; then
    printf 'n=%s, k=%s, d=%s, o=%s, a=%s, c=%s, unrank_avg_cyc=1000\n' "$n" "$k" "$d" "$o" "$a" "$c"
else
    printf 'hit_rate=0.5, size=1000, stddev=%s\n' "$v"
fi
EOF
cp "$scratch/bin/cli" "$scratch/bin/cache-access"
printf '#!/bin/sh\nexit 0\n' > "$scratch/bin/make"
chmod +x "$scratch/bin/cli" "$scratch/bin/cache-access" "$scratch/bin/make"
export BIC_DRIVER_LOG="$scratch/commands"
export PATH="$scratch/bin:$PATH"
cd "$scratch"
sh "$root/run-iwoca.sh" -s42 -i32 > report
awk -F '\t' '
    NR % 2 == 1 {
      if ($1 != "cli") exit 1
      delete args; for (i=2; i<=NF; i++) args[$i] = 1
      next
    }
    { if ($1 != "cache-access") exit 1
      for (i=2; i<=NF; i++) {
        if ($i == "-q1" || $i == "--operation=unrank") continue
        if (!($i in args)) exit 1
        delete args[$i]
      }
      for (arg in args) exit 1
    }
    END { if (NR != 240) exit 1 }' commands
for strategy in colex colex-ps colex-ab rbo spiral; do
    head -n1 report | grep -q "$strategy"_hit
done
for program in cli cache-access; do
    if BIC_DRIVER_FAIL="$program" sh "$root/run-iwoca.sh" -s42 -i32 > failed-report 2>/dev/null; then
        echo "driver ignored $program failure" >&2
        exit 1
    fi
    test ! -s failed-report
done
echo 'Cache access: 30 partitions, 120 driver pairs and failure propagation passed'
