PERF_GEN_DIR := exp/perf-ranking-generic
EXP_DIRS += $(PERF_GEN_DIR)

.INTERMEDIATE: $(PERF_GEN_DIR)/perf-unranking-none-bin.tex \
               $(PERF_GEN_DIR)/perf-unranking-cache.tex \
               $(PERF_GEN_DIR)/perf-unranking-small-cache.tex

# -o 31 is first five orders
# -a 33 is default and generic algorithms
# -c 3 is no cache and bin cache
$(PERF_GEN_DIR)/perf-unranking-none-bin.dat: bin/perf
	./$< -m 128 --k0 30 --k1 80 -c 3 -o 31 -a 33 -R 1 -p 1 -i 1024 \
		| sed 's/-default//g; s/-ucycles//g; s/-none//g' \
		> $@

# -o 63 is all six orders
# -c 4 is comb cache
$(PERF_GEN_DIR)/perf-unranking-cache.dat: bin/perf
	./$< -m 128 --k0 30 --k1 80 -c 4 -o 63 -a 33 -R 1 -p 1 -i 1024 \
		| sed 's/-default//g; s/-ucycles//g; s/-none//g' \
		> $@

# -c 40 is scomb and sacc caches
$(PERF_GEN_DIR)/perf-unranking-small-cache.dat: bin/perf
	./$< -m 128 --k0 30 --k1 80 -c 40 -o 31 -a 33 -R 1 -p 1 -i 1024 \
		| sed 's/-default//g; s/-ucycles//g; s/-none//g' \
		> $@

$(PERF_GEN_DIR)/perf-unranking-%.tex: exp/perf-ranking/perf-graph.tex \
								$(PERF_GEN_DIR)/perf-unranking-%.dat
	cp $< $@

perf-ranking-generic: $(PERF_GEN_DIR)/perf-unranking-none-bin.pdf \
				$(PERF_GEN_DIR)/perf-unranking-cache.pdf \
				$(PERF_GEN_DIR)/perf-unranking-small-cache.pdf
