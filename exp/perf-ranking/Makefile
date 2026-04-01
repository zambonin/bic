PERF_DIR := exp/perf-ranking
EXP_DIRS += $(PERF_DIR)

.INTERMEDIATE: $(PERF_DIR)/perf-unranking-none-bin.tex \
               $(PERF_DIR)/perf-unranking-cache.tex \
               $(PERF_DIR)/perf-unranking-small-cache.tex

$(PERF_DIR)/perf-unranking-none-bin.dat: bin/perf
	./$< -m 128 --k0 30 --k1 80 -c 3 -o 5 -p 1 -i 1024 \
		| sed 's/-default//g; s/-ucycles//g; s/-none//g' \
		> $@

$(PERF_DIR)/perf-unranking-cache.dat: bin/perf
	./$< -m 128 --k0 30 --k1 80 -c 20 -o 5 -p 1 -i 1024 \
		| sed 's/-default//g; s/-ucycles//g; s/-none//g' \
		> $@

$(PERF_DIR)/perf-unranking-small-cache.dat: bin/perf
	./$< -m 128 --k0 30 --k1 80 -c 40 -o 5 -p 1 -i 1024 \
		| sed 's/-default//g; s/-ucycles//g; s/-none//g' \
		> $@

$(PERF_DIR)/perf-unranking-%.tex: $(PERF_DIR)/perf-graph.tex \
								$(PERF_DIR)/perf-unranking-%.dat
	cp $< $@

perf-unranking: $(PERF_DIR)/perf-unranking-none-bin.pdf \
				$(PERF_DIR)/perf-unranking-cache.pdf \
				$(PERF_DIR)/perf-unranking-small-cache.pdf
