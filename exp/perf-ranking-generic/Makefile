PERF_GEN_DIR := exp/perf-ranking-generic
EXP_DIRS += $(PERF_GEN_DIR)

.INTERMEDIATE: $(PERF_GEN_DIR)/perf-unranking-none-bin.tex \
               $(PERF_GEN_DIR)/perf-unranking-cache.tex \
               $(PERF_GEN_DIR)/perf-unranking-small-cache.tex

# -o 31: All orders
# -a 33: Default (1) + Generic (32)
# -c 3: None (1) + Bin (2)
$(PERF_GEN_DIR)/perf-unranking-none-bin.dat: bin/perf
	./$< -m 128 --k0 30 --k1 80 -c 3 -o 31 -a 33 -R 1 -p 1 -i 1024 \
		| sed 's/-default//g; s/-ucycles//g; s/-none//g' \
		> $@

# -c 20: Comb (4) + Acc (16)
$(PERF_GEN_DIR)/perf-unranking-cache.dat: bin/perf
	./$< -m 128 --k0 30 --k1 80 -c 4 -o 63 -a 33 -R 1 -p 1 -i 1024 \
		| sed 's/-default//g; s/-ucycles//g; s/-none//g' \
		> $@

# -c 40: Small Comb (8) + Small Acc (32)
$(PERF_GEN_DIR)/perf-unranking-small-cache.dat: bin/perf
	./$< -m 128 --k0 30 --k1 80 -c 40 -o 31 -a 33 -R 1 -p 1 -i 1024 \
		| sed 's/-default//g; s/-ucycles//g; s/-none//g' \
		> $@

$(PERF_GEN_DIR)/perf-unranking-%.tex: $(PERF_GEN_DIR)/perf-graph.tex \
								$(PERF_GEN_DIR)/perf-unranking-%.dat
	cp $< $@

perf-ranking-generic: $(PERF_GEN_DIR)/perf-unranking-none-bin.pdf \
				$(PERF_GEN_DIR)/perf-unranking-cache.pdf \
				$(PERF_GEN_DIR)/perf-unranking-small-cache.pdf
