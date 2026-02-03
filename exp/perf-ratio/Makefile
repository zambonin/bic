RATIO_DIR := exp/perf-ratio
EXP_DIRS += $(RATIO_DIR)

$(RATIO_DIR)/perf-ratio.dat: bin/perf_ratio
	./$< -i 1000 > $@

$(RATIO_DIR)/perf-ratio.pdf: $(RATIO_DIR)/perf-ratio.tex $(RATIO_DIR)/perf-ratio.dat

perf-ratio: $(RATIO_DIR)/perf-ratio.pdf
