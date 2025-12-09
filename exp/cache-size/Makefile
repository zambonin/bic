CACHE_SIZE_DIR := exp/cache-size
EXP_DIRS += $(CACHE_SIZE_DIR)

.INTERMEDIATE: $(CACHE_SIZE_DIR)/grid-cache-size.tex

$(CACHE_SIZE_DIR)/grid-cache-size.dat: bin/cache-size
	./$< -i8192 -v4 | sed '1d; s/\t/,/g' | uniq | sort -n -k2,3 -t, > $@

cache-size: $(CACHE_SIZE_DIR)/grid-cache-size.png \
			$(CACHE_SIZE_DIR)/grid-cache-size.pdf
