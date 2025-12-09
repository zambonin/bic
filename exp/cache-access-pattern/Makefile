ACCESS_DIR := exp/cache-access-pattern
EXP_DIRS += $(ACCESS_DIR)

CACHE_TYPES := bin comb scomb acc sacc

.INTERMEDIATE: $(patsubst %,$(ACCESS_DIR)/%-cache-pattern.tex,$(CACHE_TYPES)) \
               $(patsubst %,$(ACCESS_DIR)/%-cache-pattern.gp,$(CACHE_TYPES)) \
               $(patsubst %,$(ACCESS_DIR)/%-cache-pattern.cairo.gp,$(CACHE_TYPES))

TARGET_PNG := $(patsubst %,$(ACCESS_DIR)/%-cache-pattern.png,$(CACHE_TYPES))
TARGET_PDF := $(TARGET_PNG:.png=.pdf)

define CACHE_RULE
# the algorithm configuration for heatmap plots makes a very slight difference
$(ACCESS_DIR)/$(1)-cache-pattern.dat: bin/cache-access
	./$$< -n 65 -k 53 -d 20 -c $(1) -a ad -i16384 > $$@

$(ACCESS_DIR)/$(1)-cache-pattern.gp: $(ACCESS_DIR)/template-$(2).gp
	cp $$< $$@

$(ACCESS_DIR)/$(1)-cache-pattern.cairo.gp: $(ACCESS_DIR)/template-$(2).cairo.gp
	cp $$< $$@
endef

$(eval $(call CACHE_RULE,bin,heatmap))
$(eval $(call CACHE_RULE,comb,heatmap))
$(eval $(call CACHE_RULE,scomb,heatmap))
$(eval $(call CACHE_RULE,acc,scatter3d))
$(eval $(call CACHE_RULE,sacc,scatter3d))

cache-access-pattern: $(TARGET_PNG) $(TARGET_PDF)
