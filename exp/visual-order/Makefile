VISUAL_DIR := exp/visual-order
EXP_DIRS += $(VISUAL_DIR)

VISUAL_ORDERS := colex gray rbo spiral bous

define VISUAL_RULE
$(VISUAL_DIR)/visual-order-comparison-$(1).dat: bin/enum
	./$$< -n6 -k5 -d4 -o $(1) \
		| perl $(VISUAL_DIR)/grow-and-transpose-enum.pl \
		> $$@
endef

$(foreach o,$(VISUAL_ORDERS),$(eval $(call VISUAL_RULE,$(o))))

$(VISUAL_DIR)/visual-order-comparison.pdf: \
	$(VISUAL_DIR)/visual-order-comparison.tex \
	$(patsubst %,$(VISUAL_DIR)/visual-order-comparison-%.dat,$(VISUAL_ORDERS))

visual-order: $(VISUAL_DIR)/visual-order-comparison.pdf
