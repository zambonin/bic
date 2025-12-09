OUTPUT_DIR := results
BACKENDS := bitint boost-fix boost-arb mpz tom
RANGE ?= 8 	# commits behind

COMMITS := $(shell git rev-list --abbrev-commit --reverse HEAD \
				| tail --lines=$(RANGE) \
				| awk '{printf "%02d-%s\n", NR, $$1}')
LOGS_SEQ := $(foreach c,$(COMMITS),$(foreach b,$(BACKENDS),\
				$(OUTPUT_DIR)/$c-$b-seq-log.txt))
LOGS_TEST := $(foreach c,$(COMMITS),$(foreach b,$(BACKENDS),\
				$(OUTPUT_DIR)/$c-$b-test-log.txt))
LOGS := $(LOGS_SEQ) $(LOGS_TEST)
LEAKS := $(LOGS_TEST:log.txt=leak.txt)

default:

test-correctness: $(LOGS)

test-memory-leaks: $(LEAKS)

$(OUTPUT_DIR):
	mkdir -p $@

define RULE
$(OUTPUT_DIR)/$1-$2: | $(OUTPUT_DIR)
	git worktree add --detach --quiet $$@ $3

$(OUTPUT_DIR)/$1-$2/libbic.a: $(OUTPUT_DIR)/$1-$2
	$$(MAKE) -C $$< --no-print-directory libbic.a BACKEND=$2

$(OUTPUT_DIR)/$1-$2/bin/seq: $(OUTPUT_DIR)/$1-$2/libbic.a
	$$(MAKE) -C $$(dir $$<) --no-print-directory bin/seq BACKEND=$2

$(OUTPUT_DIR)/$1-$2/bin/test: $(OUTPUT_DIR)/$1-$2/libbic.a
	$$(MAKE) -C $$(dir $$<) --no-print-directory bin/test BACKEND=$2

$(OUTPUT_DIR)/$1-$2-seq-log.txt: $(OUTPUT_DIR)/$1-$2/bin/seq
	./$$< -i 8 -r 0 > $$@

$(OUTPUT_DIR)/$1-$2-test-log.txt: $(OUTPUT_DIR)/$1-$2/bin/test
	./$$< -i 8 -r 0 > $$@

$(OUTPUT_DIR)/$1-$2-test-leak.txt: $(OUTPUT_DIR)/$1-$2/bin/test
	valgrind --quiet --exit-on-first-error=yes --leak-check=full \
		--errors-for-leak-kinds=all --show-leak-kinds=all --error-exitcode=1 \
		./$$< -i 1 > /dev/null
endef

$(foreach c,$(COMMITS),$(foreach b,$(BACKENDS),\
	$(eval $(call RULE,$c,$b,$(word 2,$(subst -, ,$c))))))

clean-binaries:
	$(RM) -r $(OUTPUT_DIR)
	git worktree prune

clean-results:
	$(RM) $(wildcard $(OUTPUT_DIR)/*.txt)
