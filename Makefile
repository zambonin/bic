CC = g++
CFLAGS = -Wall -Wextra -pedantic -O3 -march=native -mtune=native -Iinclude -D_XOPEN_SOURCE=500
SRC = $(wildcard src/*.c)
TARGETS = $(basename $(wildcard bin/*.c))
OBJ = $(SRC:.c=.o)
LIB = libbic.a

VALGRIND_PATH ?= /usr/bin/valgrind
PLOT_CACHE_ACCESS_SCRIPT = plot-from-dhat.py

BACKEND ?= bitint
INTWIDTH ?= 512

RANGE = $(shell seq 30 80)
IT = 128
ORDER = colex
ALG = default
CACHE = none
STRATEGY = gen
PARAMS = -o $(ORDER) -a $(ALG) -i $(IT) -c $(CACHE) -s $(STRATEGY)

ifeq ($(BACKEND), bitint)
    CC = gcc
    CFLAGS += -std=c23 -DBITINT=$(INTWIDTH)
else ifeq ($(BACKEND), boost-fix)
    CFLAGS += -x c++ -DBOOST_FIX_INT -std=c++20
else ifeq ($(BACKEND), boost-arb)
    CFLAGS += -x c++ -DBOOST_ARB_INT -std=c++20
else ifeq ($(BACKEND), mpz)
    CFLAGS += -x c++ -DBOOST_MPZ_INT -std=c++20
    LDLIBS += -lgmp
else ifeq ($(BACKEND), tom)
    CFLAGS += -x c++ -DBOOST_TOM_INT -std=c++20
    LDLIBS += -ltommath
endif

.PHONY: all clean

default:

$(LIB): $(OBJ)
	$(AR) rcs $@ $^

$(TARGETS): LDFLAGS += -L.
$(TARGETS): LDLIBS += -lbic
$(TARGETS): %: %.o $(LIB)

run-tests: bin/test
	./$< -i $(IT) -s 0

%.cli: bin/cli
	./$< -m $(subst -, -k ,$*) $(PARAMS)

cli-128: $(foreach K,$(RANGE),128-$(K).cli)
cli-192: $(foreach K,$(RANGE),192-$(K).cli)
cli-256: $(foreach K,$(RANGE),256-$(K).cli)
run-cli-tests: cli-128 cli-192 cli-256

leak: CC = gcc
leak: CFLAGS += -std=c23 -DBITINT=$(INTWIDTH) -mno-avx512f
leak: IT = 1
leak: bin/test $(VALGRIND_PATH)
	$(VALGRIND_PATH) --quiet --exit-on-first-error=yes --leak-check=full \
		--errors-for-leak-kinds=all --show-leak-kinds=all --error-exitcode=1 \
		./$< -i $(IT)

%.stats.png: CC = gcc
%.stats.png: CFLAGS += -std=c23 -DBITINT=$(INTWIDTH) -mno-avx512f -DDHAT
%.stats.png: bin/cli $(VALGRIND_PATH) $(PLOT_CACHE_ACCESS_SCRIPT)
	 $(VALGRIND_PATH) --quiet --tool=dhat --dhat-out-file=/dev/stdout \
		./$< -m $(subst -, -k ,$*) $(PARAMS) 2>/dev/null \
		| uv run --script $(PLOT_CACHE_ACCESS_SCRIPT) \
		> $@
stats-128: $(foreach K,$(RANGE),128-$(K).stats.png)
stats-192: $(foreach K,$(RANGE),192-$(K).stats.png)
stats-256: $(foreach K,$(RANGE),256-$(K).stats.png)
stats: stats-128 stats-192 stats-256

clean-stats:
	$(RM) $(wildcard *.stats.png)

clean:
	$(RM) $(OUT) $(wildcard src/*.o) $(wildcard bin/*.o) $(LIB) $(TARGETS)
