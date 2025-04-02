CC = clang++
CFLAGS = -Wall -Wextra -pedantic -O3 -ffast-math -std=gnu23 -march=native -mtune=native -x c++
LDFLAGS = -lm
SRC = src/unrank.c
OUT = $(basename $(SRC))

RANGE = $(shell seq 30 80)
IT = 128
ALG = colex
CACHE = none
PRINT = none

default:

bitint: CC = clang
bitint: CFLAGS := $(patsubst -x c++,-DBITINT,$(CFLAGS))
bitint: $(OUT)

boost-fix: CFLAGS := $(patsubst -std=gnu23,-DBOOST_FIX_INT,$(CFLAGS))
boost-fix: $(OUT)

boost-arb: CFLAGS := $(patsubst -std=gnu23,-DBOOST_ARB_INT,$(CFLAGS))
boost-arb: $(OUT)

mpz: CFLAGS := $(patsubst -std=gnu23,-DBOOST_MPZ_INT,$(CFLAGS))
mpz: LDFLAGS += -lgmp
mpz: $(OUT)

tom: CFLAGS := $(patsubst -std=gnu23,-DBOOST_TOM_INT,$(CFLAGS))
tom: LDFLAGS += -ltommath
tom: $(OUT)

%.test: $(OUT)
	./$< -m $(subst -, -k ,$*) -a $(ALG) -i $(IT) -c $(CACHE) -p $(PRINT)

test-128: $(foreach K,$(RANGE),128-$(K).test)
test-192: $(foreach K,$(RANGE),192-$(K).test)
test-256: $(foreach K,$(RANGE),256-$(K).test)
test: test-128 test-192 test-256

%.leak: $(OUT) /usr/bin/valgrind
	valgrind --quiet --exit-on-first-error=yes --leak-check=full \
		--errors-for-leak-kinds=all --show-leak-kinds=all --error-exitcode=1 \
		./$< -m $(subst -, -k ,$*) -a $(ALG) -i $(IT) -c $(CACHE) -p $(PRINT) \
		2>/dev/null

leak-128: $(foreach K,$(RANGE),128-$(K).leak)
leak-192: $(foreach K,$(RANGE),192-$(K).leak)
leak-256: $(foreach K,$(RANGE),256-$(K).leak)
leak: leak-128 leak-192 leak-256

%.stats.png: $(OUT) /usr/bin/gnuplot
	./$< -m $(subst -, -k ,$*) -a $(ALG) -i $(IT) -c $(CACHE) -p $(PRINT) \
		| gnuplot -e "set terminal png size 2560, 1440; \
			plot '/dev/stdin' matrix with image notitle" > $@

stats-128: $(foreach K,$(RANGE),128-$(K).stats.png)
stats-192: $(foreach K,$(RANGE),192-$(K).stats.png)
stats-256: $(foreach K,$(RANGE),256-$(K).stats.png)
stats: stats-128 stats-192 stats-256

clean-stats:
	$(RM) $(wildcard *.stats.png)

clean:
	$(RM) $(OUT)
