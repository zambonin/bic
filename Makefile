CC = g++
CFLAGS = -Wall -Wextra -pedantic -O3 -march=native -mtune=native -Iinclude
LDFLAGS = -lm
SRC = $(wildcard src/*.c)
OBJ = $(SRC:.c=.o)
OUT = src/main

RANGE = $(shell seq 30 80)
IT = 128
ORDER = colex
ALG = default
CACHE = none
PRINT = none
STRATEGY = gen
PARAMS = -o $(ORDER) -a $(ALG) -i $(IT) -c $(CACHE) -p $(PRINT) -s $(STRATEGY)

default:

$(OUT): $(OBJ)

bitint: CC = gcc
bitint: CFLAGS += -std=gnu23 -DBITINT
bitint: $(OUT)

boost-fix: CFLAGS += -x c++ -DBOOST_FIX_INT -std=c++20
boost-fix: $(OUT)

boost-arb: CFLAGS += -x c++ -DBOOST_ARB_INT -std=c++20
boost-arb: $(OUT)

mpz: CFLAGS += -x c++ -DBOOST_MPZ_INT -std=c++20
mpz: LDFLAGS += -lgmp
mpz: $(OUT)

tom: CFLAGS += -x c++ -DBOOST_TOM_INT -std=c++20
tom: LDFLAGS += -ltommath
tom: $(OUT)

%.test: $(OUT)
	./$< -m $(subst -, -k ,$*) $(PARAMS)

test-128: $(foreach K,$(RANGE),128-$(K).test)
test-192: $(foreach K,$(RANGE),192-$(K).test)
test-256: $(foreach K,$(RANGE),256-$(K).test)
test: test-128 test-192 test-256

%.leak: CC = gcc
%.leak: CFLAGS += -std=gnu23 -DBITINT -mno-avx512f
%.leak: IT = 1
%.leak: $(OUT) /usr/bin/valgrind
	valgrind --quiet --exit-on-first-error=yes --leak-check=full \
		--errors-for-leak-kinds=all --show-leak-kinds=all --error-exitcode=1 \
		./$< -m $(subst -, -k ,$*) $(PARAMS) 2>/dev/null

leak-128: $(foreach K,$(RANGE),128-$(K).leak)
leak-192: $(foreach K,$(RANGE),192-$(K).leak)
leak-256: $(foreach K,$(RANGE),256-$(K).leak)
leak: leak-128 leak-192 leak-256

%.stats.png: $(OUT) /usr/bin/gnuplot
	./$< -m $(subst -, -k ,$*) \
		| ( IFS=" ," read _ _ N _ _ K _ _ D _; \
			./$< -m $(subst -, -k ,$*) $(PARAMS) \
				| gnuplot -c plot-cache-access.gp $$N $$K $$D ) \
		> $@
stats-128: $(foreach K,$(RANGE),128-$(K).stats.png)
stats-192: $(foreach K,$(RANGE),192-$(K).stats.png)
stats-256: $(foreach K,$(RANGE),256-$(K).stats.png)
stats: stats-128 stats-192 stats-256

clean-stats:
	$(RM) $(wildcard *.stats.png)

clean:
	$(RM) $(OUT) $(OBJ)
