CC = clang++
CFLAGS = -Wall -Wextra -pedantic -O3 -ffast-math -std=gnu23 -march=native -mtune=native -x c++
LDFLAGS = -lm
SRC = src/unrank.c
OUT = $(basename $(SRC))
IT = 128
ALG = colex
CACHE = none

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
	./$< -m $(subst -, -k ,$*) -a $(ALG) -i $(IT) -c $(CACHE)

test-256: $(foreach K,$(shell seq 30 80),256-$(K).test)
test-512: $(foreach K,$(shell seq 60 140),512-$(K).test)
test: test-256 test-512

%.leak: $(OUT) /usr/bin/valgrind
	valgrind --quiet --exit-on-first-error=yes --leak-check=full \
		--errors-for-leak-kinds=all --show-leak-kinds=all --error-exitcode=1 \
		./$< -m $(subst -, -k ,$*) -a $(ALG) -i $(IT) -c $(CACHE) \
		2>/dev/null

leak-256: $(foreach K,$(shell seq 30 80),256-$(K).leak)
leak-512: $(foreach K,$(shell seq 60 140),512-$(K).leak)
leak: leak-256 leak-512

clean:
	$(RM) $(OUT)
