CC = clang
CFLAGS = -Wall -Wextra -pedantic -O3 -ffast-math -std=gnu23 -march=native -mtune=native
LDFLAGS = -lm
SRC = src/unrank.c
OUT = $(basename $(SRC))
IT = 128
ALG = colex
CACHE = none

all: $(OUT)

debug: CFLAGS += -g
debug: all

clean:
	$(RM) $(OUT)

%.test: $(OUT)
	@./$< -m $(subst -, -k ,$*) -a $(ALG) -i $(IT) -c $(CACHE)

test-256: $(foreach TVAL,$(shell seq 30 80),256-$(TVAL).test)
test-512: $(foreach TVAL,$(shell seq 60 140),512-$(TVAL).test)
test: test-256 test-512
