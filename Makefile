CFLAGS = -Wall -Wextra -pedantic -Ofast -g -std=gnu23 -march=native -mtune=native
LDFLAGS = -lm
SRC = src/unrank.c
OUT = $(basename $(SRC))
IT = 128

all: $(OUT)

debug: CFLAGS += -g
debug: all

clean:
	$(RM) $(OUT)

%.test: $(OUT)
	@./$< -m $(subst -, -k ,$*) -i $(IT) -c comb

test-256: $(foreach TVAL,$(shell seq 30 80),256-$(TVAL).test)
test-512: $(foreach TVAL,$(shell seq 60 140),512-$(TVAL).test)
test: test-256 test-512
