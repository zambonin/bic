CC = g++
CFLAGS = -Wall -Wextra -pedantic -O3 -march=native -mtune=native -Iinclude -D_XOPEN_SOURCE=500
SRC = $(filter-out src/extra.c, $(wildcard src/*.c))
TARGETS = $(basename $(wildcard bin/*.c))
OBJ = $(SRC:.c=.o)
LIB = libbic.a
LIBEXTRA = libbicext.a

PRIMES = include/primes.h

BACKEND ?= bitint
INTWIDTH ?= 512

ifeq ($(BACKEND), bitint)
    CFLAGS += -x c -std=c23 -DBITINT=$(INTWIDTH)
else ifeq ($(BACKEND), boost-fix)
    CFLAGS += -DBOOST_FIX_INT -std=c++20
else ifeq ($(BACKEND), boost-arb)
    CFLAGS += -DBOOST_ARB_INT -std=c++20
else ifeq ($(BACKEND), mpz)
    CFLAGS += -DBOOST_MPZ_INT -std=c++20
    LDLIBS += -lgmp
else ifeq ($(BACKEND), tom)
    CFLAGS += -DBOOST_TOM_INT -std=c++20
    LDLIBS += -ltommath
endif

default:

src/math.o: $(PRIMES)

$(LIB): $(OBJ)
	$(AR) rcs $@ $^

$(LIBEXTRA): $(OBJ) src/extra.o
	$(AR) rcs $@ $^

$(PRIMES): src/gen-primes.awk /usr/bin/awk
	awk -v width=$(INTWIDTH) -f $< > $@

$(TARGETS): %: %.o $(LIBEXTRA)

test: bin/test
	./$< -i 128 -r 0

%.leak: CFLAGS += -mno-avx512f -g
%.leak: % /usr/bin/valgrind
	valgrind --quiet --exit-on-first-error=yes --leak-check=full \
		--errors-for-leak-kinds=all --show-leak-kinds=all --error-exitcode=1 \
		./$<

clean:
	$(RM) $(OUT) $(wildcard src/*.o) $(wildcard bin/*.o) $(LIB) $(TARGETS) \
		$(LIBEXTRA) $(PRIMES)

include Makefile.figures
include Makefile.regression
