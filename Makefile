CC = g++
CFLAGS = -Wall -Wextra -pedantic -O3 -march=native -mtune=native -Iinclude -D_XOPEN_SOURCE=500
SRC = $(filter-out src/extra.c, $(wildcard src/*.c src/order/*.c src/cache/*.c src/generic/*.c))
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

$(PRIMES): src/gen-primes.awk /usr/bin/awk
	awk -v width=$(INTWIDTH) -f $< > $@

$(LIB): $(OBJ)
	$(RM) $@
	$(AR) rcs $@ $^

$(LIBEXTRA): $(OBJ) src/extra.o
	$(RM) $@
	$(AR) rcs $@ $^

$(TARGETS): %: %.o $(LIBEXTRA)

clean:
	$(RM) $(OUT) $(wildcard src/*.o src/order/*.o src/cache/*.o src/generic/*.o) $(wildcard bin/*.o) \
		$(LIB) $(LIBEXTRA) $(TARGETS) $(PRIMES)

include Makefile.figures

include tests/Makefile
