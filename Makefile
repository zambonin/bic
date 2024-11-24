CFLAGS = -Wall -Wextra -pedantic -Ofast -g -std=gnu23
LDFLAGS = -lm
SRC = src/unrank.c
OUT = $(basename $(SRC))
IT = 128

all: $(OUT)

clean:
	$(RM) $(OUT)
