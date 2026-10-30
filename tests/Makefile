TEST_PROGRAMS = $(basename $(wildcard tests/*.c))
TEST_OBJECTS = $(addsuffix .o,$(TEST_PROGRAMS))
CORE_TESTS = tests/asqrt tests/domain tests/admission tests/reflection tests/rbo tests/accumulated tests/parameters tests/errors
LEGACY_TESTS = $(filter-out $(CORE_TESTS) tests/roundtrip tests/seq tests/rte,$(TEST_PROGRAMS))

$(TEST_PROGRAMS): %: %.o $(LIBEXTRA)

.PHONY: check test check-core check-legacy check-random check-seq test-seq check-rte check-tools clean-tests
check test: check-core check-legacy check-random check-tools
check-core: $(CORE_TESTS)
	@set -e; for program in $(CORE_TESTS); do ./$$program; done
check-legacy: $(LEGACY_TESTS)
	@set -e; for program in $(LEGACY_TESTS); do ./$$program; done
check-random: tests/roundtrip
	./$< -i4 -r0
check-seq test-seq: tests/seq
	./$< -i1 -r0
check-rte: tests/rte
	./$< -i1 -r0
check-tools: bin/perf bin/cache-access bin/cli bin/enum bin/cache-hits
	sh tests/errors.sh
	sh tests/cli.sh
	sh tests/cache-access.sh
clean: clean-tests
clean-tests:
	$(RM) $(TEST_PROGRAMS) $(TEST_OBJECTS)

tests/errors: LDLIBS += -Wl,--wrap=malloc -Wl,--wrap=calloc
