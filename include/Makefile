SHELL:=/bin/bash -e
export SHELLOPTS=pipefail
.PHONY: all clean test

all:

clean:
	rm allTests

allTests: allTests.c test.sharedMaf.h ../external/CuTest.c sharedMaf.h common.h
	gcc -g -O0 -std=c99 -Wextra -Wall -Werror -pedantic -I ../external/ $^ -o $@.tmp
	mv $@.tmp $@

test: allTests
	./allTests && python2.7 test.sharedMaf.py --verbose && rm ./allTests
