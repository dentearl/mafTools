SHELL:=/bin/bash -e
export SHELLOPTS=pipefail
.PHONY: all clean test

cflags = -std=c99 -Wextra -Wall -Werror -pedantic -I ../external/

all:

clean:
	rm allTests *.o *.pyc

allTests: allTests.c test.sharedMaf.h ../external/CuTest.c common.o sharedMaf.o 
	gcc -g -O0 ${cflags} $^ -o $@.tmp
	mv $@.tmp $@

%.o: %.c %.h
	gcc -c ${cflags} $< -o $@.tmp
	mv $@.tmp $@

test: allTests
	./allTests && python2.7 test.sharedMaf.py --verbose && rm ./allTests
