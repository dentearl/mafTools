SHELL:=/bin/bash -e
export SHELLOPTS=pipefail
.PHONY: all clean test

cflags = -std=c99 -O3 -Wextra -Wall -Werror -pedantic -I ../external/

objects = common.o sharedMaf.o ../external/CuTest.a

all: ${objects}

clean:
	rm -f allTests *.o *.pyc

allTests: allTests.c test.sharedMaf.h test.sharedMaf.c ${objects}
	gcc -g -O0 ${cflags} $^ -o $@.tmp
	mv $@.tmp $@

%.o: %.c %.h
	gcc -c ${cflags} $< -o $@.tmp
	mv $@.tmp $@

test: allTests
	./allTests && python2.7 test.sharedMaf.py --verbose && rm -f ./allTests

../external/CuTest.a: ../external/CuTest.c ../external/CuTest.h
	${cc} -c ${args} $<
	ar rc CuTest.a *.o
	ranlib CuTest.a
	rm -f *.o
	mv CuTest.a $@
