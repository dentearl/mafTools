SHELL=/bin/bash
.PHONY: all clean test

cc = gcc
args = -std=c99 -O3 -Wextra -Wall -Werror -pedantic -I ../external/

objects = common.o sharedMaf.o ../external/CuTest.a
testObjects := test/sharedMaf.o test/common.o ../external/CuTest.a

all: ${objects}

clean:
	rm -f allTests *.o *.pyc

allTests: allTests.c test.sharedMaf.h test.sharedMaf.c ${testObjects}
	${cc} -g -O0 ${args} $^ -o $@.tmp
	mv $@.tmp $@

%.o: %.c %.h
	${cc} -O3 -c ${args} $< -o $@.tmp
	mv $@.tmp $@

sharedMaf.o: sharedMaf.c sharedMaf.h
	${cc} -O3 -c ${args} sharedMaf.c -o $@.tmp -lm
	mv $@.tmp $@

test/%.o: %.c %.h
	mkdir -p $(dir $@)
	${cc} -g -O0 -c ${args} $< -o $*.tmp
	mv $*.tmp $@

test/sharedMaf.o: sharedMaf.c sharedMaf.h
	${cc} -g -O0 -c ${args} sharedMaf.c -o $@.tmp -lm
	mv $@.tmp $@

test: allTests
	./allTests && python2.7 test.sharedMaf.py --verbose && rm -rf ./allTests ./test && rmdir ./tempTestDir

../external/CuTest.a: ../external/CuTest.c ../external/CuTest.h
	${cc} -c ${args} $<
	ar rc CuTest.a CuTest.o
	ranlib CuTest.a
	rm -f CuTest.o
	mv CuTest.a $@
