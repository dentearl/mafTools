SHELL:=/bin/bash
cc = gcc
sonLib = ../../sonLib/lib
args = -Wall -Wextra -Werror -std=c99 -pedantic -I ../external -I ../include -I ${sonLib} -lm
bin = ../bin
PROGS = mafTransitiveClosure
dependencies = $(wildcard ../include/common.*) $(wildcard ../include/sharedMaf.*) $(wildcard ${sonLib}/*) ${sonLib}/sonLib.a ${sonLib}/stPinchesAndCacti.a src/allTests.c
objects := ../include/common.o ../include/sharedMaf.o ${sonLib}/sonLib.a ../external/CuTest.a ${sonLib}/stPinchesAndCacti.a src/test.mafTransitiveClosure.o

.PHONY: all clean test

all: $(foreach f,${PROGS}, ${bin}/$f)

${bin}/mafTransitiveClosure: src/mafTransitiveClosure.c ${dependencies} ${objects}
	mkdir -p $(dir $@)
	${cc} ${args} -O3 $< src/allTests.c ${objects} -o $@.tmp 
	mv $@.tmp $@

test/mafTransitiveClosure: src/mafTransitiveClosure.c ${dependencies} ${objects}
	mkdir -p $(dir $@)
	${cc} ${args} -g -O0 $< src/allTests.c ${objects} -o $@.tmp 
	mv $@.tmp $@

%.o: %.c %.h
	${cc} -c ${args} $< -o $@.tmp
	mv $@.tmp $@

clean:
	rm -rf $(foreach f,${PROGS}, ${bin}/$f) src/*.o test/ 

test: test/mafTransitiveClosure
	$< --test && python2.7 src/test.mafTransitiveClosure.py --verbose && rm -rf ./test/ && rmdir ./tempTestDir

../external/CuTest.a: ../external/CuTest.c ../external/CuTest.h
	${cc} -c ${args} $<
	ar rc CuTest.a *.o
	ranlib CuTest.a
	rm -f *.o
	mv CuTest.a $@
