include ../include.mk
SHELL:=/bin/bash
cc = gcc
sonLib = ../../sonLib/lib
bin = ../bin
PROGS = mafTransitiveClosure
dependencies = $(wildcard ../include/common.*) $(wildcard ../include/sharedMaf.*) $(wildcard ${sonLib}/*) ${sonLib}/sonLib.a ${sonLib}/stPinchesAndCacti.a src/allTests.c
objects := ../include/common.o ../include/sharedMaf.o ${sonLib}/sonLib.a ../external/CuTest.a ${sonLib}/stPinchesAndCacti.a src/test.mafTransitiveClosure.o
testObjects := test/sharedMaf.o test/common.o ${sonLib}/sonLib.a ../external/CuTest.a ${sonLib}/stPinchesAndCacti.a src/test.mafTransitiveClosure.o
testFlags = -O0 -g -Wall -Werror --pedantic -lm

.PHONY: all clean test

all: $(foreach f,${PROGS}, ${bin}/$f)

${bin}/mafTransitiveClosure: src/mafTransitiveClosure.c ${dependencies} ${objects}
	mkdir -p $(dir $@)
	${cc} $< src/allTests.c ${objects} -o $@.tmp ${cflags}
	mv $@.tmp $@

test/mafTransitiveClosure: src/mafTransitiveClosure.c ${dependencies} ${testObjects}
	mkdir -p $(dir $@)
	${cc} $< src/allTests.c ${testObjects} -o $@.tmp ${testFlags}
	mv $@.tmp $@
%.o: %.c %.h
	${cc} -c $< -o $@.tmp ${cflags}
	mv $@.tmp $@
test/%.o: ../include/%.c ../include/%.h
	mkdir -p $(dir $@)
	${cc} -c $< -o $@.tmp ${testFlags}
	mv $@.tmp $@

clean:
	rm -rf $(foreach f,${PROGS}, ${bin}/$f) src/*.o test/ 

test: test/mafTransitiveClosure
	$< --test && python2.7 src/test.mafTransitiveClosure.py --verbose && rm -rf ./test/ && rmdir ./tempTestDir

../external/CuTest.a: ../external/CuTest.c ../external/CuTest.h
	${cc} -c $< ${cflags}
	ar rc CuTest.a CuTest.o
	ranlib CuTest.a
	rm -f CuTest.o
	mv CuTest.a $@
