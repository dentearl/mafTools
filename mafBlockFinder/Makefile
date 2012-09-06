SHELL:=/bin/bash
cc = gcc
args = -Wall -Wextra -std=c99 -pedantic -I ../external -I ../include
bin = ../bin
PROGS = mafBlockFinder
dependencies = $(wildcard ../include/common.*) $(wildcard ../include/sharedMaf.*)
objects = ../include/common.o ../include/sharedMaf.o ../external/CuTest.a
testObjects = test/common.o test/sharedMaf.o ../external/CuTest.a

.PHONY: all clean test

all: $(foreach f,${PROGS}, ${bin}/$f)

${bin}/mafBlockFinder: src/mafBlockFinder.c ${dependencies} ${objects}
	mkdir -p $(dir $@)
	${cc} ${args} -O3 $< ${objects} -o $@.tmp
	mv $@.tmp $@

test/mafBlockFinder: src/mafBlockFinder.c ${dependencies} ${testObjects}
	mkdir -p $(dir $@)
	${cc} ${args} -g -O0 $< ${testObjects} -o $@.tmp
	mv $@.tmp $@

%.o: %.c %.h
	${cc} -O3 -c ${args} $< -o $@.tmp
	mv $@.tmp $@
test/%.o: ../include/%.c ../include/%.h
	mkdir -p $(dir $@)
	${cc} -g -O0 -c ${args} $< -o $@.tmp
	mv $@.tmp $@
test/%.o: src/%.c src/%.h 
	mkdir -p $(dir $@)
	${cc} -c $< -o $@.tmp ${args} -g -O0
	mv $@.tmp $@

clean:
	rm -rf $(foreach f,${PROGS}, ${bin}/$f) src/*.o test/

test: test/mafBlockFinder
	python2.7 src/test.mafBlockFinder.py --verbose && rm -rf test/ && rmdir ./tempTestDir

../external/CuTest.a: ../external/CuTest.c ../external/CuTest.h
	${cc} -c ${args} $<
	ar rc CuTest.a CuTest.o
	ranlib CuTest.a
	rm -f CuTest.o
	mv CuTest.a $@
