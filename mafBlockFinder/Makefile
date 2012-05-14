SHELL:=/bin/bash
cc = gcc
args = -Wall -Wextra -std=c99 -pedantic -I ../external -I ../include
bin = ../bin
PROGS = mafBlockFinder
dependencies = $(wildcard ../include/common.*) $(wildcard ../include/sharedMaf.*)
objects = ../include/common.o ../include/sharedMaf.o

.PHONY: all clean test

all: $(foreach f,${PROGS}, ${bin}/$f)

${bin}/mafBlockFinder: src/mafBlockFinder.c ${dependencies} ${objects}
	mkdir -p $(dir $@)
	${cc} ${args} -O3 $< ${objects} -o $@.tmp
	mv $@.tmp $@

test/mafBlockFinder: src/mafBlockFinder.c ${dependencies} ${objects}
	mkdir -p $(dir $@)
	${cc} ${args} -g -O0 $< ${objects} -o $@.tmp
	mv $@.tmp $@

%.o: %.c %.h
	${cc} -c ${args} $< -o $@.tmp
	mv $@.tmp $@

clean:
	rm -rf $(foreach f,${PROGS}, ${bin}/$f) src/*.o test/

test: test/mafBlockFinder
	python2.7 src/test.mafBlockFinder.py --verbose && rm -rf test/
