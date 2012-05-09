SHELL:=/bin/bash
cc = gcc
args = -Wall -Wextra -std=c99 -pedantic
bin = ../bin
PROGS = mafBlockFinder
dependencies = ../include/common.h ../include/sharedMaf.h

.PHONY: all clean test

all: $(foreach f,${PROGS}, ${bin}/$f)

${bin}/mafBlockFinder: src/mafBlockFinder.c ${dependencies}
	mkdir -p $(dir $@)
	${cc} ${args} $< -I ../include -I ../external -o $@.tmp ${includes}
	mv $@.tmp $@

test/mafBlockFinder: src/mafBlockFinder.c ${dependencies}
	mkdir -p $(dir $@)
	${cc} ${args} -O0 $< -I ../include -I ../external -o $@.tmp ${includes}
	mv $@.tmp $@

${bin}/%: src/%.c 
	mkdir -p $(dir $@)
	${cc} ${args} $^ -o $@.tmp 
	mv $@.tmp $@

../include/common.o: ../include/common.c ../include/common.h
	mkdir -p $(dir $@)
	${cc} ${args} $< -c -o $@.tmp 
	mv $@.tmp $@

clean:
	rm -rf $(foreach f,${PROGS}, ${bin}/$f) src/*.o test/

test: test/mafBlockFinder
	python2.6 src/test.mafBlockFinder.py --verbose && rm -rf test/
