SHELL:=/bin/bash
cc = gcc
args = -g -Wall -Wextra -Werror -std=c99 -pedantic
bin = ../bin
PROGS = mafBlockSorter
dependencies = ../include/common.h

.PHONY: all clean test

all: $(foreach f,${PROGS}, ${bin}/$f)

${bin}/mafBlockSorter: src/mafBlockSorter.c ${dependencies}
	mkdir -p $(dir $@)
	${cc} ${args} $< -I ../include -I ../external -o $@.tmp 
	mv $@.tmp $@

test/mafBlockSorter: src/mafBlockSorter.c ${dependencies}
	mkdir -p $(dir $@)
	${cc} ${args} -O0 $< -I ../include -I ../external -o $@.tmp 
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

test: test/mafBlockSorter
	python2.6 src/test.mafBlockSorter.py --verbose && rm -rf test/
