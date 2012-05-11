SHELL:=/bin/bash
cc = gcc
args = -g -Wall -Wextra -Werror -std=c99 -pedantic
bin = ../bin
PROGS = mafBlockSorter
dependencies = ../include/common.h ../include/sharedMaf.h

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

clean:
	rm -rf $(foreach f,${PROGS}, ${bin}/$f) src/*.o test/

test: test/mafBlockSorter
	python2.7 src/test.mafBlockSorter.py --verbose && rm -rf test/
