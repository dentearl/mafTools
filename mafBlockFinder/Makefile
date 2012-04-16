SHELL:=/bin/bash
cc = gcc
args = -Wall -Wextra -std=c99 -pedantic
bin = ../bin
PROGS = mafBlockFinder
includes = ../include/common.o
dependencies = ${includes} ../include/common.h ../include/common.c

.PHONY: all clean

all: $(foreach f,${PROGS}, ${bin}/$f)

${bin}/mafBlockFinder: src/mafBlockFinder.c ${dependencies}
	mkdir -p $(dir $@)
	${cc} ${args} $< -I ../include -o $@.tmp ${includes}
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
	rm -rf $(foreach f,${PROGS}, ${bin}/$f) src/*.o

test: ${bin}/mafBlockFinder
	python2.6 src/test.mafBlockFinder.py -v
