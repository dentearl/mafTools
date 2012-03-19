SHELL:=/bin/bash
cc = gcc
args = -Wall -Wextra -std=c99 -pedantic
bin = ../bin
PROGS = mafBlockFinder
includes = src/common.c
dependencies = ${includes} src/common.h src/common.c

.PHONY: all clean

all: $(foreach f,${PROGS}, ${bin}/$f)

${bin}/mafBlockFinder: src/mafBlockFinder.c ${dependencies}
	mkdir -p $(dir $@)
	${cc} ${args} $< ${includes} -o $@.tmp 
	mv $@.tmp $@

${bin}/%: src/%.c 
	mkdir -p $(dir $@)
	${cc} ${args} $^ -o $@.tmp 
	mv $@.tmp $@

clean:
	rm -rf $(foreach f,${PROGS}, ${bin}/$f) src/*.o

test:
	echo "No tests."
