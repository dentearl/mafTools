SHELL:=/bin/bash
cc = gcc
args = -g -Wall -Wextra -std=c99 -pedantic
bin = ../bin
PROGS = mafBlockExtractor
INCL = src/common.c
dependencies = ${INCL} src/common.h

.PHONY: all clean

all: $(foreach f,${PROGS}, ${bin}/$f)

${bin}/mafBlockExtractor: src/mafBlockExtractor.c ${dependencies}
	mkdir -p $(dir $@)
	${cc} ${args} $< ${INCL} -o $@.tmp 
	mv $@.tmp $@

${bin}/%: src/%.c 
	mkdir -p $(dir $@)
	${cc} ${args} $^ -o $@.tmp 
	mv $@.tmp $@

clean:
	rm -rf $(foreach f,${PROGS}, ${bin}/$f) src/*.o
