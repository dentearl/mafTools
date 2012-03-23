SHELL:=/bin/bash
cc = gcc
args = -g -Wall -Wextra -std=c99 -pedantic
bin = ../bin
PROGS = mafBlockExtractor
INCL = src/common.o
dependencies = ${INCL} src/common.h src/common.c

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

src/common.o: src/common.c src/common.h
	mkdir -p $(dir $@)
	${cc} ${args} $< -c -o $@.tmp 
	mv $@.tmp $@

clean:
	rm -rf $(foreach f,${PROGS}, ${bin}/$f) src/*.o

test: ${bin}/mafBlockExtractor
	python2.6 src/test.mafBlockExtractor.py -v
