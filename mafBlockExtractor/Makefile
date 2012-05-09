SHELL:=/bin/bash
cc = gcc
args = -g -Wall -Wextra -std=c99 -pedantic
bin = ../bin
PROGS = mafBlockExtractor
dependencies = ../include/common.h ../include/sharedMaf.h

.PHONY: all clean test

all: $(foreach f,${PROGS}, ${bin}/$f)

${bin}/mafBlockExtractor: src/mafBlockExtractor.c ${dependencies}
	mkdir -p $(dir $@)
	${cc} ${args} $< -I ../include -I ../external -o $@.tmp 
	mv $@.tmp $@

test/mafBlockExtractor: src/mafBlockExtractor.c ${dependencies}
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

test: test/mafBlockExtractor
	python2.7 src/test.mafBlockExtractor.py --verbose && rm -rf ./test/
