SHELL:=/bin/bash
cc = gcc
args = -g -Wall -Wextra -Werror -std=c99 -pedantic
bin = ../bin
PROGS = mafBlockDuplicateFilter
dependencies = ../include/common.h ../include/sharedMaf.h

.PHONY: all clean test

all: $(foreach f,${PROGS}, ${bin}/$f)

${bin}/mafBlockDuplicateFilter: src/mafBlockDuplicateFilter.c ${dependencies}
	mkdir -p $(dir $@)
	${cc} ${args} $< -I ../include -I ../external -o $@.tmp 
	mv $@.tmp $@

test/mafBlockDuplicateFilter: src/mafBlockDuplicateFilter.c ${dependencies}
	mkdir -p $(dir $@)
	${cc} ${args} -O0 $< -I ../include -I ../external -o $@.tmp 
	mv $@.tmp $@

clean:
	rm -rf $(foreach f,${PROGS}, ${bin}/$f) src/*.o test/

test: test/mafBlockDuplicateFilter
	python2.7 src/test.mafBlockDuplicateFilter.py --verbose && rm -rf test/
