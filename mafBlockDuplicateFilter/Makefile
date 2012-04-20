SHELL:=/bin/bash
cc = gcc
args = -g -Wall -Wextra -Werror -std=c99 -pedantic
bin = ../bin
PROGS = mafBlockDuplicateFilter
INCL = ../include/common.o
dependencies = ${INCL} ../include/common.h ../include/common.c

.PHONY: all clean test

all: $(foreach f,${PROGS}, ${bin}/$f)

${bin}/mafBlockDuplicateFilter: src/mafBlockDuplicateFilter.c ${dependencies}
	mkdir -p $(dir $@)
	${cc} ${args} $< ${INCL} -I ../include -o $@.tmp 
	mv $@.tmp $@

test/mafBlockDuplicateFilter: src/mafBlockDuplicateFilter.c ${dependencies}
	mkdir -p $(dir $@)
	${cc} ${args} -O0 $< ${INCL} -I ../include -o $@.tmp 
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

test: test/mafBlockDuplicateFilter
	python2.6 src/test.mafBlockDuplicateFilter.py --verbose && rm -rf test/
