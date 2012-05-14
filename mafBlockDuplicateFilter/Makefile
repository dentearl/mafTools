SHELL:=/bin/bash
cc = gcc
args = -Wall -Wextra -Werror -std=c99 -pedantic -I ../external -I ../include
bin = ../bin
PROGS = mafBlockDuplicateFilter
dependencies = $(wildcard ../include/common.*) $(wildcard ../include/sharedMaf.*)
objects = ../include/common.o ../include/sharedMaf.o

.PHONY: all clean test

all: $(foreach f,${PROGS}, ${bin}/$f)

${bin}/mafBlockDuplicateFilter: src/mafBlockDuplicateFilter.c ${dependencies}
	mkdir -p $(dir $@)
	${cc} ${args} -O3 $< ${objects} -o $@.tmp 
	mv $@.tmp $@

test/mafBlockDuplicateFilter: src/mafBlockDuplicateFilter.c ${dependencies}
	mkdir -p $(dir $@)
	${cc} ${args} -g -O0 $< ${objects} -o $@.tmp 
	mv $@.tmp $@

%.o: %.c %.h
	${cc} -c ${args} $< -o $@.tmp
	mv $@.tmp $@

clean:
	rm -rf $(foreach f,${PROGS}, ${bin}/$f) src/*.o test/

test: test/mafBlockDuplicateFilter
	python2.7 src/test.mafBlockDuplicateFilter.py --verbose && rm -rf test/
