SHELL:=/bin/bash
cc = gcc
args = -Wall -Wextra -Werror -std=c99 -pedantic -I ../external -I ../include
bin = ../bin
PROGS = mafSorter
dependencies = $(wildcard ../include/common.*) $(wildcard ../include/sharedMaf.*)
objects = ../include/common.o ../include/sharedMaf.o ../external/CuTest.a src/buildVersion.o
testObjects = ./test/common.o ./test/sharedMaf.o ../external/CuTest.a test/buildVersion.o
sources = src/mafSorter.c

.PHONY: all clean test buildVersion

all: buildVersion $(foreach f,${PROGS}, ${bin}/$f)
buildVersion: src/buildVersion.c
src/buildVersion.c: ${sources} ${dependencies}
	@python ../include/createVersionSources.py

${bin}/mafSorter: src/mafSorter.c ${dependencies} ${objects}
	mkdir -p $(dir $@)
	${cc} ${args} -O3 $< ${objects} -o $@.tmp -lm
	mv $@.tmp $@

test/mafSorter: src/mafSorter.c ${dependencies} ${testObjects}
	mkdir -p $(dir $@)
	${cc} ${args} -g -O0 $< ${testObjects} -o $@.tmp -lm
	mv $@.tmp $@

%.o: %.c %.h
	${cc} -O3 -c ${args} $< -o $@.tmp -lm
	mv $@.tmp $@
test/%.o: ../include/%.c ../include/%.h
	mkdir -p $(dir $@)
	${cc} -g -O0 -c ${args} $< -o $@.tmp -lm
	mv $@.tmp $@
test/%.o: src/%.c src/%.h
	mkdir -p $(dir $@)
	${cc} -g -O0 -c ${args} $< -o $@.tmp -lm
	mv $@.tmp $@

clean:
	rm -rf $(foreach f,${PROGS}, ${bin}/$f) src/*.o test/

test: buildVersion test/mafSorter
	python2.7 src/test.mafSorter.py --verbose && rm -rf test/ && rmdir ./tempTestDir

../external/CuTest.a: ../external/CuTest.c ../external/CuTest.h
	${cc} -c ${args} $<
	ar rc CuTest.a CuTest.o
	ranlib CuTest.a
	rm -f CuTest.o
	mv CuTest.a $@
