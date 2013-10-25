SHELL:=/bin/bash
cc = gcc
args = -Wall -Wextra -Werror -std=c99 -pedantic -I ../external -I ../include
bin = ../bin
PROGS = mafStrander
dependencies = $(wildcard ../include/common.*) $(wildcard ../include/sharedMaf.*)
objects = ../include/common.o ../include/sharedMaf.o ../external/CuTest.a src/buildVersion.o
testObjects := test/sharedMaf.o test/common.o ../external/CuTest.a  test/buildVersion.o
sources = src/mafStrander.c

.PHONY: all clean test buildVersion

all: buildVersion $(foreach f,${PROGS}, ${bin}/$f)
buildVersion: src/buildVersion.c
src/buildVersion.c: ${sources} ${dependencies}
	@python ../include/createVersionSources.py

${bin}/mafStrander: src/mafStrander.c ${dependencies} ${objects}
	mkdir -p $(dir $@)
	${cc} ${args} -O3 $< ${objects} -o $@.tmp
	mv $@.tmp $@

test/mafStrander: src/mafStrander.c ${dependencies} ${testObjects}
	mkdir -p $(dir $@)
	${cc} ${args} -g -O0 $< ${testObjects} -o $@.tmp
	mv $@.tmp $@

%.o: %.c %.h
	${cc} -c ${args} $< -o $@.tmp
	mv $@.tmp $@
test/%.o: ../include/%.c ../include/%.h
	mkdir -p $(dir $@)
	${cc} -g -O0 -c ${args} $< -o $@.tmp
	mv $@.tmp $@
test/%.o: src/%.c src/%.h
	mkdir -p $(dir $@)
	${cc} -g -O0 -c ${args} $< -o $@.tmp
	mv $@.tmp $@

clean:
	rm -rf $(foreach f,${PROGS}, ${bin}/$f) src/*.o test/

test: buildVersion test/mafStrander
	python2.7 src/test.mafStrander.py --verbose && rm -rf test/ && rmdir ./tempTestDir

../external/CuTest.a: ../external/CuTest.c ../external/CuTest.h
	${cc} -c ${args} $<
	ar rc CuTest.a CuTest.o
	ranlib CuTest.a
	rm -f CuTest.o
	mv CuTest.a $@
