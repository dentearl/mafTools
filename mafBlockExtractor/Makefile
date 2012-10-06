SHELL:=/bin/bash
cc = gcc
args = -Wall -Wextra -std=c99 -pedantic -I ../external -I ../include
bin = ../bin
PROGS = mafBlockExtractor
dependencies = $(wildcard ../include/common.*) $(wildcard ../include/sharedMaf.*) src/mafBlockExtractor.h
API = ../include/common.o ../include/sharedMaf.o ../external/CuTest.a src/mafBlockExtractorAPI.o src/buildVersion.o
testAPI = test/sharedMaf.o ../external/CuTest.a test/common.o test/mafBlockExtractorAPI.o test/buildVersion.o
testObjects := test/test.mafBlockExtractor.o

.PHONY: all clean test buildVersion

all: buildVersion $(foreach f,${PROGS}, ${bin}/$f)
buildVersion:
	@python ../include/createVersionSources.py

${bin}/mafBlockExtractor: src/mafBlockExtractor.c ${dependencies} ${API}
	mkdir -p $(dir $@)
	${cc} ${args} -O3 $< ${API} -o $@.tmp -lm
	mv $@.tmp $@

test/mafBlockExtractor: src/mafBlockExtractor.c ${dependencies} ${testAPI}
	mkdir -p $(dir $@)
	${cc} ${args} -g -O0 $< ${testAPI} -o $@.tmp -lm
	mv $@.tmp $@

%.o: %.c %.h
	${cc} -O3 -c ${args} $< -o $@.tmp -lm
	mv $@.tmp $@

test/%.o: ../include/%.c ../include/%.h
	mkdir -p $(dir $@)
	${cc} -c $< -o $@.tmp ${args} -g -O0 -lm
	mv $@.tmp $@
test/%.o: src/%.c src/%.h 
	mkdir -p $(dir $@)
	${cc} -c $< -o $@.tmp ${args} -g -O0 -lm
	mv $@.tmp $@

clean:
	rm -rf $(foreach f,${PROGS}, ${bin}/$f) src/*.o test/ temTestDir/

test: buildVersion test/allTests test/mafBlockExtractor
	./test/allTests && python2.7 src/test.mafBlockExtractor.py --verbose && rmdir ./tempTestDir && rm -rf ./test/

test/allTests: src/allTests.c ${testObjects} ${testAPI}
	mkdir -p $(dir $@)
	${cc} $^ -o $@.tmp ${args} -g -O0 -lm
	mv $@.tmp $@

test/test.mafBlockExtractor.o: src/test.mafBlockExtractor.c src/test.mafBlockExtractor.h ${testAPI}
	mkdir -p $(dir $@)
	${cc} -c $< -o $@.tmp ${args} -g -O0 -lm
	mv $@.tmp $@

../external/CuTest.a: ../external/CuTest.c ../external/CuTest.h
	${cc} -c ${args} $<
	ar rc CuTest.a CuTest.o
	ranlib CuTest.a
	rm -f CuTest.o
	mv CuTest.a $@
