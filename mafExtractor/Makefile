SHELL:=/bin/bash
cc = gcc
args = -Wall -Wextra -std=c99 -pedantic -I ../external -I ../include
bin = ../bin
PROGS = mafExtractor
dependencies = $(wildcard ../include/common.*) $(wildcard ../include/sharedMaf.*) src/mafExtractor.h
API = ../include/common.o ../include/sharedMaf.o ../external/CuTest.a src/mafExtractorAPI.o src/buildVersion.o
testAPI = test/sharedMaf.o ../external/CuTest.a test/common.o test/mafExtractorAPI.o test/buildVersion.o
testObjects := test/test.mafExtractor.o

.PHONY: all clean test buildVersion

all: buildVersion $(foreach f,${PROGS}, ${bin}/$f)
buildVersion:
	@python ../include/createVersionSources.py

${bin}/mafExtractor: src/mafExtractor.c ${dependencies} ${API}
	mkdir -p $(dir $@)
	${cc} ${args} -O3 $< ${API} -o $@.tmp -lm
	mv $@.tmp $@

test/mafExtractor: src/mafExtractor.c ${dependencies} ${testAPI}
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

test: buildVersion test/allTests test/mafExtractor
	./test/allTests && python2.7 src/test.mafExtractor.py --verbose && rmdir ./tempTestDir && rm -rf ./test/

test/allTests: src/allTests.c ${testObjects} ${testAPI}
	mkdir -p $(dir $@)
	${cc} $^ -o $@.tmp ${args} -g -O0 -lm
	mv $@.tmp $@

test/test.mafExtractor.o: src/test.mafExtractor.c src/test.mafExtractor.h ${testAPI}
	mkdir -p $(dir $@)
	${cc} -c $< -o $@.tmp ${args} -g -O0 -lm
	mv $@.tmp $@

../external/CuTest.a: ../external/CuTest.c ../external/CuTest.h
	${cc} -c ${args} $<
	ar rc CuTest.a CuTest.o
	ranlib CuTest.a
	rm -f CuTest.o
	mv CuTest.a $@
