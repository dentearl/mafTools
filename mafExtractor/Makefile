include ../inc/common.mk
SHELL:=/bin/bash
bin = ../bin
inc = ../inc
lib = ../lib
PROGS = mafExtractor
dependencies = ${inc}/common.h ${inc}/sharedMaf.h ${lib}/common.c ${lib}/sharedMaf.c  src/mafExtractor.h
API = ${lib}/common.o ${lib}/sharedMaf.o ../external/CuTest.a src/mafExtractorAPI.o src/buildVersion.o
testAPI = test/sharedMaf.o ../external/CuTest.a test/common.o test/mafExtractorAPI.o test/buildVersion.o
testObjects := test/test.mafExtractor.o
sources = src/mafExtractor.c src/mafExtractor.h

.PHONY: all clean test buildVersion

all: buildVersion $(foreach f,${PROGS}, ${bin}/$f)
buildVersion: src/buildVersion.c
src/buildVersion.c: ${sources} ${dependencies}
	@python ../lib/createVersionSources.py

../lib/%.o: ../lib/%.c ../inc/%.h
	cd ../lib/ && make

${bin}/mafExtractor: src/mafExtractor.c ${dependencies} ${API}
	mkdir -p $(dir $@)
	${cxx} ${cflags} -O3 $< ${API} -o $@.tmp -lm
	mv $@.tmp $@

test/mafExtractor: src/mafExtractor.c ${dependencies} ${testAPI}
	mkdir -p $(dir $@)
	${cxx} ${cflags} -g -O0 $< ${testAPI} -o $@.tmp -lm
	mv $@.tmp $@

%.o: %.c %.h
	${cxx} -O3 -c ${cflags} $< -o $@.tmp
	mv $@.tmp $@

test/%.o: ${lib}/%.c ${inc}/%.h
	mkdir -p $(dir $@)
	${cxx} -c $< -o $@.tmp ${cflags} -g -O0
	mv $@.tmp $@
test/%.o: src/%.c src/%.h
	mkdir -p $(dir $@)
	${cxx} -c $< -o $@.tmp ${cflags} -g -O0
	mv $@.tmp $@

clean:
	rm -rf $(foreach f,${PROGS}, ${bin}/$f) src/*.o test/ temTestDir/ src/buildVersion.c src/buildVersion.h

test: buildVersion test/allTests test/mafExtractor
	./test/allTests && python2.7 src/test.mafExtractor.py --verbose && rmdir ./tempTestDir && rm -rf ./test/

test/allTests: src/allTests.c ${testObjects} ${testAPI}
	mkdir -p $(dir $@)
	${cxx} $^ -o $@.tmp ${cflags} -g -O0
	mv $@.tmp $@

test/test.mafExtractor.o: src/test.mafExtractor.c src/test.mafExtractor.h ${testAPI}
	mkdir -p $(dir $@)
	${cxx} -c $< -o $@.tmp ${cflags} -I src/ -g -O0
	mv $@.tmp $@

../external/CuTest.a: ../external/CuTest.c ../external/CuTest.h
	${cxx} -c ${cflags} $<
	ar rc CuTest.a CuTest.o
	ranlib CuTest.a
	rm -f CuTest.o
	mv CuTest.a $@
