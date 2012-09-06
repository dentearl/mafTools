include ../include.mk
SHELL:=/bin/bash
bin = ../bin
PROGS = mafStats
dependencies = $(wildcard ../include/common.*) $(wildcard ../include/sharedMaf.*)
objects := ../include/common.o ../include/sharedMaf.o ../external/CuTest.a src/test.mafStats.o ${sonLibPath}/sonLib.a src/buildVersion.o
testObjects := test/sharedMaf.o test/common.o ../external/CuTest.a src/test.mafStats.o ${sonLibPath}/sonLib.a test/buildVersion.o

.PHONY: all clean test buildVersion

all: buildVersion $(foreach f,${PROGS}, ${bin}/$f)
buildVersion:
	@python ../include/createVersionSources.py

${bin}/mafStats: src/mafStats.c ${dependencies} ${objects}
	mkdir -p $(dir $@)
	${cxx} $< ${objects} -o $@.tmp ${cflags}
	mv $@.tmp $@

test/mafStats: src/mafStats.c ${dependencies} ${testObjects}
	mkdir -p $(dir $@)
	${cxx} $< src/allTests.c ${testObjects} -o $@.tmp ${testFlags}
	mv $@.tmp $@
%.o: %.c %.h
	${cxx} -c $< -o $@.tmp ${cflags}
	mv $@.tmp $@
test/%.o: ../include/%.c ../include/%.h
	mkdir -p $(dir $@)
	${cxx} -c $< -o $@.tmp ${testFlags}
	mv $@.tmp $@
test/%.o: src/%.c src/%.h
	mkdir -p $(dir $@)
	${cxx} -c $< -o $@.tmp ${testFlags}
	mv $@.tmp $@

clean:
	rm -rf $(foreach f,${PROGS}, ${bin}/$f) src/*.o test/ 

test: test/allTests
	$< && rm -rf ./test/

test/allTests: src/allTests.c ${testObjects}
	mkdir -p $(dir $@)
	${cxx} $^ -o $@.tmp ${testFlags}
	mv $@.tmp $@

../external/CuTest.a: ../external/CuTest.c ../external/CuTest.h
	${cxx} -c $< ${cflags}
	ar rc CuTest.a CuTest.o
	ranlib CuTest.a
	rm -f CuTest.o
	mv CuTest.a $@
