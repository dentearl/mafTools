include ../inc/common.mk
SHELL:=/bin/bash
bin = ../bin
inc = ../inc
lib = ../lib
PROGS = mafTransitiveClosure
dependencies = ${inc}/common.h ${inc}/sharedMaf.h ${lib}/common.c ${lib}/sharedMaf.c $(wildcard ${sonLibPath}/*) ${sonLibPath}/stPinchesAndCacti.a ${sonLibPath}/sonLib.a src/allTests.c
objects := ${lib}/common.o ${lib}/sharedMaf.o ${sonLibPath}/stPinchesAndCacti.a  ${sonLibPath}/sonLib.a ../external/CuTest.a src/test.mafTransitiveClosure.o src/buildVersion.o
testObjects := test/sharedMaf.o test/common.o ${sonLibPath}/stPinchesAndCacti.a  ${sonLibPath}/sonLib.a ../external/CuTest.a src/test.mafTransitiveClosure.o test/buildVersion.o
sources := src/mafTransitiveClosure.c src/mafTransitiveClosure.h

.PHONY: all clean test buildVersion

all: buildVersion $(foreach f,${PROGS}, ${bin}/$f)
buildVersion: src/buildVersion.c
src/buildVersion.c: ${sources} ${dependencies}
	@python ${lib}/createVersionSources.py

${bin}/mafTransitiveClosure: src/mafTransitiveClosure.c ${dependencies} ${objects}
	mkdir -p $(dir $@)
	${cxx} $< src/allTests.c ${objects} -o $@.tmp ${cflags} -lm
	mv $@.tmp $@

test/mafTransitiveClosure: src/mafTransitiveClosure.c ${dependencies} ${testObjects}
	mkdir -p $(dir $@)
	${cxx} $< src/allTests.c ${testObjects} -o $@.tmp ${testFlags} -lm
	mv $@.tmp $@
%.o: %.c ${inc}/%.h
	${cxx} -c $< -o $@.tmp ${cflags}
	mv $@.tmp $@
%.o: %.c %.h
	${cxx} -c $< -o $@.tmp ${cflags}
	mv $@.tmp $@
test/%.o: ${lib}/%.c ${inc}/%.h
	mkdir -p $(dir $@)
	${cxx} -c $< -o $@.tmp ${testFlags}
	mv $@.tmp $@
test/%.o: src/%.c src/%.h
	mkdir -p $(dir $@)
	${cxx} -c $< -o $@.tmp ${testFlags}
	mv $@.tmp $@
test/allTests: src/allTests.c ${testObjects} ${sonLibPath}/sonLib.a
	mkdir -p $(dir $@)
	${cxx} $^ -o $@.tmp ${testFlags} ${lm}
	mv $@.tmp $@

clean:
	rm -rf $(foreach f,${PROGS}, ${bin}/$f) src/*.o test/ src/buildVersion.c src/buildVersion.h

test: buildVersion test/mafTransitiveClosure
	test/mafTransitiveClosure --test && python2.7 src/test.mafTransitiveClosure.py --verbose && rm -rf ./test/ && rmdir ./tempTestDir

../external/CuTest.a: ../external/CuTest.c ../external/CuTest.h
	${cxx} -c $< ${cflags}
	ar rc CuTest.a CuTest.o
	ranlib CuTest.a
	rm -f CuTest.o
	mv CuTest.a $@
