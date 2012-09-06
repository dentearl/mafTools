include ../include.mk
SHELL:=/bin/bash
bin = ../bin
PROGS = mafTransitiveClosure
dependencies = $(wildcard ../include/common.*) $(wildcard ../include/sharedMaf.*) $(wildcard ${sonLibPath}/*) ${sonLibPath}/sonLib.a ${sonLibPath}/stPinchesAndCacti.a src/allTests.c
objects := ../include/common.o ../include/sharedMaf.o ${sonLibPath}/sonLib.a ../external/CuTest.a ${sonLibPath}/stPinchesAndCacti.a src/test.mafTransitiveClosure.o src/buildVersion.o
testObjects := test/sharedMaf.o test/common.o ${sonLibPath}/sonLib.a ../external/CuTest.a ${sonLibPath}/stPinchesAndCacti.a src/test.mafTransitiveClosure.o test/buildVersion.o

.PHONY: all clean test buildVersion

all: buildVersion $(foreach f,${PROGS}, ${bin}/$f)
buildVersion:
	@python ../include/createVersionSources.py

${bin}/mafTransitiveClosure: src/mafTransitiveClosure.c ${dependencies} ${objects}
	mkdir -p $(dir $@)
	${cxx} $< src/allTests.c ${objects} -o $@.tmp ${cflags}
	mv $@.tmp $@

test/mafTransitiveClosure: src/mafTransitiveClosure.c ${dependencies} ${testObjects}
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

test: test/mafTransitiveClosure
	$< --test && python2.7 src/test.mafTransitiveClosure.py --verbose && rm -rf ./test/ && rmdir ./tempTestDir

../external/CuTest.a: ../external/CuTest.c ../external/CuTest.h
	${cxx} -c $< ${cflags}
	ar rc CuTest.a CuTest.o
	ranlib CuTest.a
	rm -f CuTest.o
	mv CuTest.a $@
