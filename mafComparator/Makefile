include ../include.mk
sonLibPath = ../../sonLib/lib
binPath = ../bin

extraAPI = src/cString.c src/disjointset.c ../include/sharedMaf.o ../external/CuTest.a ../include/common.o
progs =  $(foreach f,mafComparator mergeMafComparatorResults.py setDiffDroppedMissing.py setDiffStatGenerator.py getRepeatBed mfaToMaf, ${binPath}/$f)
testObjects = src/test.ComparatorAPI.o 

.PHONY: all clean test

all: ${progs}

${binPath}/mafComparator : src/mafComparator.c ${extraAPI} src/ComparatorAPI.c $(wildcard src/*.h) ${sonLibPath}/sonLib.a 
	@mkdir -p $(dir $@)
	${cxx} -I ${sonLibPath} -I ../include -I ../external -o $@.tmp $^ ${cflags}
	mv $@.tmp $@

${binPath}/mfaToMaf : src/mfaToMaf.c src/eTreeExtras.c ${sonLibPath}/sonLib.a 
	@mkdir -p $(dir $@)
	${cxx} -I ${sonLibPath} ${tokyoCabinetIncl} -I ./ -o $@.tmp $^ ${cflags}
	mv $@.tmp $@

${binPath}/getRepeatBed: src/getRepeatBed.py
	@mkdir -p $(dir $@)
	cp $< $@
	chmod +x $@

${binPath}/%.py : src/%.py
	@mkdir -p $(dir $@)
	cp $< $@.tmp
	chmod +x $@.tmp
	mv $@.tmp $@

%.o: %.c %.h
	${cxx} -O3 -c -O3 $< -I ${sonLibPath} -I ../external -I ../include -o $@.tmp ${cflags}
	mv $@.tmp $@

../external/CuTest.a: ../external/CuTest.c ../external/CuTest.h
	${cxx} -c $< ${cflags}
	ar rc CuTest.a CuTest.o
	ranlib CuTest.a
	rm -f CuTest.o
	mv CuTest.a $@

test : all test/allTests test/testRand
	./test/allTests && python2.7 src/allTests.py --verbose


test/allTests: src/allTests.c src/ComparatorAPI.c ${extraAPI} ${testObjects} ${sonLibPath}/sonLib.a
	mkdir -p $(dir $@)
	${cxx} -g -O0 $^ -I ${sonLibPath} -I ../include -I ../external -o $@.tmp ${cflags}
	mv $@.tmp $@

test/testRand: src/testRand.c src/ComparatorRandom.c ${extraAPI} ${sonLibPath}/sonLib.a
	mkdir -p $(dir $@)
	${cxx} -g -O0 $^ -I ${sonLibPath} -I ../include -I ../external -o $@.tmp ${cflags}
	mv $@.tmp $@

clean :
	rm -f *.o ${progs}
