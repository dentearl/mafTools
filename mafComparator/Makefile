include ../include.mk
sonLibPath = ../../sonLib/lib
binPath = ../bin

extraAPI = src/cString.c src/disjointset.c ../include/sharedMaf.o ../external/CuTest.a ../include/common.o src/comparatorRandom.o
testAPI = src/cString.c src/disjointset.c test/sharedMaf.o ../external/CuTest.a test/common.o test/comparatorRandom.o
progs =  $(foreach f,mafComparator mergeMafComparatorResults.py getRepeatBed mfaToMaf, ${binPath}/$f)
testObjects = src/test.comparatorAPI.o src/test.comparatorRandom.o
testArgs = -O0 -g -Wall -Werror --pedantic -lm

.PHONY: all clean test

all: ${progs}

${binPath}/mafComparator : src/mafComparator.c ${extraAPI} src/comparatorAPI.c $(wildcard src/*.h) ${sonLibPath}/sonLib.a 
	@mkdir -p $(dir $@)
	${cxx} -I ${sonLibPath} -I ../include -I ../external -o $@.tmp $^ ${cflags}
	mv $@.tmp $@

test/mafComparator : src/mafComparator.c ${testAPI} src/comparatorAPI.c $(wildcard src/*.h) ${sonLibPath}/sonLib.a
	@mkdir -p $(dir $@)
	${cxx} -I ${sonLibPath} -I ../include -I ../external -o $@.tmp $^ ${testArgs}
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
test/%.o: ../include/%.c ../include/%.h
	mkdir -p $(dir $@)
	${cxx} -g -O0 -c $< -o $@.tmp -I ../external -I ../include
	mv $@.tmp $@
test/%.o: src/%.c src/%.h
	mkdir -p $(dir $@)
	${cxx} -g -O0 -c $< -I ${sonLibPath} -o $@.tmp -I ../external -I ../include
	mv $@.tmp $@

../external/CuTest.a: ../external/CuTest.c ../external/CuTest.h
	${cxx} -c $< ${cflags}
	ar rc CuTest.a CuTest.o
	ranlib CuTest.a
	rm -f CuTest.o
	mv CuTest.a $@

test : test/allTests test/testRand test/mafComparator
#	./test/allTests && python2.7 src/test.mafComparator.py --verbose && rm -rf ./test/ && rmdir ./tempTestDir
	python2.7 src/test.mafComparator.py --verbose && rm -rf ./test/ && rmdir ./tempTestDir

test/allTests: src/allTests.c src/comparatorAPI.c ${extraAPI} ${testObjects} ${sonLibPath}/sonLib.a
	mkdir -p $(dir $@)
	${cxx} -g -O0 $^ -I ${sonLibPath} -I ../include -I ../external -o $@.tmp ${cflags}
	mv $@.tmp $@

test/testRand: src/testRand.c ${extraAPI} ${sonLibPath}/sonLib.a
	mkdir -p $(dir $@)
	${cxx} -g -O0 -pg $^ -I ${sonLibPath} -I ../include -I ../external -o $@.tmp ${cflags}
	mv $@.tmp $@

clean :
	rm -f *.o ${progs} src/*.o && rm -rf ./test/
