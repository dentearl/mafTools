include ../include.mk
sonLibPath = ../../sonLib/lib
binPath = ../bin

extraAPI = src/cString.c ../include/sharedMaf.o ../external/CuTest.a ../include/common.o src/comparatorRandom.o src/comparatorAPI.o ${sonLibPath}/sonLib.a 
testAPI = src/cString.c test/sharedMaf.o ../external/CuTest.a test/common.o test/comparatorRandom.o test/comparatorAPI.o ${sonLibPath}/sonLib.a 
progs =  $(foreach f,mafComparator, ${binPath}/$f)
testObjects = test/test.comparatorAPI.o test/test.comparatorRandom.o
testFlags = -O0 -g -Wall -Werror --pedantic -lm

.PHONY: all clean test

all: ${progs}

${binPath}/mafComparator: src/mafComparator.c ${extraAPI} $(wildcard src/*.h)
	@mkdir -p $(dir $@)
	${cxx} -I ${sonLibPath} -I ../include -I ../external -o $@.tmp $^ ${cflags}
	mv $@.tmp $@

test/mafComparator: src/mafComparator.c ${testAPI} $(wildcard src/*.h)
	@mkdir -p $(dir $@)
	${cxx} -I ${sonLibPath} -I ../include -I ../external -o $@.tmp $^ ${testFlags}
	mv $@.tmp $@

${binPath}/%.py: src/%.py
	@mkdir -p $(dir $@)
	cp $< $@.tmp
	chmod +x $@.tmp
	mv $@.tmp $@

%.o: %.c %.h
	${cxx} -c $< -I ${sonLibPath} -I ../external -I ../include -o $@.tmp ${cflags}
	mv $@.tmp $@
test/%.o: ../include/%.c ../include/%.h
	mkdir -p $(dir $@)
	${cxx} -c $< -o $@.tmp -I ${sonLibPath} -I ../external -I ../include -I ./src ${testFlags}
	mv $@.tmp $@
test/test.comparatorAPI.o: src/test.comparatorAPI.c src/test.comparatorAPI.h test/comparatorAPI.o
	mkdir -p $(dir $@)
	${cxx} -c $< -o $@.tmp -I ${sonLibPath} -I ../external -I ../include -I ./src ${testFlags}
	mv $@.tmp $@
test/%.o: src/%.c src/%.h ${sonLibPath}/sonLib.a
	mkdir -p $(dir $@)
	${cxx} -c $< -I ${sonLibPath} -o $@.tmp -I ../external -I ../include ${testFlags}
	mv $@.tmp $@

../external/CuTest.a: ../external/CuTest.c ../external/CuTest.h
	${cxx} -c $< ${cflags}
	ar rc CuTest.a CuTest.o
	ranlib CuTest.a
	rm -f CuTest.o
	mv CuTest.a $@

test: test/allTests test/mafComparator test/testRand
	./test/allTests && python2.7 src/test.mafComparator.py --verbose  && rm -rf ./test/ && rmdir ./tempTestDir

test/allTests: src/allTests.c ${testAPI} ${testObjects} ${sonLibPath}/sonLib.a
	mkdir -p $(dir $@)
	${cxx} $^ -I ${sonLibPath} -I ../include -I ../external -o $@.tmp ${testFlags}
	mv $@.tmp $@

# to actually use the testRand program, comment out the rm -rf on the "test:" rule and run "make test"
test/testRand: src/testRand.c ${testAPI} ${sonLibPath}/sonLib.a 
	mkdir -p $(dir $@)
	${cxx} $^ -I ${sonLibPath} -I ../include -I ../external -o $@.tmp ${testFlags}
	mv $@.tmp $@

clean:
	rm -f *.o ${progs} src/*.o && rm -rf ./test/
