include ../include.mk
sonLibPath = ../../sonLib/lib
binPath = ../bin

extraAPI = src/cString.c src/disjointset.c ../include/sharedMaf.o ../external/CuTest.a ../include/common.o
progs =  $(foreach f,mafComparator mergeMafComparatorResults.py setDiffDroppedMissing.py setDiffStatGenerator.py getRepeatBed mfaToMaf, ${binPath}/$f)
# ${binPath}/src/PhyloComparator

.PHONY: all clean test

all: ${progs}

${binPath}/mafComparator : src/mafComparator.c ${extraAPI} src/ComparatorAPI.c $(wildcard src/*.h) ${sonLibPath}/sonLib.a 
	@mkdir -p $(dir $@)
	${cxx} -I ${sonLibPath} -I ../include -I ../external -o $@.tmp $^ ${cflags}
	mv $@.tmp $@

${binPath}/PhyloComparator : src/PhyloComparator.c ${extraAPI} src/ComparatorAPI.c $(wildcard src/*.h) ${sonLibPath}/sonLib.a 
	@mkdir -p $(dir $@)
	${cxx} -I ${sonLibPath} -o $@.tmp $^ ${cflags}
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
	${cxx} -O3 -c -O3 $< -I ../external -o $@.tmp ${cflags}
	mv $@.tmp $@

../external/CuTest.a: ../external/CuTest.c ../external/CuTest.h
	${cxx} -c $< ${cflags}
	ar rc CuTest.a CuTest.o
	ranlib CuTest.a
	rm -f CuTest.o
	mv CuTest.a $@

test : all
	python src/allTests.py -v

clean :
	rm -f *.o ${progs}
