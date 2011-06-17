include ../include.mk
sonLibPath = ../../sonLib/lib
binPath = ../bin

extraAPI = src/cString.c src/disjointset.c
progs =  ${binPath}/mafComparator ${binPath}/mergeMafComparatorResults.py ${binPath}/setDiffDroppedMissing.py ${binPath}/setDiffStatGenerator.py ${binPath}/getRepeatBed ${binPath}/mfaToMaf
# ${binPath}/src/PhyloComparator

.PHONY: all clean test

all: ${progs}

${binPath}/mafComparator : src/mafComparator.c ${extraAPI} src/ComparatorAPI.c $(wildcard src/*.h) ${sonLibPath}/sonLib.a 
	mkdir -p $(dir $@)
	${cxx} ${cflags} -I ${sonLibPath} -o $@.tmp $^
	mv $@.tmp $@

${binPath}/PhyloComparator : src/PhyloComparator.c ${extraAPI} src/ComparatorAPI.c $(wildcard src/*.h) ${sonLibPath}/sonLib.a 
	mkdir -p $(dir $@)
	${cxx} ${cflags} -I ${sonLibPath} -o $@.tmp $^
	mv $@.tmp $@

${binPath}/mergeMafComparatorResults.py : src/mergeMafComparatorResults.py
	mkdir -p $(dir $@)
	cp $< $@
	chmod +x $@

${binPath}/mfaToMaf : src/mfaToMaf.c src/eTreeExtras.c ${sonLibPath}/sonLib.a 
	mkdir -p $(dir $@)
	${cxx} ${cflags} -I ${sonLibPath} ${tokyoCabinetIncl} -I ./ -o $@.tmp $^
	mv $@.tmp $@

${binPath}/getRepeatBed : src/getRepeatBed.py
	mkdir -p $(dir $@)
	cp $< $@
	chmod +x $@

${binPath}/setDiffDroppedMissing.py : src/setDiffDroppedMissing.py
	mkdir -p $(dir $@)
	cp $< $@
	chmod +x $@

${binPath}/setDiffStatGenerator.py : src/setDiffStatGenerator.py
	mkdir -p $(dir $@)
	cp $< $@
	chmod +x $@

test :
	python src/allTests.py -v

clean :
	rm -f *.o ${progs}
