include include.mk
libPath = ../sonLib/lib
binPath = bin

extraAPI = src/cString.c src/disjointset.c
progs =  ${binPath}/MAFComparator ${binPath}/mergeMAFComparatorResults.py ${binPath}/setDiffDroppedMissing.py ${binPath}/setDiffStatGenerator.py ${binPath}/getRepeatBed ${binPath}/MFAToMAF
# ${binPath}/src/PhyloComparator

.PHONY: all clean test

all: ${progs}

${binPath}/MAFComparator : src/MAFComparator.c ${extraAPI} src/ComparatorAPI.c $(wildcard src/*.h) ${libPath}/sonLib.a 
	mkdir -p $(dir $@)
	${cxx} ${cflags} -I ${libPath} -o $@.tmp ${extraAPI} src/MAFComparator.c src/ComparatorAPI.c ${libPath}/sonLib.a
	mv $@.tmp $@

${binPath}/PhyloComparator : src/PhyloComparator.c ${extraAPI} src/ComparatorAPI.c $(wildcard src/*.h) ${libPath}/sonLib.a 
	mkdir -p $(dir $@)
	${cxx} ${cflags} -I ${libPath} -o $@.tmp ${extraAPI} src/PhyloComparator.c src/ComparatorAPI.c ${libPath}/sonLib.a
	mv $@.tmp $@

${binPath}/mergeMAFComparatorResults.py : src/mergeMAFComparatorResults.py
	mkdir -p $(dir $@)
	cp $< $@
	chmod +x $@

${binPath}/MFAToMAF : src/MFAToMAF.c src/eTreeExtras.c ${libPath}/sonLib.a 
	mkdir -p $(dir $@)
	${cxx} ${cflags} -I ${libPath} ${tokyoCabinetIncl} -I ./ -o $@.tmp src/MFAToMAF.c src/eTreeExtras.c ${libPath}/sonLib.a 
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
