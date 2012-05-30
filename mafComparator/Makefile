include ../include.mk
sonLibPath = ../../sonLib/lib
binPath = ../bin

extraAPI = src/cString.c src/disjointset.c
progs =  $(foreach f,mafComparator mergeMafComparatorResults.py setDiffDroppedMissing.py setDiffStatGenerator.py getRepeatBed mfaToMaf, ${binPath}/$f)
# ${binPath}/src/PhyloComparator

.PHONY: all clean test

all: ${progs}

${binPath}/mafComparator : src/mafComparator.c ${extraAPI} src/ComparatorAPI.c $(wildcard src/*.h) ${sonLibPath}/sonLib.a 
	@mkdir -p $(dir $@)
	${cxx} ${cflags} -I ${sonLibPath} -o $@.tmp $^ -lm
	mv $@.tmp $@

${binPath}/PhyloComparator : src/PhyloComparator.c ${extraAPI} src/ComparatorAPI.c $(wildcard src/*.h) ${sonLibPath}/sonLib.a 
	@mkdir -p $(dir $@)
	${cxx} ${cflags} -I ${sonLibPath} -o $@.tmp $^
	mv $@.tmp $@

${binPath}/mfaToMaf : src/mfaToMaf.c src/eTreeExtras.c ${sonLibPath}/sonLib.a 
	@mkdir -p $(dir $@)
	${cxx} ${cflags} -I ${sonLibPath} ${tokyoCabinetIncl} -I ./ -o $@.tmp $^ -lm
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

test : all
	python src/allTests.py -v

clean :
	rm -f *.o ${progs}
