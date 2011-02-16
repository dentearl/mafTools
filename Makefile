include include.mk
libPath = ../sonLib/lib
binPath = bin

extraAPI = cString.c disjointset.c
progs =  ${binPath}/eval_MAFComparator ${binPath}/eval_mergeMAFComparatorResults.py ${binPath}/eval_setDiffDroppedMissing.py ${binPath}/eval_setDiffStatGenerator.py ${binPath}/eval_getRepeatBed ${binPath}/eval_MFAToMAF
# ${binPath}/eval_PhyloComparator

all: ${progs}

${binPath}/eval_MAFComparator : eval_MAFComparator.c ${extraAPI} eval_ComparatorAPI.c *.h ${libPath}/sonLib.a 
	mkdir -p $(dir $@)
	${cxx} ${cflags} -I ${libPath} -o ${binPath}/eval_MAFComparator ${extraAPI} eval_MAFComparator.c eval_ComparatorAPI.c ${libPath}/sonLib.a

${binPath}/eval_PhyloComparator : eval_PhyloComparator.c ${extraAPI} eval_ComparatorAPI.c *.h ${libPath}/sonLib.a 
	mkdir -p $(dir $@)
	${cxx} ${cflags} -I ${libPath} -o ${binPath}/eval_PhyloComparator ${extraAPI} eval_PhyloComparator.c eval_ComparatorAPI.c ${libPath}/sonLib.a

${binPath}/eval_mergeMAFComparatorResults.py : eval_mergeMAFComparatorResults.py
	mkdir -p $(dir $@)
	cp eval_mergeMAFComparatorResults.py ${binPath}/eval_mergeMAFComparatorResults.py
	chmod +x ${binPath}/eval_mergeMAFComparatorResults.py

${binPath}/eval_MFAToMAF : eval_MFAToMAF.c eTreeExtras.c ${libPath}/sonLib.a 
	mkdir -p $(dir $@)
	${cxx} ${cflags} -I ${libPath} ${tokyoCabinetIncl} -I ./ -o ${binPath}/eval_MFAToMAF eval_MFAToMAF.c eTreeExtras.c ${libPath}/sonLib.a 

${binPath}/eval_getRepeatBed : eval_getRepeatBed.py
	mkdir -p $(dir $@)
	cp eval_getRepeatBed.py ${binPath}/eval_getRepeatBed
	chmod +x ${binPath}/eval_getRepeatBed

${binPath}/eval_setDiffDroppedMissing.py : eval_setDiffDroppedMissing.py
	mkdir -p $(dir $@)
	cp $< $@
	chmod +x $@

${binPath}/eval_setDiffStatGenerator.py : eval_setDiffStatGenerator.py
	mkdir -p $(dir $@)
	cp $< $@
	chmod +x $@

clean :
	rm -f *.o ${progs}
