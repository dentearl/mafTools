include ../../../include.mk
libPath = ../../../lib
binPath = ../../../bin

extraAPI = cString.c disjointset.c

all: ${binPath}/eval_PhyloComparator ${binPath}/eval_MAFComparator ${binPath}/eval_mergeMAFComparatorResults.py ${binPath}/eval_setDiffDroppedMissing.py ${binPath}/eval_setDiffStatGenerator.py ${binPath}/eval_getRepeatBed ${binPath}/eval_MFAToMAF

${binPath}/eval_MAFComparator : eval_MAFComparator.c ${extraAPI} eval_ComparatorAPI.c *.h ${libPath}/sonLib.a 
	${cxx} ${cflags} -I ${libPath} -o ${binPath}/eval_MAFComparator ${extraAPI} eval_MAFComparator.c eval_ComparatorAPI.c ${libPath}/sonLib.a

${binPath}/eval_PhyloComparator : eval_PhyloComparator.c ${extraAPI} eval_ComparatorAPI.c *.h ${libPath}/sonLib.a 
	${cxx} ${cflags} -I ${libPath} -o ${binPath}/eval_PhyloComparator ${extraAPI} eval_PhyloComparator.c eval_ComparatorAPI.c ${libPath}/sonLib.a

${binPath}/eval_mergeMAFComparatorResults.py : eval_mergeMAFComparatorResults.py
	cp eval_mergeMAFComparatorResults.py ${binPath}/eval_mergeMAFComparatorResults.py
	chmod +x ${binPath}/eval_mergeMAFComparatorResults.py

${binPath}/eval_MFAToMAF : eval_MFAToMAF.c ../mAFComparison/eTreeExtras.c ${libPath}/sonLib.a 
	${cxx} ${cflags} -I ${libPath} ${tokyoCabinetIncl} -I ../mAFComparison -o ${binPath}/eval_MFAToMAF eval_MFAToMAF.c ../mAFComparison/eTreeExtras.c ${libPath}/sonLib.a 

${binPath}/eval_getRepeatBed : eval_getRepeatBed.py
	cp eval_getRepeatBed.py ${binPath}/eval_getRepeatBed
	chmod +x ${binPath}/eval_getRepeatBed

${binPath}/eval_setDiffDroppedMissing.py : eval_setDiffDroppedMissing.py
	cp $< $@
	chmod +x $@

${binPath}/eval_setDiffStatGenerator.py : eval_setDiffStatGenerator.py
	cp $< $@
	chmod +x $@

clean :
	rm -f *.o
	rm -f ${binPath}/eval_MAFComparator
	rm -f ${binPath}/eval_PhyloComparator
	rm -f ${binPath}/eval_MFAToMAF
	rm -f ${binPath}/eval_getRepeatBed
