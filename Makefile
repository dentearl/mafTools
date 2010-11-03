include ../../../include.mk
libPath = ../../../lib
binPath = ../../../bin

extraAPI = cString.c disjointset.c

all: ${binPath}/eval_PhyloComparator ${binPath}/eval_MAFComparator ${binPath}/eval_mergeMAFComparatorResults.py

${binPath}/eval_MAFComparator : eval_MAFComparator.c ${extraAPI} eval_ComparatorAPI.c *.h ${libPath}/sonLib.a 
	${cxx} ${cflags} -I ${libPath} -o ${binPath}/eval_MAFComparator ${extraAPI} eval_MAFComparator.c eval_ComparatorAPI.c ${libPath}/sonLib.a

${binPath}/eval_PhyloComparator : eval_PhyloComparator.c ${extraAPI} eval_ComparatorAPI.c *.h ${libPath}/sonLib.a 
	${cxx} ${cflags} -I ${libPath} -o ${binPath}/eval_PhyloComparator ${extraAPI} eval_PhyloComparator.c eval_ComparatorAPI.c ${libPath}/sonLib.a

${binPath}/eval_mergeMAFComparatorResults.py : eval_mergeMAFComparatorResults.py
	cp eval_mergeMAFComparatorResults.py ${binPath}/eval_mergeMAFComparatorResults.py
	chmod +x ${binPath}/eval_mergeMAFComparatorResults.py

${binPath}/eval_intersectDroppedMissing.py : eval_intersectDroppedMissing.py
	cp eval_intersectDroppedMissing.py ${binPath}/eval_intersectDroppedMissing.py
	chmod +x ${binPath}/eval_intersectDroppedMissing.py

clean :
	rm -f *.o
	rm -f ${binPath}/eval_MAFComparator
	rm -f ${binPath}/eval_PhyloComparator
