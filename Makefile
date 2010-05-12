include ../../../include.mk
libPath = ../../../lib
binPath = ../../../bin

all: ${binPath}/eval_PhyloComparator ${binPath}/eval_MAFComparator ${binPath}/eval_mergeMAFComparatorResults.py

${binPath}/eval_MAFComparator : eval_MAFComparator.c eval_ComparatorAPI.c disjointset.c ${libPath}/sonLib.a 
	${cxx} ${cflags} -I ${libPath} -o ${binPath}/eval_MAFComparator eval_MAFComparator.c eval_ComparatorAPI.c disjointset.c ${libPath}/sonLib.a

${binPath}/eval_PhyloComparator : eval_PhyloComparator.c eval_ComparatorAPI.c disjointset.c ${libPath}/sonLib.a 
	${cxx} ${cflags} -I ${libPath} -o ${binPath}/eval_PhyloComparator eval_PhyloComparator.c disjointset.c eval_ComparatorAPI.c ${libPath}/sonLib.a

${binPath}/eval_mergeMAFComparatorResults.py : eval_mergeMAFComparatorResults.py
	cp eval_mergeMAFComparatorResults.py ${binPath}/eval_mergeMAFComparatorResults.py
	chmod +x ${binPath}/eval_mergeMAFComparatorResults.py
	
clean :
	rm -f *.o
	rm -f ${binPath}/eval_MAFComparator
	rm -f ${binPath}/eval_PhyloComparator
