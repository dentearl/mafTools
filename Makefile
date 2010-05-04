include ../../../include.mk
libPath = ../../../lib
binPath = ../../../bin

all: eval_PhyloComparator eval_MAFComparator

eval_MAFComparator : eval_MAFComparator.c eval_ComparatorAPI.c ${libPath}/sonLib.a 
	${cxx} ${cflags} -I ${libPath} -o ${binPath}/eval_MAFComparator eval_MAFComparator.c  eval_ComparatorAPI.c disjointset.o ${libPath}/sonLib.a

eval_PhyloComparator : eval_PhyloComparator.c eval_ComparatorAPI.o disjointset.o ${libPath}/sonLib.a 
	${cxx} ${cflags} -I ${libPath} -o ${binPath}/eval_PhyloComparator eval_PhyloComparator.c disjointset.o eval_ComparatorAPI.o ${libPath}/sonLib.a

disjointset.o : disjointset.c
	${cxx} ${cflags} -I ${libPath} disjointset.c -c -o disjointset.o

eval_ComparatorAPI.o : eval_ComparatorAPI.c ${libPath}/sonLib.a
	${cxx} ${cflags} -I ${libPath} -c -o eval_ComparatorAPI.o eval_ComparatorAPI.c

clean :
	rm -f *.o
	rm -f ${binPath}/eval_MAFComparator
	rm -f ${binPath}/eval_PhyloComparator
