
include ../../../include.mk
binPath = ../../../bin
libPath = ../../../lib

#cflags = ${cflags_opt}
cflags = ${cflags_dbg}

cflags += -I${libPath}
ifneq ($(wildcard ${kentLibWeb}),)
mafJoinObjs = mafJoin.o jkmaf.o genome.o mafTree.o malnComp.o malnBlk.o malnBlkCursor.o malnBlkSet.o malnSet.o malnJoinBlks.o malnJoinDups.o malnJoinSets.o malnMergeComps.o malnMultiParents.o malnCompCompMap.o
mafOverlapObjs = mafOverlap.o jkmaf.o
cflags += -I ${kentInc}
progs = ${binPath}/mafJoin ${binPath}/mafOverlap
endif

# holy didn't read the make manual batman
CFLAGS=${cflags} -std=c99 -pedantic

all: ${mafJoinObjs} ${progs}

${binPath}/mafJoin: ${mafJoinObjs}
	${CC} ${cflags} -I ${libPath} -I ${kentInc} -o $@ $^ ${kentLibWeb} ${libPath}/sonLib.a

${binPath}/mafOverlap: ${mafOverlapObjs}
	${CC} ${cflags} -I ${libPath} -I ${kentInc} -o $@ $^ ${kentLibWeb} ${libPath}/sonLib.a


clean: 
	rm -f ${mafJoinObjs} ${mafOverlapObjs}

${objs}: *.h

savebak:
	savebak mafJoin *.[ch] simulation_mafJoinTest.py Makefile 
