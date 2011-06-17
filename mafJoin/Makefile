include ../../../include.mk
sonLibPath = ../../sonLib/lib
binPath = ../bin
incPath = inc/

cflags = ${cflags_opt} -DNDEBUG
#cflags = ${cflags_dbg}
#cflags = ${cflags_dbg} -DNDEBUG
#cflags = ${cflags_dbg} -DASSERT_SLOW

cflags += -I${sonLibPath}
ifneq ($(wildcard ${kentLibWeb}),)
mafJoinObjs = mafJoin.o jkmaf.o genome.o mafTree.o malnComp.o malnBlk.o malnBlkCursor.o malnBlkSet.o malnSet.o malnJoinBlks.o malnJoinWithinSet.o malnJoinSets.o malnMergeComps.o malnMultiParents.o malnCompCompMap.o mafIO.o
mafOverlapObjs = mafOverlap.o jkmaf.o
mafAdjustObjs = mafAdjust.o jkmaf.o genome.o mafTree.o malnComp.o malnBlk.o malnBlkCursor.o malnBlkSet.o malnSet.o malnJoinBlks.o malnJoinWithinSet.o malnJoinSets.o malnMergeComps.o malnCompCompMap.o mafIO.o
mafMultiParentCheckObjs = mafMultiParentCheck.o jkmaf.o genome.o mafTree.o malnComp.o malnBlk.o malnBlkCursor.o malnBlkSet.o malnSet.o malnCompCompMap.o malnMultiParents.o mafIO.o
cflags += -I ${kentInc}
progs = ${binPath}/mafJoin ${binPath}/mafOverlap ${binPath}/mafMultiParentCheck ${binPath}/mafAdjust
endif

# holy didn't read the make manual batman
CFLAGS=${cflags} -std=c99 -pedantic


all: ${mafJoinObjs} ${progs}

${mafJoinObjs} ${mafOverlapObjs}: $(wildcard *.h)


${binPath}/mafJoin: ${mafJoinObjs}
	${CC} ${cflags} -I ${incPath} -I ${sonLibPath} -I ${kentInc} -o $@ $^ ${kentLibWeb} ${sonLibPath}/sonLib.a

${binPath}/mafMultiParentCheck: ${mafMultiParentCheckObjs}
	${CC} ${cflags} -I ${incPath} -I ${sonLibPath} -I ${kentInc} -o $@ $^ ${kentLibWeb} ${sonLibPath}/sonLib.a

${binPath}/mafOverlap: ${mafOverlapObjs}
	${CC} ${cflags} -I ${incPath} -I ${sonLibPath} -I ${kentInc} -o $@ $^ ${kentLibWeb} ${sonLibPath}/sonLib.a

${binPath}/mafAdjust: ${mafAdjustObjs}
	${CC} ${cflags} -I ${incPath} -I ${sonLibPath} -I ${kentInc} -o $@ $^ ${kentLibWeb} ${sonLibPath}/sonLib.a

clean: 
	rm -f ${mafJoinObjs} ${mafOverlapObjs} ${mafAdjustObjs} ${mafMultiParentCheckObjs} ${progs}

${objs}: ${incPath}/*.h

test:
	python test/mafJoinTest.py

savebak:
	savebak mafJoin *.[ch] mafJoinTest.py Makefile 
