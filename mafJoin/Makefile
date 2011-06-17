include ../include.mk
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
	@mkdir -p $(dir $@)
	${CC} ${cflags} -I ${incPath} -I ${sonLibPath} -I ${kentInc} -o $@.tmp $^ ${kentLibWeb} ${sonLibPath}/sonLib.a
	mv $@.tmp $@

${binPath}/mafMultiParentCheck: ${mafMultiParentCheckObjs}
	@mkdir -p $(dir $@)
	${CC} ${cflags} -I ${incPath} -I ${sonLibPath} -I ${kentInc} -o $@.tmp $^ ${kentLibWeb} ${sonLibPath}/sonLib.a
	mv $@.tmp $@

${binPath}/mafOverlap: ${mafOverlapObjs}
	@mkdir -p $(dir $@)
	${CC} ${cflags} -I ${incPath} -I ${sonLibPath} -I ${kentInc} -o $@.tmp $^ ${kentLibWeb} ${sonLibPath}/sonLib.a
	mv $@.tmp $@

${binPath}/mafAdjust: ${mafAdjustObjs}
	@mkdir -p $(dir $@)
	${CC} ${cflags} -I ${incPath} -I ${sonLibPath} -I ${kentInc} -o $@.tmp $^ ${kentLibWeb} ${sonLibPath}/sonLib.a
	mv $@.tmp $@

clean: 
	rm -f ${mafJoinObjs} ${mafOverlapObjs} ${mafAdjustObjs} ${mafMultiParentCheckObjs} ${progs}

${objs}: ${incPath}/*.h

test:
	python test/mafJoinTest.py

savebak:
	savebak mafJoin *.[ch] mafJoinTest.py Makefile 
