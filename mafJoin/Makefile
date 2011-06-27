include ../include.mk
sonLibPath = ../../sonLib/lib
binPath = ../bin
srcPath = ./src
objPath = ./obj

cflags = ${cflags_opt} -DNDEBUG
#cflags = ${cflags_dbg}
#cflags = ${cflags_dbg} -DNDEBUG
#cflags = ${cflags_dbg} -DASSERT_SLOW
cflags += -I ${sonLibPath} -I ${srcPath}

# required kent library stuff (git clone git://genome-source.cse.ucsc.edu/kent.git)
# if you have installed the kent library, you can change the kentDir variable
# on the command line by using $ make kentDir=/path/to/kent/src
kentDir = /hive/groups/recon/local/kent/src
ifneq ($(wildcard ${kentDir}),)
kentInc = ${kentDir}/inc
kentLib = ${kentDir}/lib
kentLibWeb = ${kentLib}/${MACH}/jkweb.a
mafJoinObjs = $(foreach f,mafJoin jkmaf genome mafTree malnComp malnBlk malnBlkCursor malnBlkSet malnSet malnJoinBlks malnJoinWithinSet malnJoinSets malnMergeComps malnMultiParents malnCompCompMap mafIO, ${objPath}/$f.o)
mafOverlapObjs = $(foreach f,mafOverlap jkmaf, ${objPath}/$f.o)
mafAdjustObjs = $(foreach f,mafAdjust jkmaf genome mafTree malnComp malnBlk malnBlkCursor malnBlkSet malnSet malnJoinBlks malnJoinWithinSet malnJoinSets malnMergeComps malnCompCompMap mafIO, ${objPath}/$f.o)
mafMultiParentCheckObjs = $(foreach f,mafMultiParentCheck jkmaf genome mafTree malnComp malnBlk malnBlkCursor malnBlkSet malnSet malnCompCompMap malnMultiParents mafIO, ${objPath}/$f.o)
cflags += -I ${kentInc}
CFLAGS = ${cflags} -std=c99 -pedantic
progs = ${binPath}/mafJoin ${binPath}/mafOverlap ${binPath}/mafMultiParentCheck ${binPath}/mafAdjust
endif

.PHONY: all clean test

all: ${progs} ${mafJoinObjs} 

${objPath}/%.o: ${srcPath}/%.c $(wildcard ${srcPath}/*.h)
	@mkdir -p $(dir $@)            
	${CC} ${CFLAGS} -c -o $@.tmp $<
	mv $@.tmp $@                   

${binPath}/mafJoin: ${mafJoinObjs}
	@mkdir -p $(dir $@)
	${CC} ${CFLAGS} -o $@.tmp $^ ${kentLibWeb} ${sonLibPath}/sonLib.a
	mv $@.tmp $@

${binPath}/mafMultiParentCheck: ${mafMultiParentCheckObjs}
	@mkdir -p $(dir $@)
	${CC} ${CFLAGS} -o $@.tmp $^ ${kentLibWeb} ${sonLibPath}/sonLib.a
	mv $@.tmp $@

${binPath}/mafOverlap: ${mafOverlapObjs}
	@mkdir -p $(dir $@)
	${CC} ${CFLAGS} -o $@.tmp $^ ${kentLibWeb} ${sonLibPath}/sonLib.a
	mv $@.tmp $@

${binPath}/mafAdjust: ${mafAdjustObjs}
	@mkdir -p $(dir $@)
	${CC} ${CFLAGS} -o $@.tmp $^ ${kentLibWeb} ${sonLibPath}/sonLib.a
	mv $@.tmp $@

clean: 
	rm -f ${progs}
	rm -rf ${objPath}

test:
	python test/mafJoinTest.py
