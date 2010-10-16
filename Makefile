
include ../../../include.mk
binPath = ../../../bin
libPath = ../../../lib

#cflags = ${cflags_opt}
cflags = ${cflags_dbg}

cflags += -I${libPath}
ifneq ($(wildcard ${kentLibWeb}),)
objs = jkmaf.o genome.o mafTree.o malnComp.o malnBlk.o malnBlkCursor.o malnBlkMap.o malnSet.o malnJoinBlks.o malnJoinDups.o malnJoinSets.o malnMergeComps.o malnMultiParents.o malnCompCompMap.o mafJoin.o
cflags += -I ${kentInc}
progs = ${binPath}/mafJoin
endif

# holy didn't read the make manual batman
CFLAGS=${cflags} -std=c99 -pedantic

all: ${objs} ${progs}

${binPath}/mafJoin: ${objs}
	${CC} ${cflags} -I ${libPath} -I ${kentInc} -o $@ $^ ${kentLibWeb} ${libPath}/sonLib.a


clean: 
	rm -f ${objs}

${objs}: *.h

savebak:
	savebak mafJoin *.[ch] simulation_mafJoinTest.py Makefile 
