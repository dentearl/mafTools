
include ../../../include.mk
binPath = ../../../bin
libPath = ../../../lib

cflags = ${cflags_opt}
cflags += -I../../sonLib/inc
ifneq ($(wildcard ${kentLibWeb}),)
objs = jkmaf.o genome.o mafTree.o malnComp.o malnBlk.o malnBlkCursor.o malnSet.o malnJoinBlks.o malnJoinDups.o malnJoinSets.o mafJoin.o
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
