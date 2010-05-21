
include ../../../include.mk
binPath = ../../../bin
libPath = ../../../lib

ifneq ($(wildcard ${kentLibWeb}),)
objs = genome.o mafInvert.o mafInvertJoin.o mafRevert.o mafJoin.o
cflags += -I ${kentInc}
progs = ${binPath}/mafJoin
endif

# holy didn't read the make manual batman
CFLAGS=${cflags} -std=c99

all: ${objs} ${progs}

${binPath}/mafJoin: ${objs}
	${CC} ${cflags} -I ${libPath} -I ${kentInc} -o $@ $^ ${kentLibWeb}


clean: 
	rm -f ${objs}

