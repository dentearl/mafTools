
include ../../../include.mk
binPath = ../../../bin
libPath = ../../../lib

cflags += -I../../sonLib
ifneq ($(wildcard ${kentLibWeb}),)
objs = jkmaf.o genome.o mafInvert.o mafInvertJoin.o mafInvertWrite.o mafJoin.o
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

