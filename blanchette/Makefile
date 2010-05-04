include ../../../../include.mk
binPath = ../../../../bin
libPath = ../../../../lib

all : ${binPath}/eval_MFAToMAF

${binPath}/eval_MFAToMAF : eval_MFAToMAF.c ${libPath}/sonLib.a 
	${cxx} ${cflags} -I ${libPath} ${tokyoCabinetIncl} -o ${binPath}/eval_MFAToMAF eval_MFAToMAF.c  ${libPath}/sonLib.a 

clean :
	rm -f *.o ${binPath}/eval_MFAToMAF
