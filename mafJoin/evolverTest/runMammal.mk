 #
# override simOutDir to save experiments
#

# simDir copied from /hive/users/dearl/simTreeWorking/MARK_BENCODE_2M/
simRoot = evolver/mammal
simDir = ${simRoot}/simulation
simOutDir = ${simRoot}/output/expr1

include defs.mk

mafJoinOpts = -maxBlkWidth=10000

createdMafs = \
	${simOutDir}/sMouse-sRat/sMouse-sRat.maf \
	${simOutDir}/sMouse-sRat/sH-sM-sR.maf \
	${simOutDir}/sCow-sDog/sCow-sDog.maf \
	${simOutDir}/sCow-sDog/root.maf \
	${simOutDir}/sH-sM-sR/sH-sM-sR.maf \
	${simOutDir}/sH-sM-sR/root.maf \
	${simOutDir}/root/root.maf



join: ${createdMafs}
${simOutDir}/sMouse-sRat/sMouse-sRat.maf: ${simDir}/simRat/sMouse-sRat.maf ${simDir}/simMouse/sMouse-sRat.maf
	./runMafJoin ${mafJoinOpts} -treelessRoot1=sMouse-sRat -treelessRoot2=sMouse-sRat 'sMouse-sRat' ${simDir}/simRat/sMouse-sRat.maf ${simDir}/simMouse/sMouse-sRat.maf ${simOutDir}/sMouse-sRat/sMouse-sRat.maf

${simOutDir}/sMouse-sRat/sH-sM-sR.maf: ${simOutDir}/sMouse-sRat/sMouse-sRat.maf ${simDir}/sMouse-sRat/sH-sM-sR.tmp.maf 
	./runMafJoin ${mafJoinOpts} -treelessRoot2=sH-sM-sR 'sMouse-sRat' ${simOutDir}/sMouse-sRat/sMouse-sRat.maf ${simDir}/sMouse-sRat/sH-sM-sR.tmp.maf ${simOutDir}/sMouse-sRat/sH-sM-sR.maf

${simOutDir}/sCow-sDog/sCow-sDog.maf: ${simDir}/simDog/sCow-sDog.maf ${simDir}/simCow/sCow-sDog.maf
	./runMafJoin ${mafJoinOpts} -treelessRoot1=sCow-sDog -treelessRoot2=sCow-sDog 'sCow-sDog' ${simDir}/simDog/sCow-sDog.maf ${simDir}/simCow/sCow-sDog.maf ${simOutDir}/sCow-sDog/sCow-sDog.maf

${simOutDir}/sCow-sDog/root.maf: ${simOutDir}/sCow-sDog/sCow-sDog.maf ${simDir}/sCow-sDog/root.tmp.maf
	./runMafJoin ${mafJoinOpts} -treelessRoot2=root 'sCow-sDog' ${simOutDir}/sCow-sDog/sCow-sDog.maf ${simDir}/sCow-sDog/root.tmp.maf ${simOutDir}/sCow-sDog/root.maf

${simOutDir}/sH-sM-sR/sH-sM-sR.maf: ${simOutDir}/sMouse-sRat/sH-sM-sR.maf ${simDir}/simHuman/sH-sM-sR.maf
	./runMafJoin ${mafJoinOpts} -treelessRoot2=sH-sM-sR 'sH-sM-sR' ${simOutDir}/sMouse-sRat/sH-sM-sR.maf ${simDir}/simHuman/sH-sM-sR.maf ${simOutDir}/sH-sM-sR/sH-sM-sR.maf

${simOutDir}/sH-sM-sR/root.maf: ${simOutDir}/sH-sM-sR/sH-sM-sR.maf ${simDir}/sH-sM-sR/root.tmp.maf
	./runMafJoin ${mafJoinOpts} -treelessRoot2=root 'sH-sM-sR' ${simOutDir}/sH-sM-sR/sH-sM-sR.maf ${simDir}/sH-sM-sR/root.tmp.maf ${simOutDir}/sH-sM-sR/root.maf

${simOutDir}/root/root.maf: ${simOutDir}/sH-sM-sR/root.maf ${simOutDir}/sCow-sDog/root.maf
	./runMafJoin ${mafJoinOpts} 'root' ${simOutDir}/sH-sM-sR/root.maf ${simOutDir}/sCow-sDog/root.maf ${simOutDir}/root/root.maf


compareBaseNames = \
    simRat/sMouse-sRat \
    simMouse/sMouse-sRat \
    sMouse-sRat/sH-sM-sR.tmp \
    simDog/sCow-sDog \
    simCow/sCow-sDog \
    sCow-sDog/root.tmp \
    simHuman/sH-sM-sR \
    sH-sM-sR/root.tmp

include rules.mk

