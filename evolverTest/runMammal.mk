#
# override simOutDir to save experiments
#

# simDir copied from /hive/users/dearl/simTreeWorking/MARK_BENCODE_2M/
simRoot = evolver/mammal
simDir = ${simRoot}/simulation
simOutDir = ${simRoot}/output/expr1

export PATH := ../../../bin:/hive/groups/recon/local/bin:${PATH}
export PYTHONPATH = ../../../src

mafJoinOpts = -maxBlkWidth=10000

createdMafs = \
	${simOutDir}/sMouse-sRat/sMouse-sRat.maf \
	${simOutDir}/sMouse-sRat/sH-sM-sR.maf \
	${simOutDir}/sCow-sDog/sCow-sDog.maf \
	${simOutDir}/sCow-sDog/root.maf \
	${simOutDir}/sH-sM-sR/sH-sM-sR.maf \
	${simOutDir}/sH-sM-sR/root.maf \
	${simOutDir}/root/root.maf


mafJoin = time -p mafJoin

join: ${createdMafs}
${simOutDir}/sMouse-sRat/sMouse-sRat.maf: ${simDir}/simRat/sMouse-sRat.maf ${simDir}/simMouse/sMouse-sRat.maf
	@mkdir -p $(dir $@)
	${mafJoin} ${mafJoinOpts} -treelessRoot1=sMouse-sRat -treelessRoot2=sMouse-sRat 'sMouse-sRat' \
	    -inMaf1Dump=${simDir}/simRat/sMouse-sRat.mafdmp -inMaf2Dump=${simDir}/simMouse/sMouse-sRat.mafdmp -outMafDump=${simOutDir}/sMouse-sRat/sMouse-sRat.mafdmp \
	    ${simDir}/simRat/sMouse-sRat.maf ${simDir}/simMouse/sMouse-sRat.maf ${simOutDir}/sMouse-sRat/sMouse-sRat.maf >& $(basename $@).err

${simOutDir}/sMouse-sRat/sH-sM-sR.maf: ${simOutDir}/sMouse-sRat/sMouse-sRat.maf ${simDir}/sMouse-sRat/sH-sM-sR.tmp.maf 
	@mkdir -p $(dir $@)
	${mafJoin} ${mafJoinOpts} -treelessRoot2=sH-sM-sR 'sMouse-sRat' \
	    -inMaf1Dump=${simOutDir}/sMouse-sRat/sMouse-sRat.mafdmp -inMaf2Dump=${simDir}/sMouse-sRat/sH-sM-sR.tmp.mafdmp -outMafDump=${simOutDir}/sMouse-sRat/sH-sM-sR.mafdmp \
	    ${simOutDir}/sMouse-sRat/sMouse-sRat.maf ${simDir}/sMouse-sRat/sH-sM-sR.tmp.maf ${simOutDir}/sMouse-sRat/sH-sM-sR.maf >& $(basename $@).err

${simOutDir}/sCow-sDog/sCow-sDog.maf: ${simDir}/simDog/sCow-sDog.maf ${simDir}/simCow/sCow-sDog.maf
	@mkdir -p $(dir $@)
	${mafJoin} ${mafJoinOpts} -treelessRoot1=sCow-sDog -treelessRoot2=sCow-sDog 'sCow-sDog' \
	    -inMaf1Dump=${simDir}/simDog/sCow-sDog.mafdmp -inMaf2Dump=${simDir}/simCow/sCow-sDog.mafdmp -outMafDump=${simOutDir}/sCow-sDog/sCow-sDog.mafdmp \
	    ${simDir}/simDog/sCow-sDog.maf ${simDir}/simCow/sCow-sDog.maf ${simOutDir}/sCow-sDog/sCow-sDog.maf >& $(basename $@).err

${simOutDir}/sCow-sDog/root.maf: ${simOutDir}/sCow-sDog/sCow-sDog.maf ${simDir}/sCow-sDog/root.tmp.maf
	@mkdir -p $(dir $@)
	${mafJoin} ${mafJoinOpts} -treelessRoot2=root 'sCow-sDog' \
	    -inMaf1Dump=${simOutDir}/sCow-sDog/sCow-sDog.mafdmp -inMaf2Dump=${simDir}/sCow-sDog/root.tmp.mafdmp -outMafDump=${simOutDir}/sCow-sDog/root.mafdmp \
	    ${simOutDir}/sCow-sDog/sCow-sDog.maf ${simDir}/sCow-sDog/root.tmp.maf ${simOutDir}/sCow-sDog/root.maf >& $(basename $@).err

${simOutDir}/sH-sM-sR/sH-sM-sR.maf: ${simOutDir}/sMouse-sRat/sH-sM-sR.maf ${simDir}/simHuman/sH-sM-sR.maf
	@mkdir -p $(dir $@)
	${mafJoin} ${mafJoinOpts} -treelessRoot2=sH-sM-sR 'sH-sM-sR' \
	    -inMaf1Dump=${simOutDir}/sMouse-sRat/sH-sM-sR.mafdmp -inMaf2Dump=${simDir}/simHuman/sH-sM-sR.mafdmp -outMafDump=${simOutDir}/sH-sM-sR/sH-sM-sR.mafdmp \
	    ${simOutDir}/sMouse-sRat/sH-sM-sR.maf ${simDir}/simHuman/sH-sM-sR.maf ${simOutDir}/sH-sM-sR/sH-sM-sR.maf >& $(basename $@).err

${simOutDir}/sH-sM-sR/root.maf: ${simOutDir}/sH-sM-sR/sH-sM-sR.maf ${simDir}/sH-sM-sR/root.tmp.maf
	@mkdir -p $(dir $@)
	${mafJoin} ${mafJoinOpts} -treelessRoot2=root 'sH-sM-sR' \
	    -inMaf1Dump=${simOutDir}/sH-sM-sR/sH-sM-sR.mafdmp -inMaf2Dump=${simDir}/sH-sM-sR/root.tmp.mafdmp -outMafDump=${simOutDir}/sH-sM-sR/root.mafdmp \
	    ${simOutDir}/sH-sM-sR/sH-sM-sR.maf ${simDir}/sH-sM-sR/root.tmp.maf ${simOutDir}/sH-sM-sR/root.maf >& $(basename $@).err

${simOutDir}/root/root.maf: ${simOutDir}/sH-sM-sR/root.maf ${simOutDir}/sCow-sDog/root.maf
	@mkdir -p $(dir $@)
	${mafJoin} ${mafJoinOpts} 'root' \
	    -inMaf1Dump=${simOutDir}/sH-sM-sR/root.mafdmp -inMaf2Dump=${simOutDir}/sCow-sDog/root.mafdmp -outMafDump=${simOutDir}/root/root.mafdmp \
	    ${simOutDir}/sH-sM-sR/root.maf ${simOutDir}/sCow-sDog/root.maf ${simOutDir}/root/root.maf >& $(basename $@).err


compareBaseNames = \
    simRat/sMouse-sRat \
    simMouse/sMouse-sRat \
    sMouse-sRat/sH-sM-sR.tmp \
    simDog/sCow-sDog \
    simCow/sCow-sDog \
    sCow-sDog/root.tmp \
    simHuman/sH-sM-sR \
    sH-sM-sR/root.tmp

include runCommon.mk

