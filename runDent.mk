#
# This target  
#

# simDir copied from /hive/users/dearl/simTreeWorking/MARK_MAFJOIN
simDir = MARK_MAFJOIN
export PATH := ../../../bin:/hive/groups/recon/local/bin:${PATH}
export PYTHONPATH = ../../../src

maxBlkWidth=10000

step1Mafs = \
	${simDir}/sHuman-sChimp/sHuman-sChimp.maf \
	${simDir}/sHuman-sChimp/sG-sH-sC.maf

step2Mafs = \
	${simDir}/sG-sH-sC/sG-sH-sC.maf \
	${simDir}/sG-sH-sC/root.maf
step3Mafs = \
	${simDir}/root/root.maf

createdMafs = ${step1Mafs} ${step2Mafs} ${step3Mafs}

mafJoin = time mafJoin

all: ${createdMafs}
## step 1
${simDir}/sHuman-sChimp/sHuman-sChimp.maf: ${simDir}/simChimp/sHuman-sChimp.maf ${simDir}/simHuman/sHuman-sChimp.maf
	${mafJoin} -treelessRoot1=sHuman-sChimp -treelessRoot2=sHuman-sChimp 'sHuman-sChimp' ${simDir}/simChimp/sHuman-sChimp.maf ${simDir}/simHuman/sHuman-sChimp.maf ${simDir}/sHuman-sChimp/sHuman-sChimp.maf

${simDir}/sHuman-sChimp/sG-sH-sC.maf:  ${simDir}/sHuman-sChimp/sHuman-sChimp.maf ${simDir}/sHuman-sChimp/sG-sH-sC.tmp.maf 
	${mafJoin} -treelessRoot2=sG-sH-sC 'sHuman-sChimp' ${simDir}/sHuman-sChimp/sHuman-sChimp.maf ${simDir}/sHuman-sChimp/sG-sH-sC.tmp.maf ${simDir}/sHuman-sChimp/sG-sH-sC.maf

## step 2
${simDir}/sG-sH-sC/sG-sH-sC.maf:  ${simDir}/sHuman-sChimp/sG-sH-sC.maf ${simDir}/simGorilla/sG-sH-sC.maf
	${mafJoin} -treelessRoot2=sG-sH-sC 'sG-sH-sC' ${simDir}/sHuman-sChimp/sG-sH-sC.maf ${simDir}/simGorilla/sG-sH-sC.maf ${simDir}/sG-sH-sC/sG-sH-sC.maf

${simDir}/sG-sH-sC/root.maf:  ${simDir}/sG-sH-sC/sG-sH-sC.maf ${simDir}/sG-sH-sC/root.tmp.maf
	${mafJoin} -maxBlkWidth=${maxBlkWidth} -treelessRoot2=root 'sG-sH-sC' ${simDir}/sG-sH-sC/sG-sH-sC.maf ${simDir}/sG-sH-sC/root.tmp.maf ${simDir}/sG-sH-sC/root.maf

## step 3
${simDir}/root/root.maf: ${simDir}/simOrang/root.maf ${simDir}/sG-sH-sC/root.maf 
	${mafJoin} -maxBlkWidth=${maxBlkWidth} -treelessRoot1=root 'root' ${simDir}/simOrang/root.maf ${simDir}/sG-sH-sC/root.maf ${simDir}/root/root.maf


clean:
	rm -f ${createdMafs}

##
# this part used to discover the next commands needed give the current
# state.  We then add the result to the make file and do another round
##
getNextCmds:
	python ../simControl/simCtrl_postSimMAFextractor.py --simDir ${simDir} --debug --mergeStep

