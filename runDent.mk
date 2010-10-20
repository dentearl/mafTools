#
# This target  
#

# simDir copied from /hive/users/dearl/simTreeWorking/MARK_MAFJOIN
simDir = MARK_MAFJOIN
export PATH := ../../../bin:/hive/groups/recon/local/bin:${PATH}
export PYTHONPATH = ../../../src

all:	${simDir}/sHuman-sChimp/sHuman-sChimp.maf \
	${simDir}/sHuman-sChimp/sG-sH-sC.maf
## step  one
${simDir}/sHuman-sChimp/sHuman-sChimp.maf: ${simDir}/simChimp/sHuman-sChimp.maf ${simDir}/simHuman/sHuman-sChimp.maf
	mafJoin -treelessRoot1=sHuman-sChimp -treelessRoot2=sHuman-sChimp 'sHuman-sChimp' ${simDir}/simChimp/sHuman-sChimp.maf ${simDir}/simHuman/sHuman-sChimp.maf ${simDir}/sHuman-sChimp/sHuman-sChimp.maf

${simDir}/sHuman-sChimp/sG-sH-sC.maf:  ${simDir}/sHuman-sChimp/sHuman-sChimp.maf ${simDir}/sHuman-sChimp/sG-sH-sC.tmp.maf 
	mafJoin -treelessRoot2=sG-sH-sC 'sHuman-sChimp' ${simDir}/sHuman-sChimp/sHuman-sChimp.maf ${simDir}/sHuman-sChimp/sG-sH-sC.tmp.maf ${simDir}/sHuman-sChimp/sG-sH-sC.maf



##
# this part used to discover the next commands needed give the current
# state.  We then add the result to the make file and do another round
##
getNextCmds:
	python simControl/simCtrl_postSimMAFextractor.py --simDir ${simDir} --debug --mergeStep
