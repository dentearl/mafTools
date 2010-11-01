#
# override simOutDir to save experiments
#

# simDir copied from /hive/users/dearl/simTreeWorking/MARK_BENCODE_2P/
simRoot = evolver/primate
simDir = ${simRoot}/simulation
simOutDir = ${simRoot}/output/expr1

export PATH := ../../../bin:/hive/groups/recon/local/bin:${PATH}
export PYTHONPATH = ../../../src

mafJoinOpts = -maxBlkWidth=10000
#mafJoinOpts = -maxInputBlkWidth=2500  # doesn't work very well

createdMafs = \
	${simOutDir}/sHuman-sChimp/sHuman-sChimp.maf \
	${simOutDir}/sHuman-sChimp/sG-sH-sC.maf \
	${simOutDir}/sG-sH-sC/sG-sH-sC.maf \
	${simOutDir}/sG-sH-sC/root.maf \
	${simOutDir}/root/root.maf


mafJoin = time -p mafJoin

join: ${createdMafs}
${simOutDir}/sHuman-sChimp/sHuman-sChimp.maf: ${simDir}/simChimp/sHuman-sChimp.maf ${simDir}/simHuman/sHuman-sChimp.maf
	@mkdir -p $(dir $@)
	${mafJoin} ${mafJoinOpts} -treelessRoot1=sHuman-sChimp -treelessRoot2=sHuman-sChimp sHuman-sChimp ${simDir}/simChimp/sHuman-sChimp.maf ${simDir}/simHuman/sHuman-sChimp.maf ${simOutDir}/sHuman-sChimp/sHuman-sChimp.maf >& $(basename $@).err

${simOutDir}/sHuman-sChimp/sG-sH-sC.maf:  ${simOutDir}/sHuman-sChimp/sHuman-sChimp.maf ${simDir}/sHuman-sChimp/sG-sH-sC.tmp.maf 
	@mkdir -p $(dir $@)
	${mafJoin} ${mafJoinOpts} -treelessRoot2=sG-sH-sC sHuman-sChimp ${simOutDir}/sHuman-sChimp/sHuman-sChimp.maf ${simDir}/sHuman-sChimp/sG-sH-sC.tmp.maf ${simOutDir}/sHuman-sChimp/sG-sH-sC.maf >& $(basename $@).err

${simOutDir}/sG-sH-sC/sG-sH-sC.maf: ${simOutDir}/sHuman-sChimp/sG-sH-sC.maf ${simDir}/simGorilla/sG-sH-sC.maf
	@mkdir -p $(dir $@)
	${mafJoin} ${mafJoinOpts} -treelessRoot2=sG-sH-sC sG-sH-sC ${simOutDir}/sHuman-sChimp/sG-sH-sC.maf ${simDir}/simGorilla/sG-sH-sC.maf ${simOutDir}/sG-sH-sC/sG-sH-sC.maf >& $(basename $@).err

${simOutDir}/sG-sH-sC/root.maf:  ${simOutDir}/sG-sH-sC/sG-sH-sC.maf ${simDir}/sG-sH-sC/root.tmp.maf
	@mkdir -p $(dir $@)
	${mafJoin} ${mafJoinOpts} -treelessRoot2=root sG-sH-sC ${simOutDir}/sG-sH-sC/sG-sH-sC.maf ${simDir}/sG-sH-sC/root.tmp.maf ${simOutDir}/sG-sH-sC/root.maf >& $(basename $@).err

${simOutDir}/root/root.maf: ${simDir}/simOrang/root.maf ${simOutDir}/sG-sH-sC/root.maf 
	@mkdir -p $(dir $@)
	${mafJoin} ${mafJoinOpts} -treelessRoot1=root root ${simDir}/simOrang/root.maf ${simOutDir}/sG-sH-sC/root.maf ${simOutDir}/root/root.maf >& $(basename $@).err


compareBaseNames = \
    simChimp/sHuman-sChimp \
    simHuman/sHuman-sChimp \
    sHuman-sChimp/sG-sH-sC.tmp \
    simGorilla/sG-sH-sC \
    sG-sH-sC/root.tmp \
    simOrang/root

include runCommon.mk

##
# this part used to discover the next commands needed give the current
# state.  We then add the result to the make file and do another round
##
getNextCmds:
	python ../simControl/simCtrl_postSimMAFextractor.py --simDir ${simDir} --debug --mergeStep

