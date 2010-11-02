#
# common targets and rules
#
compare: $(compareBaseNames:%=${simOutDir}/%.cmpout) $(compareBaseNames:%=${simOutDir}/%.probcmpout)

${simOutDir}/%.cmpout: ${simDir}/%.maf ${simOutDir}/root/root.maf
	@mkdir -p $(dir $@)
	eval_MAFComparator --mAFFile1 $< --mAFFile2 ${simOutDir}/root/root.maf --outputFile $(basename $@).xml  --ultraVerbose >& $(basename $@).cmpout.tmp
	mv -f $(basename $@).cmpout.tmp $(basename $@).cmpout

${simOutDir}/%.probcmpout: ${simOutDir}/%.cmpout
	awk '!/Missing:/ || $$2!=$$4' $< >$@

clean:
	rm -f ${simOutDir}

