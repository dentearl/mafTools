include ../include.mk
binPath := ../bin
srcPath := src
.PHONY := all test clean
scripts := mafCoveragePickleCreator.py mafCoveragePicklePlotter.py mafCoveragePickleSubsetExtractor.py mafCoveragePickleGapAnalysis.py

${binPath}/%.py: ${srcPath}/%.py
	mkdir -p $(dir $@)
	cp $< $@.tmp
	chmod 755 $@.tmp
	mv $@.tmp $@

all: $(foreach f,${scripts},${binPath}/$f)

test:
	python2.6 ${srcPath}/test.mafCoveragePickleCreator.py --verbose
	python2.6 ${srcPath}/test.mafCoveragePickleGapAnalysis.py --verbose

clean:
	rm -f $(foreach f,${scripts},${binPath}/$f)
