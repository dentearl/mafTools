include ../include.mk
binPath := ../bin
srcPath := src
.PHONY := all test clean
scripts := mafCoveragePickleCreator.py mafCoveragePicklePlotter.py mafCoveragePickleSubsetExtractor.py mafCoveragePickleAnalysis.py

${binPath}/%.py: ${srcPath}/%.py
	mkdir -p $(dir $@)
	cp $< $@.tmp
	chmod 755 $@.tmp
	mv $@.tmp $@

all: $(foreach f,${scripts},${binPath}/$f)

test:
	python2.7 ${srcPath}/test.mafCoveragePickleCreator.py --verbose
	python2.7 ${srcPath}/test.mafCoveragePickleAnalysis.py --verbose

clean:
	rm -f $(foreach f,${scripts},${binPath}/$f)
