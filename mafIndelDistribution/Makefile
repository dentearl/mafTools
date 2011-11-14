include ../include.mk
binPath := ../bin
srcPath := src/
.PHONY := all test clean

${binPath}/mafIndelDistribution.py: ${srcPath}/mafIndelDistribution.py
	mkdir -p $(dir $@)
	cp $< $@.tmp
	chmod 755 $@.tmp
	mv $@.tmp $@

all: ${binPath}/mafIndelDistribution.py

test:
	python2.6 ${srcPath}/test.mafIndelDistribution.py --verbose

clean:
	rm -f ${binPath}/mafIndelDistribution.py ${srcPath}/mafIndelDistribution.pyc
