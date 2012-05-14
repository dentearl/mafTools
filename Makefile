include include.mk
SHELL:=/bin/bash -e
export SHELLOPTS=pipefail

##############################
# These modules are dependent and are
# only included if their depedencies exist!
ifeq ($(wildcard ../sonLib/Makefile),)
	Comparator = 
	TransitiveClosure = 
$(warning Because dependency ../sonLib is missing mafComparator and mafTransitiveClosure will not be built / tested / cleaned. See README.md for information about dependencies.)
else
	Comparator = mafComparator
ifeq ($(wildcard ../pinchesAndCacti/Makefile),)
	TransitiveClosure = 
$(warning Because dependency ../pinchesAndCacti is missing mafTransitiveClosure will not be built / tested / cleaned. See README.md for information about dependencies.)
else
	TransitiveClosure = mafTransitiveClosure
endif
endif
##############################
dependentModules= ${Comparator} ${TransitiveClosure}

modules = ${dependentModules} mafCoveragePickles mafValidator mafBlockFinder mafBlockExtractor mafBlockSorter mafBlockDuplicateFilter include 

.PHONY: all %.all clean %.clean test %.test

all: ${modules:%=%.all}

%.all:
	cd $* && make all

clean: ${modules:%=%.clean}

%.clean:
	cd $* && make clean

test : ${modules:%=%.test} ${Warnings:%=%.warn}

%.test:
	cd $* && make test
