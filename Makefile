include include.mk
SHELL:=/bin/bash -e
export SHELLOPTS=pipefail

modules= mafComparator

.PHONY: all %.all clean %.clean test %.test

all: ${modules:%=all.%}

all.%:
	cd $* && make all

clean: ${modules:%=clean.%}

clean.%:
	cd $* && make clean

test : ${modules:%=test.%}

test.%:
	cd $* && make test
