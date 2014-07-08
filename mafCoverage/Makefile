# Copyright (C) 2009-2013 by
# Dent Earl (dearl@soe.ucsc.edu, dentearl@gmail.com)
# Benedict Paten (benedict@soe.ucsc.edu, benedictpaten@gmail.com)
# Mark Diekhans (markd@soe.ucsc.edu)
# ... and other members of the Reconstruction Team of David Haussler's
# lab (BME Dept. UCSC).
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

include ../inc/common.mk
SHELL:=/bin/bash
bin = ../bin
inc = ../inc
lib = ../lib
PROGS = mafCoverage
dependencies = ${inc}/common.h ${inc}/sharedMaf.h ${lib}/common.c ${lib}/sharedMaf.c $(wildcard ${sonLibPath}/*) ${sonLibPath}/sonLib.a src/allTests.c
extraAPI := ${lib}/common.o ${lib}/sharedMaf.o ../external/CuTest.a src/mafCoverageAPI.o ${sonLibPath}/sonLib.a src/buildVersion.o
testAPI := test/sharedMaf.o test/common.o ../external/CuTest.a test/mafCoverageAPI.o ${sonLibPath}/sonLib.a test/buildVersion.o
testObjects := test/test.mafCoverageAPI.o
sources := src/mafCoverage.c src/mafCoverage.h

.PHONY: all clean test buildVersion

all: buildVersion $(foreach f,${PROGS}, ${bin}/$f)
buildVersion: src/buildVersion.c
src/buildVersion.c: ${sources} ${dependencies}
	@python ../lib/createVersionSources.py

../lib/%.o: ../lib/%.c ../inc/%.h
	cd ../lib/ && make

${bin}/mafCoverage: src/mafCoverage.c ${dependencies} ${extraAPI}
	mkdir -p $(dir $@)
	${cxx} $< ${extraAPI} -o $@.tmp ${cflags} ${lm}
	mv $@.tmp $@
%.o: %.c %.h
	${cxx} -c $< -o $@.tmp ${cflags}
	mv $@.tmp $@
%/mafCoverageAPI.o: src/mafCoverageAPI.c src/mafCoverageAPI.h
	${cxx} -c $< -o $@.tmp ${cflags}
	mv $@.tmp $@

test: buildVersion test/allTests test/mafCoverage
	./test/allTests && python2.7 src/test.mafCoverage.py --verbose  && rm -rf ./test/ && rmdir ./tempTestDir
test/allTests: src/allTests.c ${testAPI} ${testObjects} ${sonLibPath}/sonLib.a
	mkdir -p $(dir $@)
	${cxx} $^ -o $@.tmp ${testFlags} ${lm}
	mv $@.tmp $@
test/mafCoverage: src/mafCoverage.c ${dependencies} ${testAPI}
	mkdir -p $(dir $@)
	${cxx} $< ${testAPI} -o $@.tmp ${testFlags} ${lm}
	mv $@.tmp $@
test/%.o: ${lib}/%.c ${inc}/%.h
	mkdir -p $(dir $@)
	${cxx} -c $< -o $@.tmp ${testFlags}
	mv $@.tmp $@
test/test.mafCoverageAPI.o: src/test.mafCoverageAPI.c src/test.mafCoverageAPI.h test/mafCoverageAPI.o
	mkdir -p $(dir $@)
	${cxx} -c $< -o $@.tmp ${testFlags}
	mv $@.tmp $@
test/%.o: src/%.c src/%.h ${sonLibPath}/sonLib.a
	mkdir -p $(dir $@)
	${cxx} -c $< -o $@.tmp ${testFlags}
	mv $@.tmp $@

clean:
	rm -rf $(foreach f,${PROGS}, ${bin}/$f) src/*.o test/ src/buildVersion.c src/buildVersion.h
../external/CuTest.a: ../external/CuTest.c ../external/CuTest.h
	${cxx} -c $< ${cflags}
	ar rc CuTest.a CuTest.o
	ranlib CuTest.a
	rm -f CuTest.o
	mv CuTest.a $@
