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

include ../include.mk
SHELL:=/bin/bash
bin = ../bin
PROGS = mafPairCoverage
dependencies = $(wildcard ../include/common.*) $(wildcard ../include/sharedMaf.*) $(wildcard ${sonLibPath}/*) ${sonLibPath}/sonLib.a src/allTests.c
extraAPI := ../include/common.o ../include/sharedMaf.o ../external/CuTest.a src/mafPairCoverageAPI.o ${sonLibPath}/sonLib.a src/buildVersion.o 
testAPI := test/sharedMaf.o test/common.o ../external/CuTest.a test/mafPairCoverageAPI.o ${sonLibPath}/sonLib.a test/buildVersion.o
testObjects := test/test.mafPairCoverageAPI.o 
sources := src/mafPairCoverage.c src/mafPairCoverage.h

.PHONY: all clean test buildVersion

all: buildVersion $(foreach f,${PROGS}, ${bin}/$f)
buildVersion: src/buildVersion.c
src/buildVersion.c: ${sources} ${dependencies}
	@python ../include/createVersionSources.py

${bin}/mafPairCoverage: src/mafPairCoverage.c ${dependencies} ${extraAPI}
	mkdir -p $(dir $@)
	${cxx} $< ${extraAPI} -o $@.tmp ${cflags}
	mv $@.tmp $@
%.o: %.c %.h
	${cxx} -c $< -o $@.tmp ${cflags}
	mv $@.tmp $@
%/mafPairCoverageAPI.o: src/mafPairCoverageAPI.c src/mafPairCoverageAPI.h
	${cxx} -c $< -o $@.tmp ${cflags}
	mv $@.tmp $@

test: buildVersion test/allTests test/mafPairCoverage
	./test/allTests && python2.7 src/test.mafPairCoverage.py --verbose && rm -rf ./test/ && rmdir ./tempTestDir
test/allTests: src/allTests.c ${testAPI} ${testObjects} ${sonLibPath}/sonLib.a
	mkdir -p $(dir $@)
	${cxx} $^ -o $@.tmp ${testFlags}
	mv $@.tmp $@
test/mafPairCoverage: src/mafPairCoverage.c ${dependencies} ${testAPI}
	mkdir -p $(dir $@)
	${cxx} $< ${testAPI} -o $@.tmp ${testFlags}
	mv $@.tmp $@
test/%.o: ../include/%.c ../include/%.h
	mkdir -p $(dir $@)
	${cxx} -c $< -o $@.tmp ${testFlags}
	mv $@.tmp $@
test/test.mafPairCoverageAPI.o: src/test.mafPairCoverageAPI.c src/test.mafPairCoverageAPI.h test/mafPairCoverageAPI.o
	mkdir -p $(dir $@)
	${cxx} -c $< -o $@.tmp ${testFlags}
	mv $@.tmp $@
test/%.o: src/%.c src/%.h ${sonLibPath}/sonLib.a
	mkdir -p $(dir $@)
	${cxx} -c $< -o $@.tmp ${testFlags}
	mv $@.tmp $@

clean:
	rm -rf $(foreach f,${PROGS}, ${bin}/$f) src/*.o test/ 
../external/CuTest.a: ../external/CuTest.c ../external/CuTest.h
	${cxx} -c $< ${cflags}
	ar rc CuTest.a CuTest.o
	ranlib CuTest.a
	rm -f CuTest.o
	mv CuTest.a $@
