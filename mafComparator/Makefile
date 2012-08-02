# Copyright (C) 2009-2012 by 
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
binPath = ../bin
extraAPI = src/cString.c ../include/sharedMaf.o ../external/CuTest.a ../include/common.o src/comparatorRandom.o src/comparatorAPI.o ${sonLibPath}/sonLib.a src/buildVersion.o
testAPI = src/cString.c test/sharedMaf.o ../external/CuTest.a test/common.o test/comparatorRandom.o test/comparatorAPI.o ${sonLibPath}/sonLib.a test/buildVersion.o
progs =  $(foreach f, mafComparator mafPairCounter, ${binPath}/$f)
testObjects = test/test.comparatorAPI.o test/test.comparatorRandom.o

.PHONY: all clean test buildVersion

all: buildVersion ${progs}
buildVersion:
	@python src/createVersionSources.py

${binPath}/%: src/%.c ${extraAPI} $(wildcard src/*.h)
	@mkdir -p $(dir $@)
	${cxx} -o $@.tmp $^ ${cflags}
	mv $@.tmp $@

test/%: src/%.c ${testAPI} $(wildcard src/*.h)
	@mkdir -p $(dir $@)
	${cxx} -o $@.tmp $^ ${testFlags}
	mv $@.tmp $@

${binPath}/%.py: src/%.py
	@mkdir -p $(dir $@)
	cp $< $@.tmp
	chmod +x $@.tmp
	mv $@.tmp $@

%.o: %.c %.h
	${cxx} -c $< -o $@.tmp ${cflags}
	mv $@.tmp $@
test/%.o: ../include/%.c ../include/%.h
	mkdir -p $(dir $@)
	${cxx} -c $< -o $@.tmp ${testFlags}
	mv $@.tmp $@
test/test.comparatorAPI.o: src/test.comparatorAPI.c src/test.comparatorAPI.h test/comparatorAPI.o
	mkdir -p $(dir $@)
	${cxx} -c $< -o $@.tmp ${testFlags}
	mv $@.tmp $@
test/%.o: src/%.c src/%.h ${sonLibPath}/sonLib.a
	mkdir -p $(dir $@)
	${cxx} -c $< -o $@.tmp ${testFlags}
	mv $@.tmp $@

../external/CuTest.a: ../external/CuTest.c ../external/CuTest.h
	${cxx} -c $< ${cflags}
	ar rc CuTest.a CuTest.o
	ranlib CuTest.a
	rm -f CuTest.o
	mv CuTest.a $@

test: test/allTests test/mafComparator test/testRand
	./test/allTests && python2.7 src/test.mafComparator.py --verbose  && rm -rf ./test/ && rmdir ./tempTestDir

test/allTests: src/allTests.c ${testAPI} ${testObjects} ${sonLibPath}/sonLib.a
	mkdir -p $(dir $@)
	${cxx} $^ -o $@.tmp ${testFlags}
	mv $@.tmp $@

# to actually use the testRand program, comment out the rm -rf on the "test:" rule and run "make test", 
# then you may run test/testRand 
test/testRand: src/testRand.c ${testAPI} ${sonLibPath}/sonLib.a 
	mkdir -p $(dir $@)
	${cxx} $^ -o $@.tmp ${testFlags}
	mv $@.tmp $@

clean:
	rm -f *.o ${progs} src/*.o && rm -rf ./test/
