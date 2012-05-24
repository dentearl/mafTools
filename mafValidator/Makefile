include ../include.mk
binPath = ../bin

progs =  $(foreach f,mafValidator.py, ${binPath}/$f)

.PHONY: all clean test

all: ${progs}

${binPath}/%.py : src/%.py
	@mkdir -p $(dir $@)
	cp $< $@.tmp
	chmod +x $@.tmp
	mv $@.tmp $@

test :
	python src/test.mafValidator.py -v && rmdir tempTestDir

clean :
	rm -f ${progs}
