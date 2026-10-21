
export msSrc ?= $(abspath ../../src)
MS_INST ?= ../../..


all:  ;
	$(MAKE) -f  ${msSrc}/script/Makefile.in  recurseMake USE_msRecurse=1
	mkdir -p ${MS_INST}/bin/ContAngle
	cp -a tutorial  ${MS_INST}/bin/ContAngle/
	cp AllRunContAngle  ${MS_INST}/bin/

clean:; $(MAKE) -f  ${msSrc}/script/Makefile.in  recurseClean USE_msRecurse=1

tsts=test.py
USE_msTEST=1
export msSrc

include  ${msSrc}/script/Makefile.in
