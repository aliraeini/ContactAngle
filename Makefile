msSrc  ?= ${abspath ..}
msInst ?= ${abspath ../..}

ifneq "${OPT}" ".exe"
all:  ;
	$(MAKE) -f  ${msSrc}/script/Makefile.in  recurseMake USE_msRecurse=1
	cp AllRunContAngle  ${msInst}/bin/foamx4m/
else
all: ; @echo skipping ${CURDIR}
endif

clean:; $(MAKE) -f  ${msSrc}/script/Makefile.in  recurseClean USE_msRecurse=1

tsts=test.sh
USE_msTEST=1
include  ${msSrc}/script/Makefile.in
