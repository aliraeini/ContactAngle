# warning type clean, dependencies are not handelled yet
msSrc ?= $(abspath ../../../src)

mains =   writeCAPoreInds  writepAreaPoreInds
all: $(mains)
$(info mains: $(mains),   tsts: $(tsts) )

USE_ZLIB=1
USE_TIFF=1
USE_VOX_IMG=1
MS_VERSION=2
USE_POCKETPY=1
USE_OMP=1
USE_SINGLECPP=1
USE_msTEST=1
msCFLAGS += -I${msSrc}/../msmkit/voxutils -I${msSrc}/../msmkit/primitives -I${msSrc}/../msmkit/optim
include  ${msSrc}/script/Makefile.in
