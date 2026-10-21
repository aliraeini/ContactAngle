

all: msApp

msSrc ?= $(abspath ../../../src)

USE_MINIFOAMX3=1
OFCs=$(shell find | grep -F .C | sed s:./::)
OFincs =   surfMesh meshTools
OFlibs =   surfMesh meshTools

USE_ZLIB=1
USE_TIFF=1
USE_VOX_IMG=1
USE_POCKETPY=1
USE_msMAKE=1
msCFLAGS += -I${msSrc}/../msmkit/voxutils -I${msSrc}/../msmkit/primitives -I${msSrc}/../msmkit/optim
include  ${msSrc}/script/Makefile.in
