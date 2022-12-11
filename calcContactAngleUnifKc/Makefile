
all: msApp

msSrc ?= ${abspath ../..}

USE_MINIFOAMX3=1
OFCs=$(shell find | grep -F .C | sed s:./::)
OFincs = surfMesh
OFlibs = surfMesh

msCFLAGS += -IincompressibleTwoPhaseMixture
USE_msMAKE=1
include  ${msSrc}/script/Makefile.in
