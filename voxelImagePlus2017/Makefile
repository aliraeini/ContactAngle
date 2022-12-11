all_: ; @echo not maintained -- build manually by running 'make all' in voxelImagePlus2017

all: writeRaCLInds writeRaPoreInds writepAreaPoreInds writeCAThroatInds writeRcThroatInds writeCAPoreInds writeRcPoreInds writeCACornerInds writeRcCornerInds
# warning type clean, dependencies are not handelled yet

myBinDir ?= ../../../bin




clean:
	rm  -f *.o


writepAreaPoreInds.o:	writepAreaPoreInds.cpp
	g++ -O3 -Wall -std=c++0x -pedantic -c writepAreaPoreInds.cpp
writepAreaPoreInds:	writepAreaPoreInds.o
	g++ -O3 -Wall    writepAreaPoreInds.o  -o $(myBinDir)/writepAreaPoreInds

writeCAPoreInds.o:	writeCAPoreInds.cpp
	g++ -O3 -Wall -std=c++0x -pedantic -c writeCAPoreInds.cpp
writeCAPoreInds:	writeCAPoreInds.o
	g++ -O3 -Wall    writeCAPoreInds.o  -o $(myBinDir)/writeCAPoreInds

writeRcPoreInds.o:	writeRcPoreInds.cpp
	g++ -O3 -Wall -std=c++0x -pedantic -c writeRcPoreInds.cpp
writeRcPoreInds:	writeRcPoreInds.o
	g++ -O3 -Wall    writeRcPoreInds.o  -o $(myBinDir)/writeRcPoreInds

writeRaPoreInds.o:	writeRaPoreInds.cpp
	g++ -O3 -Wall -std=c++0x -pedantic -c writeRaPoreInds.cpp
writeRaPoreInds:	writeRaPoreInds.o
	g++ -O3 -Wall    writeRaPoreInds.o  -o $(myBinDir)/writeRaPoreInds


writeRaCLInds.o:	writeRaCLInds.cpp
	g++ -O3 -Wall -std=c++0x -pedantic -c writeRaCLInds.cpp
writeRaCLInds:	writeRaCLInds.o
	g++ -O3 -Wall    writeRaCLInds.o  -o $(myBinDir)/writeRaCLInds


writeCACornerInds.o:	writeCACornerInds.cpp
	g++ -O3 -Wall -std=c++0x -pedantic -c writeCACornerInds.cpp
writeCACornerInds:	writeCACornerInds.o
	g++ -O3 -Wall    writeCACornerInds.o  -o $(myBinDir)/writeCACornerInds

writeRcCornerInds.o:	writeRcCornerInds.cpp
	g++ -O3 -Wall -std=c++0x -pedantic -c writeRcCornerInds.cpp
writeRcCornerInds:	writeRcCornerInds.o
	g++ -O3 -Wall    writeRcCornerInds.o  -o $(myBinDir)/writeRcCornerInds



writeCAThroatInds.o:	writeCAThroatInds.cpp
	g++ -O3 -Wall -std=c++0x -pedantic -c writeCAThroatInds.cpp
writeCAThroatInds:	writeCAThroatInds.o
	g++ -O3 -Wall    writeCAThroatInds.o  -o $(myBinDir)/writeCAThroatInds

writeRcThroatInds.o:	writeRcThroatInds.cpp
	g++ -O3 -Wall -std=c++0x -pedantic -c writeRcThroatInds.cpp
writeRcThroatInds:	writeRcThroatInds.o
	g++ -O3 -Wall    writeRcThroatInds.o  -o $(myBinDir)/writeRcThroatInds
