/*-------------------------------------------------------------------------*\
You can redistribute this code and/or modify this code under the
terms of the GNU General Public License (GPL) as published by the
Free Software Foundation, either version 3 of the License, or (at
your option) any later version. see <http://www.gnu.org/licenses/>.


The code has been developed by Ahmed AlRatrout as a part his PhD
at Imperial College London, under the supervision of Branko Bijeljic
and Martin Blunt.


For further information please contact us by email:
Ahmed AlRarout:  a.alratrout14@imperial.ac.uk
Branko Bijeljic: b.bijeljic@imperial.ac.uk
Martin J Blunt:  m.blunt@imperial.ac.uk
\*-------------------------------------------------------------------------*/

#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>

#include <cassert>

#include "voxelImageI.h"


int usage()
{
	std::cout<<"Error: wrong number of arguments"
		<< std::endl;
	return -1;
}


int main(int argc, char *argv[])
{
	if(argc!=4) usage();


	VoxelImageT<int> vimage(argv[1], readOpt::justRead); // (fileName_VElems.mhd)
	ensure(vimage.nx(), "cannot read image "+std::string(argv[1]),-1);

	int3 nn=vimage.size3();
	dbl3 X0=vimage.X0();
	dbl3 dx=vimage.dx();

	std::string filNam = argv[2];  // file name
	std::cout<<"\nOpenning Ra file: "<<filNam<<std::endl;
	std::ifstream infile(filNam);  ///Ra_x.txt
	ensure(infile, "cannot read "+filNam,-1);


	std::cout<<"\n\n Openning VElems file: "<<argv[2]<<std::endl;
	std::ifstream inMSM(argv[3]); ///(fileName_ms.xmf)
	assert(inMSM);

	std::string tmpc;
	std::string ttmpc;
	std::string line;
	int  nPors;



	{ std::string tmpc; getline(infile,tmpc); }// skip first line
	for ( unsigned int i=0; i<5;i++)   inMSM>>ttmpc, std::cout<<" "<<ttmpc;
	inMSM>>nPors, std::cout<<"\n  nPores = "<<nPors<<std::endl;
	inMSM>>tmpc;

	std::vector<double> PRs(nPors+2,0.0);

	std::vector<int> i(nPors+2,0);
	std::vector<int> j(nPors+2,0);
	std::vector<int> k(nPors+2,0);


	for (int countPI = 2;countPI<nPors+2; ++countPI)
	{
		std::string lline;
		getline(inMSM, lline);
		char g;
		inMSM >> i[countPI]>>j[countPI]>>k[countPI]>>g>>PRs[countPI]>>tmpc;
		std::cout<< "\n";
		std::cout<< " " << countPI  << " " <<i[countPI]<< " " <<j[countPI]<< " " <<k[countPI]<<" "<<g<<" "<<PRs[countPI]<<"  \n";
	}

	int countLines = 0;
	std::ofstream outRa("outRa.txt");
	outRa<<"Number of pores = "<<nPors<<std::endl;
	while(infile.good()) {
		++countLines;
		double xx,yy,zz, Ra;
		infile >> xx>>yy>>zz>>Ra;

		int ii = std::max(1,std::min(int((xx/dx.x)+1.49-X0.x),nn.x-2)); // VElem has +-one extra voxels
		int jj = std::max(1,std::min(int((yy/dx.y)+1.49-X0.y),nn.y-2));
		int kk = std::max(1,std::min(int((zz/dx.z)+1.49-X0.z),nn.z-2));

		std::cout<< " " <<kk<< " " <<jj<< " " <<ii<< " : " <<std::endl;
		int EI = vimage(ii,jj,kk);

		if(EI<=0)	EI = maxNei(vimage,ii,jj,kk,-1,1);
		if(EI<=0)	EI = maxNei(vimage,ii,jj,kk,-2,2);

		if(PI>0)
		{
		  outRa << xx<<" "<<yy<<" "<<zz<<" "<<" "<<Ra<<" "<<PI<<" "<<(PI>=0 ? PRs[PI] : 0) <<" "<<i[PI]<<" "<<j[PI]<<" "<<k[PI]<<"\n";
		}
	}
	std::cout<< "\nTotal PI = " << nPors;

	std::cout<< "\nTotal Ra = " << countLines;


	std::cout<< "\n ***** END ***** \n\n\n" << std::endl;

	return 0;
}
