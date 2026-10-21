/**-------------------------------------------------------------------------*\
# writeCAPoreInds: Map curvatures and areas to pore-network

Developed by:
 - Ahmed AlRatrout (2014-2018) as a part his PhD at Imperial College
   London, under the supervision of Branko Bijeljic and Martin Blunt.

 - Ali Q Raeini (2021)
   write pore-averaged contact angles -- file outPoresAvgCAngle.txt

For contact and references visit:
  https://www.imperial.ac.uk/earth-science/research/research-groups/pore-scale-modelling

Copyright:
  You can redistribute this code and/or modify this code under the
  terms of the GNU General Public License (GPL) as published by the
  Free Software Foundation, either version 3 of the License, or (at
  your option) any later version. see <http://www.gnu.org/licenses/>.

\*------------------------------------------------------------------------**/

#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>

#include <cassert>

#include "voxelImageI.h"


int usage()
{
	std::cout<<
		"Usage:"
		"  writeAreaPoreInds  Image  \n"
		" or \n"
		"  writeAreaPoreInds  Image  contactAngles_x.tsv  \n"
		"where 'Image'_VElems.mhd should be present relative to current directory.\n"
		"Outputs: \n"
		"  outCAnglePoreIndex.tsv # contact angle and pore index of contact line points\n"
		"  outPoresAvgCAngle.tsv  # average contact angles of pores\n"
		<< std::endl;
	return -1;
}


int main(int argc, char *argv[])
{
	if(argc<2) return usage();


	VoxelImageT<int> vimage(std::string(argv[1])+"_VElems.mhd", readOpt::justRead); // (fileName_VElems.mhd)
	ensure(vimage.nx(), "cannot read image "+std::string(argv[1]),-1);

	int3 nn=vimage.size3();
	dbl3 X0=vimage.X0();
	dbl3 dx=vimage.dx();

	std::string filNam = argc>2 ? argv[2] : "point_Kc_A.tsv";  // file name
	std::cout<<"\nOpenning CA file: "<<filNam<<std::endl;
	std::ifstream infile(filNam); ///pA_x.txt
	ensure(infile, "cannot read "+filNam,-1);


	int  nPors=accumulate(vimage,std::max,0);

	std::string tmpc;
	getline(infile,tmpc);// skip first line
	dbls Areas(nPors,0.), Kcs(nPors,0.),  sumWs(nPors,0.);


	int countLines = 0;
	std::ofstream outCA("outAsKcPoreIndex.tsv");
	outCA<<"Number of pores = \t"<<nPors<<std::endl;
	outCA << "X"<<"\t"<<"Y"<<"\t"<<"Z"<<"\t"<<"cosCA"<<"\t"<<"ContAngle"<<"\t"<<"PoreI"<<"\n";
	while(infile.good()) {
		++countLines;
		double xx,yy,zz, pA,Kc;
		infile >> xx>>yy>>zz>>pA>>Kc;

		int ii = std::max(1,std::min(int((xx/dx.x)+1.49-X0.x),nn.x-2)); // VElem has +-one extra voxels
		int jj = std::max(1,std::min(int((yy/dx.y)+1.49-X0.y),nn.y-2));
		int kk = std::max(1,std::min(int((zz/dx.z)+1.49-X0.z),nn.z-2));

		std::cout<< " " <<kk<< " " <<jj<< " " <<ii<< " : " <<std::endl;
		int EI = vimage(ii,jj,kk);

		if(EI<=0)	EI = maxNei(vimage,ii,jj,kk,-1,1);

		if(EI>0)
		{
			outCA <<xx<<"	"<<yy<<"	"<<zz<<"	"<<pA<<"	"<<Kc<<"	"<<EI<<"\n";
			Areas[EI] += pA;  Kcs[EI] += Kc;  sumWs[EI] += 1.;
		}
	}
	std::cout<< "n_Faces = " << countLines <<", valid:"<<int(sum(sumWs)/countLines*1000)/10<< "%\n :/ \n\n" << std::endl;

	std::ofstream outPorCA("outPoresAvgAreaKc.tsv");
	outPorCA<<"PoreIndx 	 area 	curvature	 Nfaces";
	for (int EI=0; EI<nPors; ++EI)
		if(sumWs[EI]>0.)
			outPorCA <<Areas[EI]<<'\t'<<Kcs[EI]/(sumWs[EI]+1.0e-36)<<'\t'<<sumWs[EI]<<'\n';

	std::cout<< "\n ***** END ***** \n\n\n" << std::endl;

	return 0;
}
