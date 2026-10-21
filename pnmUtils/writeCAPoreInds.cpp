/**-------------------------------------------------------------------------*\
# writeCAPoreInds: Map contact angles to pore-network


Developed by:

 - Ahmed AlRatrout 2014-2018
   as a part his PhD at Imperial College London, under the supervision
   of Branko Bijeljic and Martin Blunt.

 - Ali Q Raeini 2021
   write pore-averaged contact angles -- file outPoresAvgCAngle.txt

Contact and references:
 - visit https://www.imperial.ac.uk/earth-science/research/research-groups/pore-scale-modelling

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

using namespace std;

int usage()  {
	cout<<
		"Usage:"
		"  writeCAPoreInds  Image  \n"
		" or \n"
		"  writeCAPoreInds  Image  contactAngles_x.tsv  \n"
		"where 'Image'_VElems.mhd should be present relative to current directory.\n"
		"Outputs: \n"
		"  outCAnglePoreIndex.tsv # contact angle and pore index of contact line points\n"
		"  outPoresAvgCAngle.tsv  # average contact angles of pores\n"
		<< endl;
	return -1;
}
//#define DEVEL


int main(int argc, char *argv[])
{
	if(argc<2) return usage();

	string name = argv[1];
	if(hasExt(name,".mhd"))  name.resize(name.size()-4);
	if(hasExt(name,"_VElems"))  name.resize(name.size()-7);

	string CANam = argc>2 ? argv[2] : "contactAngles_x.txt";  // file name
	cout<<"\n Openning CA file: "<<CANam<<endl;
	ifstream CAFil(CANam);     ///contactAngles_x.tsv
	if(!CAFil) { CAFil.clear();  CAFil.open(name+"/"+CANam); } ///contactAngles_x.tsv
	ensure(CAFil, "cannot read "+CANam,-1);

	VoxelImageT<int> vimage(name+"_VElems.mhd", readOpt::justRead); // (fileName_VElems.mhd)
	ensure(vimage.nx(), "cannot read image "+string(argv[1]),-1);

	int3 nn=vimage.size3();
	dbl3 X0=vimage.X0();


	dbl3 dx=vimage.dx();
	ensure(mag(dx)<0.02, "ElementSize too large, assuming Unit is um");
	if(mag(dx)<0.02) {
		dx*=1e6; //COnAngle code use um as unit
		X0*=1e6; //COnAngle code use um as unit
	}

	ensure(nn.z>6, " !!! read image size "+_s(nn));
	ensure(mag(X0)<dx.x*3, "X0 (offset) should have been set to zero");


	int  nPors=accumulate(vimage,max,0)+1;

	string tmpc;
	getline(CAFil,tmpc);// skip first line
	dbls cosCAs(nPors,0.),  sumWs(nPors,0.);

	cout<<"dx: "<<dx<<",  N_img: "<<nn<<" "<<endl;
	int countLines = 0,  countCAs = 0;
	ofstream outCA("outCAnglePoreIndex.tsv");
	outCA << "Number of pores = \t"<<nPors<<endl;
	outCA << "X"<<'\t'<<"Y"<<'\t'<<"Z"<<'\t'<<"cosCA"<<'\t'<<"ContAngle"<<'\t'<<"PoreI"<<'\n';
	while(CAFil.good()) {
		++countLines;
		double xx,yy,zz, cosCA, CA;//, KcSm2, KCL, KCL; //w->wall
		CAFil >> xx>>yy>>zz>>cosCA>>CA;
		#ifdef DEVEL //EXP_0:
		double As, Aw, Kc, Kcw, KcSm;//, KcSm2, KCL, KCL; //w->wall
		CAFil >>As>>Aw>>Kc>>Kcw>>KcSm; getline(CAFil,tmpc);
		#endif  //EXP_0;

		int ii = max(2,min(int((xx-X0.x)/dx.x+1.49),nn.x-3)); // VElem has +-one extra voxels
		int jj = max(2,min(int((yy-X0.y)/dx.y+1.49),nn.y-3));
		int kk = max(2,min(int((zz-X0.z)/dx.z+1.49),nn.z-3));
		ensure(int((zz-X0.z)/dx.z+1.49) < nn.z+10,_s(ii)+" "+_s(jj)+" "+_s(kk)
			+" "+_s((zz-X0.z)/dx.z+1.49)+" dx: "+_s(dx)+" xo: "+_s(X0)+" n: "+_s(nn),-1);

		int EI = vimage(ii,jj,kk);

		if(EI<0)  EI = maxNei(vimage,ii,jj,kk,-1,1);
		if(EI<0)  EI = maxNei(vimage,ii,jj,kk,-2,2);

		outCA << xx<<'\t'<<yy<<'\t'<<zz<<'\t'<<cosCA<<'\t'<<CA<<'\t'<<EI<<'\t'<<'\n';

		if(EI>=0 && EI<nPors) {
			cosCAs[EI] += cosCA;  sumWs[EI] += 1.;
			++countCAs;
		}
	}
	cosCAs/=sumWs+1e-36;
	cout<< "n_CAs = " << countLines <<", valid:"<<int(sum(sumWs)/(countLines+1e-36)*1000)/10<< "%\n :/ \n\n" << endl;

	ofstream outPorCA("outPoresAvgCAngle.tsv");
	outPorCA<<"PoreIndx\t avgCosCA\t avgContAngle\t NContactLines\n";
	for (int EI=0; EI<nPors; ++EI) {
		double CA = sumWs[EI]>0.01 ? acos(cosCAs[EI])*180/3.1415926535 : -1;
		outPorCA<<EI <<'\t'<<cosCAs[EI]<<'\t'<<CA<<'\t'<<sumWs[EI]<<'\n';
	}

	cout<< "\n ***** END ***** \n\n\n" << endl;

	return 0;
}
