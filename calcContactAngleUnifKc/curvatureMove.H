/*-------------------------------------------------------------------------*\
You can redistribute this code and/or modify this code under the
terms of the GNU General Public License (GPL) as published by the
Free Software Foundation, either version 3 of the License, or (at
your option) any later version. see <http://www.gnu.org/licenses/>.


The code has been developed by Ahmed AlRatrout as a part his PhD
at Imperial College London, under the supervision of Dr. Branko Bijeljic
and Prof. Martin Blunt.

Please see our website for relavant literature:
AlRatrout et al, AWR, 2017 - https://www.sciencedirect.com/science/article/pii/S0309170817303342
AlRatrout et al, WRR, 2018 - https://agupubs.onlinelibrary.wiley.com/doi/pdf/10.1029/2017WR022124
AlRatrout et al, PNAS, 2018 - http://www.pnas.org/

For further information please contact us by email:
Ahmed AlRarout:  a.alratrout14@imperial.ac.uk
Branko Bijeljic: b.bijeljic@imperial.ac.uk
Martin J Blunt:  m.blunt@imperial.ac.uk

Description
    Volume-preserving curvature unifrom smoothing

\*---------------------------------------------------------------------------*/

/// Numenclature fI: face index,  pI:point index, pJ: adjacent-point index

#include "calcCurvatureNormals.H"



vectorField previousPoints(newPoints);

for_(pNeips, pI) {  //curvature uniformization

	label markI = pMarks[pI];

    if( markI == 4  ) { ///. contact line curvature uniformization
                scalar DispCL = 1.*relaxCL*( -1.*(pKcSmooth[pI]-pKcSmooth2[pI]) + 0.9*(pKc[pI]-pKcSmooth[pI])+ 1.*(pKcCL[pI]-pKcCLSmooth[pI]))
                *sqrt(pAreas[pI]/ (pKcSmooth[pI]*pKcSmooth[pI] +0.1/pAreasS[pI]));
                newPoints[pI] +=  DispCL*(pNs[pI] - (pNs[pI] & pNw[pI])*pNw[pI]);
            }
            else  if(markI == 2) {  ///. curvature uniformization
                scalar DispCL = relax*(0.5*(pKc[pI]-pKcSmooth[pI])+8.*(pKcSmooth[pI]-pKcSmooth2[pI]))
                            *Foam::sqrt(pAreas[pI]/ (pKcSmooth[pI]*pKcSmooth[pI]  +0.1/pAreasS[pI]));
                newPoints[pI] +=  pNs[pI] * DispCL;
            }



	//if( markI == 4  )
	//{ ///. contact line curvature uniformization
		//vector DispCL =  relaxCL*(-0.05*(pKcSmooth[pI]-pKc[pI])+ 5.*(pKcCL[pI]-pKcCLSmooth[pI]))* pAreas[pI]*pNs[pI];
		//displacementsRc[pI] =  (DispCL - (DispCL & pNw[pI])*pNw[pI]);
	//}
	//else  if(markI == 2)
	//{  ///. curvature uniformization
		 //displacementsRc[pI] =  relax*(pKc[pI]-pKcSmooth[pI]) * pAreas[pI]*pNs[pI];
	//}
}


{
	#include "./write.H"
	//#include "./calcCurvature.H"
}
//===========================================================================
//===========================================================================
//================= Gausian like smoothing during curvature move ===========
//================== to smooth point positions laterally ====================
//===========================================================================


vectorField previousPoints2(newPoints);



vectorField displacements(newPoints.size(), dbl3(0,0,0));
for_(pNeips, pI)  { // Gaussian
	const labelLoop& neiPoints = pNeips[pI];
	label markI = pMarks[pI];

	dbl3 avgPos(dbl3::zero);
	scalar sumWeights(1e-28);
	if ( pMarks[pI] == 4 ) {
		for(int pJ: neiPoints) {

			scalar weight=(pAreas[pJ]+0.1*pAreas[pI]);
			weight*=weight*pWeights[pJ];
			if( pMarks[pJ]!=markI ) weight*=0.1;  ///.  TODO change 0.5,   , affects convergence
			point neiPos=newPoints[pJ];
			if( pMarks[pJ]==2 ) {
				neiPos -= 1.01*((neiPos-newPoints[pI])&pNw[pI])*pNw[pI];    ///.  TODO change 1.05
			}
			else if( pMarks[pJ] == 4) weight *= 8.; ///. TODO: Change 3.333, remove ...
			avgPos += weight * neiPos;
			sumWeights += weight;
		}

	}
	else
	{
		for(int pJ: neiPoints) {

			scalar weight=(pAreas[pJ]+0.1*pAreas[pI]);//*max(NormalsOld[pJ]&NormalsOld[pI],0.001);// /mag(newPoints[pJ]-newPoints[pI]);
			weight*=weight*pWeights[pJ];
			//if( pMarks[pJ]!=markI ) weight*=0.9;     ///.  TODO change 0.7, smaller value makes CL faces bigger
			//if( pMarks[pJ]==3 ) weight*=2.;     ///.  TODO change 0.7, smaller value makes CL faces bigger
			avgPos += weight * newPoints[pJ];
			sumWeights += weight;

		}
	}

	if (sumWeights>1e-18)	{
		avgPos /= sumWeights;//myEdges.size();
		displacements[pI] = avgPos-newPoints[pI]; //(0.8*relax)*((avgPos-newPoints[pI])&pNw[pI])*pNw[pI] + (0.2*relax)*(avgPos-newPoints[pI]);
	}

}
newPoints += 0.5*relax*(displacements - 0.5*(displacements& pNw)*pNw);
//newPoints += 0.03*relaxCL*displacementsCL;



/// moving contact line back to the solid wall (alleviate the crunch effect)

displacements=newPoints-previousPoints2;

for (int iKern=0; iKern<2; ++iKern)  {  // volume preserving
  vectorField displacementsTmp = displacements;
	for_(pNeips, pI)  {
		label markI = pMarks[pI];
		bool markIEq2 = markI==2;
		//if( markI != 2  )
		{ ///. contact line smoothing Vol preserve
			const labelLoop& neiPoints = pNeips[pI];
			dbl3 avgDisp(0.,0.,0.);
			scalar sumWeights(1e-28);
			for(int pJ: neiPoints) {
				//if( pMarks[pJ] == markI  )
				if( pMarks[pJ]==markI || !( (pMarks[pJ]==2) ^ markIEq2)  || pMarks[pJ] == 4  ) {
					scalar weight=pAreas[pJ];
					avgDisp += weight * displacementsTmp[pJ];
					sumWeights += weight;
				}
			}

			displacements[pI] = avgDisp/sumWeights;

		}
	}
}

displacements = 0.3*displacements+0.7*(displacements& pNw)*pNw;
displacements = 0.3*displacements+0.7*(displacements& pNs)*pNs;
newPoints -= 1.09*displacements;     ///.  TODO change 0.5 0.5,  sum should be 1

Info  <<"CL "; cout.flush();




surf123.movePoints(newPoints);



displacements = newPoints-previousPoints;
Info<< "  A:  ["<<min(pAreas)<<" "<<max(pAreas)<<"] "; cout.flush();
Info<< "  k:  ["<<min(pKc); cout.flush();
Info<< " "<<max(pKc); cout.flush();
Info<< "],  DispMag  max:  "<<max(mag(displacements)); cout.flush();
Info<< "  avg:  "<<average(mag(displacements)); cout.flush();
Info<<endl;
