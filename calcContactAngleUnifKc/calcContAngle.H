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
    Contact angle calculation

\*---------------------------------------------------------------------------*/

#include <fstream>
 #include <iostream>

    Info<< " calculating contact angles" << endl;


#include "calcNormalVectors.H"


	///. point-normal vectors correction at contact line (BC), using linear extrapolation TODO: Change 3.333, remove ...
	forAll(pNeips, vertI)
	{
		const labelLoop& neiPoints = pNeips[vertI];
		if( pMarks[vertI] == 4  ) ///. only for contact line
		{
			vector  pNsNei(0,0,0);
			forAll(neiPoints, neiPointI)
			{
				const label neiVertI = neiPoints[neiPointI];
				if( pMarks[neiVertI] == 2  )
					pNsNei += pNs[neiVertI];
			}
			if (mag(pNsNei)>0.1)
			{
				pNsNei /= mag(pNsNei);
				pNs[vertI]=1.333*pNs[vertI]-0.333*pNsNei;
				pNs[vertI] /= mag(pNs[vertI]);
			}
		}
   }



/// /////////////////////////////////////////////


  for (int iKern=0; iKern<3; ++iKern)
  {
	const vectorField pAwTmp(pAw);

	forAll(pNeips, vertI)
	{
		vector avgDisp(vector::zero);
		const labelLoop& neiPoints = pNeips[vertI];
		if(pMarks[vertI]!=2)
		forAll(neiPoints, neiPointI)
		{
			const label neiVertI = neiPoints[neiPointI];
			if( pMarks[neiVertI]!=2 )
				avgDisp += pAwTmp[neiVertI];
		}
		pAw[vertI] = avgDisp;//myEdges.size();
	  }
	}
  vectorField  pNw=pAw/(mag(pAw)+1e-12);



  for (int iKern=0; iKern<5; ++iKern)
  {
	const vectorField pNsTmp(pNs);
	forAll(pNeips, vertI)
	{
		vector avgDisp(vector::zero);
		const labelLoop& neiPoints = pNeips[vertI];
		if(pMarks[vertI]==4)
		forAll(neiPoints, neiPointI)
		{
			const label neiVertI = neiPoints[neiPointI];
			if( pMarks[neiVertI]==2 ) 	avgDisp += 0.5*pNsTmp[neiVertI];
			if( pMarks[neiVertI]==4 )	avgDisp += pNsTmp[neiVertI];
		}
		pNs[vertI] = avgDisp;//myEdges.size();
	  }
	}
   pNs=pNs/(mag(pNs)+1e-12);







//~ exit(0);
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------

	scalarField cosContactAngle(pNeips.size(),2.);


	int countCA = 0;
	forAll(pNeips, vertI)
	{

		if(pMarks[vertI] == 4)
		{
			cosContactAngle[vertI] = pNs[vertI] & pNw[vertI] ;
			++countCA;
	    }

    }


 	boundBox BB(points);
	vector BMin=BB.min();
	vector BMax=BB.max();
	scalar vxlsize=Foam::sqrt(average(pAreas));


	vector    zz(0.,0.,1.);
	vector BMid=0.5*(BB.max()+BB.min());
	scalar radius=0.25*(BB.max()[0]-BB.min()[0]  +  BB.max()[1]-BB.min()[1]);

	Info<< "vxlsize = "<<vxlsize<<endl;
	Info<< "image min    = "<<BMin<<endl;
	Info<< "image max    = "<<BMax<<endl;
	Info<< "image mid    = "<<BMid<<endl;
	Info<< "image radius = "<<radius<< " NOT USED! "<<endl;
	(Info<< " *  ").flush();

		std::ofstream  offContAngle("contactAngles.txt");
		std::ofstream  offContAngleX("contactAngles_x.txt");

		offContAngle<<"\n\n CA 1 "<<pNeips.size()<<" float"<<std::endl;
		offContAngleX<<"CA 1 "<<pNeips.size()<<" "<<countCA<<" float"<<std::endl;

	double CAsum = 0., CAmin = 1e16, CAmax = -1e16;
	int    CLnum = 1e-36;
	constexpr double rad2Degree = 180/3.1415926535;


	for (int i=0; i<pMarks.size(); ++i)
	{
		if(pMarks[i] == 4) {
			double cosCA = -cosContactAngle[i];
			offContAngleX<<" "<<points[i][0]<<" "<<points[i][1]<<" "<<points[i][2] << " "<<cosCA<<" "<<std::acos(cosCA)*rad2Degree <<"\n";

			CAmin = std::min(CAmin, cosCA);
			CAmax = std::max(CAmax, cosCA);
			CAsum += cosCA;
			++CLnum;
		}
	}
	// Used in regression tests:
	std::cout<<"\n  ContactAngMax: "<< rad2Degree*std::acos(CAmin)
	         <<"\n  ContactAngMin: "<< rad2Degree*std::acos(CAmax)
	         <<"\n  ContactAngAvg: "<< rad2Degree*std::acos(CAsum/CLnum) <<std::endl;


        forAll(pMarks, vertI)
        {
           if(pMarks[vertI] == 4)
            {
				offContAngle<<180-std::acos(cosContactAngle[vertI])*180/3.1415926535 <<" "<<"\n";
				//~ offContAngleX<<" "<<points[i][0]<<" "<<points[i][1]<<" "<<points[i][2] << "   "<<-(cosContactAngle[i])<<"  "<<std::acos(cosContactAngle[i])*180/3.1415926535 <<"\n";
			}
			if(pMarks[vertI] != 4)
			{
				offContAngle<<  0 <<" "<<"\n";
			}
       }




		offContAngle.close();
		offContAngleX.close();
