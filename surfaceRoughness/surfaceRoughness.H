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
    Rock surface roughness calculation

\*---------------------------------------------------------------------------*/

#include "calcNormalVectors.H"

      scalarField pAreasS(pAreas);
      vectorField pNwS(pAw);
      forAll(pointPoints, pI)
      {
         const labelLoop& neiPoints = pointPoints[pI];
         { ///. calculating smooth weight
             scalar avgAKcp(0.);
             vector avgpNorm(0.,0.,0.);
             scalar sumWeightsKc(1e-18);
             for(label pJ: neiPoints) {
                  scalar weight=1.;
                  if( pMarks[pJ]==4 ) weight *= 0.1;
                  avgAKcp += weight * pAreas[pJ];
                  avgpNorm += pAw[pJ];
                  sumWeightsKc += 1.;
             }
             if (sumWeightsKc>1e-15)
             {
                 pAreasS[pI] = 0.5*avgAKcp/sumWeightsKc + (1.-0.5)*pAreasS[pI];
                 pNwS[pI]    =  avgpNorm/sumWeightsKc;
             }
         }
     }

     pAreasS=0.7*pAreasS+0.3*average(pAreas);
     pAw+=0.0001*pNwS;
     pNwS/=(mag(pNwS)+1e-18);
     vectorField  pNw=pAw/(mag(pAw)+1e-18);

     vectorField previousPoints(newPoints);


         ///. point-normal curvature
		 scalarField Kcw(pAreasS.size());
         scalarField pKcw(pNw.size(),0.);
         forAll(faces,fI)
         {
         label corLabel=fMarks[fI];
             if ( corLabel==1 || corLabel==3 || corLabel==4)
             {///. compute Kcw
                 const face & f=faces[fI];
                 forAll(f, pI)
                 {

                     edge e=f.faceEdge(pI);
					 vector L_E=newPoints[e[1]] - newPoints[e[0]];
                     vector Nex2=(pNs[e[0]] + pNs[e[1]]);
                     vector Ce=0.5*(newPoints[e[0]] + newPoints[e[1]]);

                     vector Le_se = 0.5*(0.5*(Nex2/*+fNorms[fI]*/) ^ (Ce-Cf[fI])) + 0.5*mag(Ce-Cf[fI])/(mag(L_E)+1e-18)*(L_E);

					 pKcw[e[0]] += Le_se & pNw[e[0]];
					 pKcw[e[1]] -= Le_se & pNw[e[1]];
                 }
             }
         }
         pKcw /= (pAreas);
         Kcw = mag(pKcw);


      Info <<"     Roughness original,"; cout.flush();
	  scalarField Ra(Kcw);
	  for (int iKern=0; iKern<kernelRadius; ++iKern)
      {
        scalarField pKcwTmp(Kcw);
        forAll(pointPoints, pI)
        {
            const labelLoop& neiPoints = pointPoints[pI];
			{
				scalar R(pKcwTmp[pI]);
				scalar sumWeights(1e-16);
				scalar sumpAreas(1e-16);
					if ( pMarks[pI] == 1 || pMarks[pI] == 3 || pMarks[pI] == 4)
					{
					  for(label pJ: neiPoints) {
						if( pMarks[pJ]!=2)
						{
							scalar weight=1.;
							R += weight * pKcwTmp[pJ]* pAreas[pJ];
							sumWeights += weight;
							sumpAreas += pAreas[pJ];
						}
					  }

					}

				if (sumWeights>0.05)
				{
					Ra[pI] = mag(R)/sumWeights;//myEdges.size();
					//Ra[pI] = mag(R)/sumpAreas;//myEdges.size();
					//Ra[pI] = (R*avgAreaR)/sumWeights;//myEdges.size();
				}

			}

		}

    }


	boundBox BB(points);
	vector BMin=BB.min();
	vector BMax=BB.max();
	scalar vxlsize=Foam::sqrt(average(pAreas));


	vector    zz(0.,0.,1.);
	vector BMid=0.5*(BB.max()+BB.min());
	scalar radius=0.25*(BB.max()[0]-BB.min()[0]  +  BB.max()[1]-BB.min()[1]);

	Info<< "vxlsize = "<<vxlsize<<endl;;
	Info<< "BMin = "<<BMin<<endl;;
	Info<< "BMax = "<<BMax<<endl;;
	Info<< "BMid = "<<BMid<<endl;;
	(Info<< " *  ").flush();


		///Ra_x.txt
		std::ofstream  offRaX("Ra_x.txt");
		offRaX<<"Ra 1 "<<pointPoints.size()<<" float"<<std::endl;

		forAll(pointPoints, pI)
		{
			label currentLabel = pMarks[pI];
			if (currentLabel==1 || currentLabel==3 || currentLabel==4)
			{
				if(
				///cylindrical image
					mag(points[pI][2]-BMin[2])>5.*vxlsize
				 && mag(points[pI][2]-BMax[2])>5.*vxlsize
				 &&	mag(BMid - points[pI] - ((BMid - points[pI])&zz)*zz)<(radius-5.*vxlsize)

				///cubical image
					//mag(points[pI][0]-BMin[0])>5.*vxlsize
				 //&& mag(points[pI][0]-BMax[0])>5.*vxlsize
				 //&& mag(points[pI][1]-BMin[1])>5.*vxlsize
				 //&& mag(points[pI][1]-BMax[1])>5.*vxlsize
				 //&& mag(points[pI][2]-BMin[2])>5.*vxlsize
				 //&& mag(points[pI][2]-BMax[2])>5.*vxlsize
				  )
				{
					offRaX<<"";
					offRaX<<points[pI][0]<<" "<<points[pI][1]<<" "<<points[pI][2]<<" "<<Ra[pI]<<" \n";
				}
			}
		}
		offRaX.close();
		Info  <<"Ra_x.txt is written "; cout.flush();
		///pKcw_x.txt
		std::ofstream  ofpKcw("pKcw_x.txt");
		ofpKcw<<"ppKcw 1 "<<pointPoints.size()<<" float"<<std::endl;
        forAll(pointPoints, pI)
        {
			label currentLabel = pMarks[pI];
			if (currentLabel==1||currentLabel==3|| currentLabel==4)
			{
				if(
				///cylindrical image
					mag(points[pI][2]-BMin[2])>5.*vxlsize
				 && mag(points[pI][2]-BMax[2])>5.*vxlsize
				 &&	mag(BMid - points[pI] - ((BMid - points[pI])&zz)*zz)<(radius-5.*vxlsize)

				 ///cubical image
					//mag(points[pI][0]-BMin[0])>5.*vxlsize
				 //&& mag(points[pI][0]-BMax[0])>5.*vxlsize
				 //&& mag(points[pI][1]-BMin[1])>5.*vxlsize
				 //&& mag(points[pI][1]-BMax[1])>5.*vxlsize
				 //&& mag(points[pI][2]-BMin[2])>5.*vxlsize
				 //&& mag(points[pI][2]-BMax[2])>5.*vxlsize
				  )
				  {
					ofpKcw<<"";
					ofpKcw<<points[pI][0]<<" "<<points[pI][1]<<" "<<points[pI][2]<<" "<<pKcw[pI]<<" \n";
				  }
			}
       }
       ofpKcw.close();
       Info  <<"pKcw.txt is written "; cout.flush();

       ///pA_x.txt
		std::ofstream  ofpA("pA_x.txt");
		ofpA<<"pArea 1 "<<pointPoints.size()<<" float"<<std::endl;
        forAll(pointPoints, pI)
        {
			label currentLabel = pMarks[pI];
			if (currentLabel==1||currentLabel==3)
			{
				if(
				///cylindrical image
					mag(points[pI][2]-BMin[2])>5.*vxlsize
				 && mag(points[pI][2]-BMax[2])>5.*vxlsize
				 &&	mag(BMid - points[pI] - ((BMid - points[pI])&zz)*zz)<(radius-5.*vxlsize)

				 ///cubical image
					//mag(points[pI][0]-BMin[0])>5.*vxlsize
				 //&& mag(points[pI][0]-BMax[0])>5.*vxlsize
				 //&& mag(points[pI][1]-BMin[1])>5.*vxlsize
				 //&& mag(points[pI][1]-BMax[1])>5.*vxlsize
				 //&& mag(points[pI][2]-BMin[2])>5.*vxlsize
				 //&& mag(points[pI][2]-BMax[2])>5.*vxlsize
				  )
				  {
					ofpA<<"";
					ofpA<<points[pI][0]<<" "<<points[pI][1]<<" "<<points[pI][2]<<" "<<pAreas[pI]<<" \n";
				  }
			}
       }
       ofpA.close();
       Info  <<"pA_x.txt is written "; cout.flush();


		///Ra.txt
		std::ofstream  offRa("Ra.txt");
		offRa<<"POINT_DATA "<<pointPoints.size()<<std::endl;
		offRa<<"FIELD attributes 2 "<<std::endl;
		offRa<<"Ra 1 "<<pointPoints.size()<<" float"<<std::endl;
		forAll(pMarks, pI)
		{
			label currentLabel = pMarks[pI];
			if (currentLabel==1 || currentLabel==3 )
			{
				offRa<< Ra[pI] <<" \n";
			}else{
				offRa<<  0 <<"\n";
			}
		}
		offRa<<"pKcw 1 "<<pointPoints.size()<<" float"<<std::endl;
		forAll(pMarks, pI)
		{
			label currentLabel = pMarks[pI];
			if (currentLabel==1 || currentLabel==3 )
			{
				offRa<< pKcw[pI] <<" \n";
			}else{
				offRa<<  0 <<"\n";
			}
		}
		offRa.close();
		Info  <<"Ra.txt is written "; cout.flush();

		Info<< "  maxRa:  "<<max(Ra); cout.flush();
		Info<< "  avgRa:  "<<average(Ra); cout.flush();
		Info<<endl;
