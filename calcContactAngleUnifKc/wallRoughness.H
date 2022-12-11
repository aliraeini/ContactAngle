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
      forAll(pointPoints, vertI)
      {
         const labelLoop& neiPoints = pointPoints[vertI];
         //~ label currentLabel = pMarks[vertI];
         //~ if( currentLabel == 2 || currentLabel == 4  )
         { ///. calculating smooth weight
             scalar avgAKcp(0.);
             vector avgpNorm(0.,0.,0.);
             scalar sumWeightsKc(1e-18);
             forAll(neiPoints, neiPointI)
             {
                 const label neiVertI = neiPoints[neiPointI];
                 {    scalar weight=1.;
                     if( pMarks[neiVertI]==4 ) weight *= 0.1;
                     avgAKcp += weight * pAreas[neiVertI];
                     avgpNorm += pAw[neiVertI];
                     sumWeightsKc += 1.;
                 }
             }
             if (sumWeightsKc>1e-15)
             {
                 pAreasS[vertI]   =    0.5*avgAKcp/sumWeightsKc + (1.-0.5)*pAreasS[vertI];
                 pNwS[vertI]   =    avgpNorm/sumWeightsKc;
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
         forAll(faces,faceI)
         {
         label corLabel=fMarks[faceI];
             if ( corLabel==1 || corLabel==3)
             {///. compute Kcw
                 const face & f=faces[faceI];
                 forAll(f, pI)
                 {

                     edge e=f.faceEdge(pI);
					 vector L_E=newPoints[e[1]] - newPoints[e[0]];
                     vector Nex2=(pNs[e[0]] + pNs[e[1]]);
                     vector Ce=0.5*(newPoints[e[0]] + newPoints[e[1]]);

                     vector Le_se = 0.5*(0.5*(Nex2/*+fNorms[faceI]*/) ^ (Ce-Cf[faceI])) + 0.5*mag(Ce-Cf[faceI])/(mag(L_E)+1e-18)*(L_E);

					 pKcw[e[0]] += Le_se & pNw[e[0]];
					 pKcw[e[1]] -= Le_se & pNw[e[1]];
                 }
             }
         }
         pKcw /= (pAreas);
         Kcw = mag(pKcw);

         //Info << "pKcw = "<<pKcw;cout.flush();
         //Info << "Kcw = "<<Kcw;cout.flush();

//#include "curvatureMove.H"
		//Info <<"     Kcw: "<<Kcw; cout.flush();
      Info <<"     Roughness original,"; cout.flush();
	  scalarField Ra(Kcw);
	  for (int iKern=0; iKern<kernelRadius; ++iKern)
      {
        scalarField pKcwTmp(Kcw);
        forAll(pointPoints, vertI)
        {
            const labelLoop& neiPoints = pointPoints[vertI];
			{
				//vector Ra(vector::zero);
				scalar R(pKcwTmp[vertI]);
				scalar sumWeights(1e-16);
				if ( pMarks[vertI] == 1 || pMarks[vertI] == 3 || pMarks[vertI] == 4)
				{
				  forAll(neiPoints, neiPointI)
				  {
					const label neiVertI = neiPoints[neiPointI];
					if( pMarks[neiVertI]!=2 && pKcwTmp[neiVertI]>0)
					{
						scalar weight=1.;
						R += weight * pKcwTmp[neiVertI];
						sumWeights += weight;
					}
					else{
						R = 0;
					}


				  }

			    }
				if (sumWeights>0.05)
				{
					Ra[vertI] = R/sumWeights;//myEdges.size();
				}
				//if (Ra[vertI]!=0)
				//{
					//Ra[vertI] = 1/Ra[vertI];
				//}
			}

		}

    }



	  scalarField pKcwSmooth(pKcw);
      for (int iKern=0; iKern<kernelRadius; ++iKern)
      {
			scalarField pKcwSmoothTmp(pKcwSmooth);
			forAll(pointPoints, vertI)
			{

				const labelLoop& neiPoints = pointPoints[vertI];
				label currentLabel = pMarks[vertI];
				if( currentLabel == 1 || currentLabel == 3 || currentLabel == 4  )
				{ ///. calculating smooth curvature
					scalar avgAKcp(pKcwSmoothTmp[vertI]);
					scalar sumWeightsKc(1.);
					forAll(neiPoints, neiPointI)
					{
						const label neiVertI = neiPoints[neiPointI];
						if( currentLabel == 1 || currentLabel == 3 || pMarks[neiVertI] == 4  )
						{
							scalar weight=pWeights[neiVertI];
							//weight*=weight;weight*=weight;
							//if( pMarks[neiVertI] == 4) weight*=0.01;
							//~ if( pMarks[neiVertI] == 4 && currentLabel == 4) weight=10.;
							avgAKcp += weight * pKcwSmoothTmp[neiVertI];
							sumWeightsKc += weight;
						}
					}
					if (sumWeightsKc>0.5)
						pKcwSmooth[vertI] =     0.5*avgAKcp/sumWeightsKc + (1.-0.5)*pKcwSmoothTmp[vertI];
					else
						pKcwSmooth[vertI] = pKcwSmoothTmp[vertI];
				}
		   }
	  }



	  scalarField pKcwSmooth2(pKcwSmooth);
      for (int iKern=0; iKern<kernelRadius+2; ++iKern)
      {
			scalarField pKcwSmooth2Tmp(pKcwSmooth2);
			forAll(pointPoints, vertI)
			{
				const labelLoop& neiPoints = pointPoints[vertI];
				label currentLabel = pMarks[vertI];
				if( currentLabel == 1 || currentLabel == 3 || currentLabel == 4  )
				{ ///. calculating smooth curvature
					scalar avgAKcp(pKcwSmooth2Tmp[vertI]);
					scalar sumWeightsKc(1.);
					forAll(neiPoints, neiPointI)
					{
						const label neiVertI = neiPoints[neiPointI];
						if( currentLabel == 1 || currentLabel == 3 || pMarks[neiVertI] == 4  )
						{
							scalar weight=pWeights[neiVertI];
							//weight*=weight;weight*=weight;
							//if( pMarks[neiVertI] == 4) weight*=0.01;  ///.  TODO change 0.01 : 0-1
							//~ if( pMarks[neiVertI] == 4 && currentLabel == 4) weight=10.;
							avgAKcp += weight * pKcwSmooth2Tmp[neiVertI];
							sumWeightsKc += weight;
						}
					}
					if (sumWeightsKc>0.5)
						pKcwSmooth2[vertI] =     0.5*avgAKcp/sumWeightsKc + (1.-0.5)*pKcwSmooth2Tmp[vertI];
					else
						pKcwSmooth2[vertI] = pKcwSmooth2Tmp[vertI];
				}

		   }
	  }




        vectorField displacementsRc(newPoints.size(), vector(0,0,0));
        forAll(pointPoints, vertI)
        {


            label currentLabel = pMarks[vertI];

             if(currentLabel == 1 || currentLabel == 3 || currentLabel == 4)
             {  ///. curvature uniformization
                 newPoints[vertI] +=  pNw[vertI] * relax*(0.5*(pKcw[vertI]-pKcwSmooth[vertI])+8.*(pKcwSmooth[vertI]-pKcwSmooth2[vertI]))
												*Foam::sqrt(pAreas[vertI]/ (pKcwSmooth[vertI]*pKcwSmooth[vertI]  +0.1/pAreasS[vertI]));  ///.  TODO change -2. , etc
             }

		}


	//===========================================================================
	//===========================================================================
	//================= Gausian like smoothing during curvature move ===========
	//================== to smooth point positions laterally ====================
	//===========================================================================


     vectorField previousPoints2(newPoints);



	vectorField displacements(newPoints.size(), vector(0,0,0));
	//~ vectorField displacementsCL(newPoints.size(), vector(0,0,0));
	forAll(pointPoints, vertI)
	{

		const labelLoop& neiPoints = pointPoints[vertI];
		label currentLabel = pMarks[vertI];


		vector avgPos(vector::zero);
		scalar sumWeights(1e-28);
		//if ( pMarks[vertI] == 4 )
		//{
			//forAll(neiPoints, neiPointI)
			//{
				//const label neiVertI = neiPoints[neiPointI];

				//scalar weight=(pAreas[neiVertI]+0.1*pAreas[vertI]);//*max(NormalsOld[neiVertI]&NormalsOld[vertI],0.001);// /mag(newPoints[neiVertI]-newPoints[vertI]);
				//weight*=weight*pWeights[neiVertI];
				//if( pMarks[neiVertI]!=currentLabel ) weight*=0.1;  ///.  TODO change 0.5,   , affects convergence
				//point neiPos=newPoints[neiVertI];
				//if( pMarks[neiVertI]==2 )
				//{
					//neiPos -= 1.01*((neiPos-newPoints[vertI])&pNw[vertI])*pNw[vertI];    ///.  TODO change 1.05
				//}
				//else if( pMarks[neiVertI] == 4) weight *= 8.; ///. TODO: Change 3.333, remove ...
				//avgPos += weight * neiPos;
				//sumWeights += weight;
			//}

		//}
		if ( pMarks[vertI] == 1 || pMarks[vertI] == 3 || pMarks[vertI] == 4)
		{
			forAll(neiPoints, neiPointI)
			{
				const label neiVertI = neiPoints[neiPointI];

				scalar weight=(pAreas[neiVertI]+0.1*pAreas[vertI]);//*max(NormalsOld[neiVertI]&NormalsOld[vertI],0.001);// /mag(newPoints[neiVertI]-newPoints[vertI]);
				weight*=weight*pWeights[neiVertI];
				//~ if( pMarks[neiVertI]!=currentLabel ) weight*=0.9;     ///.  TODO change 0.7, smaller value makes CL faces bigger
				//~ if( pMarks[neiVertI]==3 ) weight*=2.;     ///.  TODO change 0.7, smaller value makes CL faces bigger
				avgPos += weight * newPoints[neiVertI];
				sumWeights += weight;

		  }
		}

			if (sumWeights>1e-18)
			{
				avgPos /= sumWeights;//myEdges.size();
				displacements[vertI] = avgPos-newPoints[vertI]; //(0.8*relax)*((avgPos-newPoints[vertI])&pNw[vertI])*pNw[vertI] + (0.2*relax)*(avgPos-newPoints[vertI]);
			}


			//~ if( currentLabel == 4  )
			//~ { ///. contact line smoothing
				//~ vector avgPos(vector::zero);
				//~ scalar sumWeights(1e-16);
				//~ forAll(neiPoints, iii)
				//~ {
					//~ const label neiVertI = neiPoints[iii];
					//~ if( pMarks[neiVertI] == 4  )
					//~ {
						//~ scalar weight=1.;
						//~ avgPos += weight * newPoints[neiVertI];
						//~ sumWeights += weight;
					//~ }
				//~ }
				//~ vector DispCL =  (avgPos/sumWeights-newPoints[vertI]); // relax*((avgPos/sumWeights-newPoints[vertI])&pNs[vertI])*pNs[vertI];
				//~ DispCL =  (DispCL & pNs[vertI])*pNs[vertI];
				//~ displacementsCL[vertI] =  (DispCL - (DispCL & pNw[vertI])*pNw[vertI]);
			//~ }

      }
	  newPoints += 0.5*relax*(displacements - 0.5*(displacements& pNw)*pNw);     ///.  TODO change 0.2
	  //~ newPoints += 0.03*relaxCL*displacementsCL;


/// moving contact line back to the solid wall (alleviate the crunch effect)

	  displacements=newPoints-previousPoints2;

	 scalarField D = mag(displacements);
     for (int iKern=0; iKern<2; ++iKern)
     {
 	  vectorField displacementsTmp = displacements;
      forAll(pointPoints, vertI)
      {
		label currentLabel = pMarks[vertI];
		bool curentLabelEq2 = currentLabel==2;
		{ ///. contact line smoothing Vol preserve
			const labelLoop& neiPoints = pointPoints[vertI];
			vector avgDisp(vector::zero);
			scalar sumWeights(1e-28);
			forAll(neiPoints, iii)
			{
				const label neiVertI = neiPoints[iii];
				if( pMarks[neiVertI]==currentLabel || !( (pMarks[neiVertI]==2) ^ curentLabelEq2)  || pMarks[neiVertI] == 4  )
				{
					scalar weight=pAreas[neiVertI];
					avgDisp += weight * displacementsTmp[neiVertI];
					sumWeights += weight;
				}
			}

			displacements[vertI] = avgDisp/sumWeights;//myEdges.size();

		}
      }
     }

	  //~ newPoints -= 0.3*displacements+0.79*(displacements& pNw)*pNw;     ///.  TODO change 0.5 0.5,  sum should be 1

	displacements = 0.3*displacements+0.7*(displacements& pNw)*pNw;
	//displacements = 0.3*displacements+0.7*(displacements& pNs)*pNs;
	newPoints -= 1.05*displacements;     ///.  TODO change 0.5 0.5,  sum should be 1

		Info  <<"CL "; cout.flush();


	displacements = newPoints-previousPoints;

	Info<< "  A:  ["<<min(pAreas)<<" "<<max(pAreas)<<"] "; cout.flush();
	Info<< "  Kw:  ["<<min(pKcw); cout.flush();
	//~ Info<< " "<< 0.004+(double(average(max(min((pKc-0.004)*1e8,1e6),-1e6)*max(pMarks-3,0))))/(0.0001+double(average(100000*max(pMarks-3,0))))/1000.; cout.flush();
	Info<< " "<<max(pKcw); cout.flush();
	Info<< "],  DispMag  max:  "<<max(mag(displacements)); cout.flush();
	Info<< "  avg:  "<<average(mag(displacements)); cout.flush();
	Info<<endl;



      Info <<"     Roughness smooth,"; cout.flush();
	  //scalarField RaS(Ra);
	  for (int iKern=0; iKern<kernelRadius; ++iKern)
      {
        //scalarField pKcwTmp(pKcwSmooth2);
        forAll(pointPoints, vertI)
        {
            const labelLoop& neiPoints = pointPoints[vertI];
			{
				scalar RS(D[vertI]);
				scalar sumWeights(1e-16);
				if ( pMarks[vertI] == 1 || pMarks[vertI] == 3 || pMarks[vertI] == 4)
				{
				  forAll(neiPoints, neiPointI)
				  {
					const label neiVertI = neiPoints[neiPointI];
					if( pMarks[neiVertI]!=2 && D[neiVertI]>0)
					{
						scalar weight=1.;
						RS += weight * D[neiVertI];
						sumWeights += weight;
					}
				  }

			    }
				if (sumWeights>0.05)
				{
					//RaS[vertI] = RS/sumWeights;
					Ra[vertI] = RS/sumWeights;
				}
				//if (Ra[vertI]!=0)
				//{
					//Ra[vertI] = 1/Ra[vertI];
				//}
			}

		}
    }
    //Ra = RaS-Ra;


	///Ra_x.txt
	std::ofstream  offRaX("Ra_x.txt");
	offRaX<<"Ra 1 "<<pointPoints.size()<<" float"<<std::endl;
	forAll(pointPoints, vertI)
	{
		label currentLabel = pMarks[vertI];
		if (currentLabel==1 || currentLabel==3 || pMarks[vertI] == 4)
		{
			offRaX<<"";
			offRaX<<points[vertI][0]<<" "<<points[vertI][1]<<" "<<points[vertI][2]<<" "<<Ra[vertI]<<" \n";
		}
	}
	offRaX.close();

   ///Ra.txt
	std::ofstream  offRa("Ra.txt");
	offRa<<"Ra 1 "<<pointPoints.size()<<" float"<<std::endl;
	forAll(pMarks, vertI)
	{
		label currentLabel = pMarks[vertI];
		if (currentLabel==1 || currentLabel==3 || pMarks[vertI] == 4)
		{
			offRa<< Ra[vertI] <<" \n";
		}else{
			offRa<<  0 <<"\n";
		}
   }
   offRa.close();



	Info<< "  maxRa:  "<<max(Ra); cout.flush();
	Info<< "  avgRa:  "<<average(Ra); cout.flush();
	Info<<endl;
