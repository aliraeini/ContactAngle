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

#include "calcNormalVectors.H"

      scalarField pAreasS(pAreas);
      vectorField pNwS(pAw);
      forAll(pNeips, pI)
      {
         const labelLoop& neiPoints = pNeips[pI];

         { ///. calculating smooth weight
             scalar avgAKcp(0.);
             vector avgpNorm(0.,0.,0.);
             scalar sumWeightsKc(1e-18);
             forAll(neiPoints, neiPointI)
             {
                 const label pJ = neiPoints[neiPointI];
                 {    scalar weight=1.;
                     if( pMarks[pJ]==4 ) weight *= 0.1;
                     avgAKcp += weight * pAreas[pJ];
                     avgpNorm += pAw[pJ];
                     sumWeightsKc += 1.;
                 }
             }
             if (sumWeightsKc>1e-15)
             {
                 pAreasS[pI]   =    0.5*avgAKcp/sumWeightsKc + (1.-0.5)*pAreasS[pI];
                 pNwS[pI]   =    avgpNorm/sumWeightsKc;
             }
         }
     }

     pAreasS=0.7*pAreasS+0.3*average(pAreas);
     pAw+=0.0001*pNwS;
     pNwS/=(mag(pNwS)+1e-18);
     vectorField  pNw=pAw/(mag(pAw)+1e-18);

     vectorField previousPoints(newPoints);



         Info <<iter<<"     Curvature,"; cout.flush();

         ///. point-normal curvature
         scalarField pKc(pNs.size(),0.);
         scalarField pKcw(pNw.size(),0.);
         scalarField pKcCL(pNs.size(),0.);
         scalarField pLCL(pNs.size(),1e-12);

         forAll(faces,fI)
         {/// compute Kc
         label corLabel=fMarks[fI];
             if( corLabel == 2 )
             {
                 const face & f=faces[fI];
                 forAll(f, pI)
                 {/// compute Kc

                     edge e=f.faceEdge(pI);
					 vector Lsx2=newPoints[e[1]] - newPoints[e[0]];
                     vector Nex2=(pNs[e[0]] + pNs[e[1]]);
					 vector Le=0.5*(newPoints[e[0]] + newPoints[e[1]])-Cf[fI];

                     /// big impact...!!
                     //vector Le_se = 0.5*(0.5*Nex2+fNorms[fI]) ^ Le; ///. wrong
                     //vector Le_se = 0.33333*(Nex2+fNorms[fI]) ^ Le;
                     vector Le_se = 0.5*(0.5*(Nex2/*+fNorms[fI]*/) ^ Le) + 0.5*mag(Le)/(mag(Lsx2)+1e-18)*(Lsx2);



					 pKc[e[0]] += Le_se & pNs[e[0]];
					 pKc[e[1]] -= Le_se & pNs[e[1]];
					 if (pMarks[e[0]]==4 && pMarks[e[1]]==4)
					 {
						Le_se = (Lsx2)/(mag(Lsx2)+1e-18);
						pKcCL[e[0]] += Le_se & pNs[e[0]];
						pKcCL[e[1]] -= Le_se & pNs[e[1]];
						pLCL[e[0]] += (mag(Lsx2));
					 }
                 }
             //}else if ( corLabel!=4 && corLabel!=2)
             }else if ( corLabel==1 || corLabel==3)
             {///. compute Kcw
                 const face & f=faces[fI];
                 forAll(f, pI)
                 {

                     edge e=f.faceEdge(pI);
					 vector Lsx2=newPoints[e[1]] - newPoints[e[0]];
                     vector Nex2=(pNs[e[0]] + pNs[e[1]]);
                     vector Le=0.5*(newPoints[e[0]] + newPoints[e[1]])-Cf[fI];

                     //vector Le_se = 0.5*(0.5*Nex2+fNorms[fI]) ^ Le; ///. wrong
                     //vector Le_se = 0.33333*(Nex2+fNorms[fI]) ^ Le;
                     vector Le_se = 0.5*(0.5*((Nex2/*+fNorms[fI]*/) ^ Le) + (mag(Le)/(mag(Lsx2)+1e-18))*Lsx2);

					 pKcw[e[0]] += Le_se & pNw[e[0]];
					 pKcw[e[1]] -= Le_se & pNw[e[1]];
                 }
             }
         }
         pKc /= (pAreas);
         pKcw /= (pAreas);
         pKcCL /= pLCL;
		pAreas=pAreastmp;












	  scalarField pKcCLSmooth(pKcCL);
      for (int iKern=0; iKern<kernelRadius; ++iKern)
      {
			scalarField pKcCLSmoothTmp(pKcCLSmooth);
			forAll(pNeips, pI)
			{
				const labelLoop& neiPoints = pNeips[pI];
				label markI = pMarks[pI];
				if( markI == 4 )
				{ ///. calculating smooth curvature
					scalar avgAKcp(pKcCLSmoothTmp[pI]);
					scalar sumWeightsKc(1.);
					forAll(neiPoints, neiPointI)
					{
						const label pJ = neiPoints[neiPointI];
						if( pMarks[pJ] == 4  )
						{
							scalar weight=1.;
							avgAKcp += weight * pKcCLSmoothTmp[pJ];
							sumWeightsKc += weight;
						}
					}
					if (sumWeightsKc>0.5)
						pKcCLSmooth[pI] =     0.5*avgAKcp/sumWeightsKc + (1.-0.5)*pKcCLSmoothTmp[pI];
					else
						pKcCLSmooth[pI] = pKcCLSmoothTmp[pI];
				}
		   }
	  }


	  scalarField pKcSmooth(pKc);
      for (int iKern=0; iKern<kernelRadius; ++iKern)
      {
			scalarField pKcSmoothTmp(pKcSmooth);
			forAll(pNeips, pI)
			{

				const labelLoop& neiPoints = pNeips[pI];
				label markI = pMarks[pI];
				if( markI == 2 || markI == 4  )
				{ ///. calculating smooth curvature
					scalar avgAKcp(pKcSmoothTmp[pI]);
					scalar sumWeightsKc(1.);
					forAll(neiPoints, neiPointI)
					{
						const label pJ = neiPoints[neiPointI];
						if( pMarks[pJ] == 2 || pMarks[pJ] == 4  )
						{
							scalar weight=pWeights[pJ];  weight*=weight;weight*=weight;
							if( pMarks[pJ] == 4) weight*=0.01;
							//if( pMarks[pJ] == 4 && markI == 4) weight=10.;
							avgAKcp += weight * pKcSmoothTmp[pJ];
							sumWeightsKc += weight;
						}
					}
					if (sumWeightsKc>0.5)
						pKcSmooth[pI] =     0.5*avgAKcp/sumWeightsKc + (1.-0.5)*pKcSmoothTmp[pI];
					else
						pKcSmooth[pI] = pKcSmoothTmp[pI];
				}
		   }
	  }



	  scalarField pKcSmooth2(pKcSmooth);
      for (int iKern=0; iKern<kernelRadius+2; ++iKern)
      {
			scalarField pKcSmooth2Tmp(pKcSmooth2);
			forAll(pNeips, pI)
			{
				const labelLoop& neiPoints = pNeips[pI];
				label markI = pMarks[pI];
				if( markI == 2 || markI == 4  )
				{ ///. calculating smooth curvature
					scalar avgAKcp(pKcSmooth2Tmp[pI]);
					scalar sumWeightsKc(1.);
					forAll(neiPoints, neiPointI)
					{
						const label pJ = neiPoints[neiPointI];
						if( pMarks[pJ] == 2 || pMarks[pJ] == 4  )
						{
							scalar weight=pWeights[pJ];  weight*=weight;weight*=weight;
							if( pMarks[pJ] == 4) weight*=0.01;
							avgAKcp += weight * pKcSmooth2Tmp[pJ];
							sumWeightsKc += weight;
						}
					}
					if (sumWeightsKc>0.5)
						pKcSmooth2[pI] =     0.5*avgAKcp/sumWeightsKc + (1.-0.5)*pKcSmooth2Tmp[pI];
					else
						pKcSmooth2[pI] = pKcSmooth2Tmp[pI];
				}

		   }
	  }







        vectorField displacementsRc(newPoints.size(), vector(0,0,0));
        forAll(pNeips, pI)
        {


            label markI = pMarks[pI];

             if( markI == 4  )
             { ///. contact line curvature uniformization
                scalar DispCL = 1.*relaxCL*( -1.*(pKcSmooth[pI]-pKcSmooth2[pI]) + 0.9*(pKc[pI]-pKcSmooth[pI])+ 1.*(pKcCL[pI]-pKcCLSmooth[pI]))
                *Foam::sqrt(pAreas[pI]/ (pKcSmooth[pI]*pKcSmooth[pI] +0.1/pAreasS[pI]));
                newPoints[pI] +=  DispCL*(pNs[pI] - (pNs[pI] & pNw[pI])*pNw[pI]);
             }
             else  if(markI == 2)
             {  ///. curvature uniformization
                 newPoints[pI] +=  pNs[pI] * relax*(0.5*(pKc[pI]-pKcSmooth[pI])+8.*(pKcSmooth[pI]-pKcSmooth2[pI]))
												*Foam::sqrt(pAreas[pI]/ (pKcSmooth[pI]*pKcSmooth[pI]  +0.1/pAreasS[pI]));
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








		///K_x.txt
		std::ofstream  offKcX("Kc_x.txt");
		offKcX<<"K 1 "<<pNeips.size()<<" float"<<std::endl;
        forAll(pNeips, pI)
        {
			label markI = pMarks[pI];
            if(iter==nIters-1)
            {
				// if (markI==2||markI==4)
				if (markI==2 && markI!=4)
				{
					offKcX<<"";
					offKcX<<points[pI][0]<<" "<<points[pI][1]<<" "<<points[pI][2]<<" "<<pKc[pI]<<" \n";
				}
			}
       }
       offKcX.close();

       	///pA_x.txt
		std::ofstream  ofpA("pA_x.txt");
		ofpA<<"pArea 1 "<<pNeips.size()<<" float"<<std::endl;
        forAll(pNeips, pI)
        {
			label markI = pMarks[pI];
            if(iter==nIters-1)
            {
				if (markI==1||markI==3 ||markI==4)
				{
					ofpA<<"";
					ofpA<<points[pI][0]<<" "<<points[pI][1]<<" "<<points[pI][2]<<" "<<pAreas[pI]<<" \n";
				}
			}
       }
       ofpA.close();

		std::ofstream  offRc("Kc.txt");
		offRc<<"POINT_DATA "<<pNeips.size()<<std::endl;
		offRc<<"FIELD attributes 6 "<<std::endl;
		offRc<<"Kc 1 "<<pNeips.size()<<" float"<<std::endl;
        forAll(pNeips, pI)
        {
            if(iter==nIters-1)
            {
				if ((pI+1)%1==0)	offRc<<"";
				offRc<<pKc[pI]<<" "<<"  \n";
			}
       }
       offRc<<"Kcw 1 "<<pNeips.size()<<" float"<<std::endl;
        forAll(pNeips, pI)
        {
            if(iter==nIters-1)
            {
				if ((pI+1)%1==0)	offRc<<"";
				offRc<<pKcw[pI]<<" "<<"  \n";
			}
       }
		offRc<<"\n\n RcCL 1 "<<pNeips.size()<<" float"<<std::endl;
        forAll(pNeips, pI)
        {
            if(iter==nIters-1)
            {
				if ((pI+1)%1==0)	offRc<<"";
				offRc<<pKcCL[pI]<<" "<<"  \n";
			}
       }
		offRc<<"\n\n InterfNorms 3 "<<pNeips.size()<<" float"<<std::endl;
        forAll(pNeips, pI)
        {
            if(iter==nIters-1)
            {
				if ((pI+1)%1==0)	offRc<<"";
				offRc<<pNs[pI][0]<<" "<<pNs[pI][1]<<" "<<pNs[pI][2]<<" "<<" \n";
			}
       }
		offRc<<"\n\n wallNorms 3 "<<pNeips.size()<<" float"<<std::endl;
        forAll(pNeips, pI)
        {
            if(iter==nIters-1)
            {
				if ((pI+1)%1==0)	offRc<<"";
				offRc<<pNw[pI][0]<<" "<<pNw[pI][1]<<" "<<pNw[pI][2]<<" "<<"  \n";
			}
       }

		offRc.close();

		//newPoints += displacementsRc;



	//===========================================================================
	//===========================================================================
	//================= Gausian like smoothing during curvature move ===========
	//================== to smooth point positions laterally ====================
	//===========================================================================


     vectorField previousPoints2(newPoints);



	vectorField displacements(newPoints.size(), vector(0,0,0));
	forAll(pNeips, pI)
	{
		const labelLoop& neiPoints = pNeips[pI];
		label markI = pMarks[pI];

		vector avgPos(vector::zero);
		scalar sumWeights(1e-28);
		if ( pMarks[pI] == 4 )
		{
			forAll(neiPoints, neiPointI)
			{
				const label pJ = neiPoints[neiPointI];

				scalar weight=(pAreas[pJ]+0.1*pAreas[pI]);
				weight*=weight*pWeights[pJ];
				if( pMarks[pJ]!=markI ) weight*=0.1;  ///.  TODO change 0.5,   , affects convergence
				point neiPos=newPoints[pJ];
				if( pMarks[pJ]==2 )
				{
					neiPos -= 1.01*((neiPos-newPoints[pI])&pNw[pI])*pNw[pI];    ///.  TODO change 1.05
				}
				else if( pMarks[pJ] == 4) weight *= 8.; ///. TODO: Change 3.333, remove ...
				avgPos += weight * neiPos;
				sumWeights += weight;
			}

		}
		else
		{
			forAll(neiPoints, neiPointI)
			{
				const label pJ = neiPoints[neiPointI];

				scalar weight=(pAreas[pJ]+0.1*pAreas[pI]);//*max(NormalsOld[pJ]&NormalsOld[pI],0.001);// /mag(newPoints[pJ]-newPoints[pI]);
				weight*=weight*pWeights[pJ];
				//if( pMarks[pJ]!=markI ) weight*=0.9;     ///.  TODO change 0.7, smaller value makes CL faces bigger
				//if( pMarks[pJ]==3 ) weight*=2.;     ///.  TODO change 0.7, smaller value makes CL faces bigger
				avgPos += weight * newPoints[pJ];
				sumWeights += weight;

			}
		}

		if (sumWeights>1e-18)
		{
			avgPos /= sumWeights;//myEdges.size();
			displacements[pI] = avgPos-newPoints[pI]; //(0.8*relax)*((avgPos-newPoints[pI])&pNw[pI])*pNw[pI] + (0.2*relax)*(avgPos-newPoints[pI]);
		}

      }
	  newPoints += 0.5*relax*(displacements - 0.5*(displacements& pNw)*pNw);
	  //newPoints += 0.03*relaxCL*displacementsCL;


/// moving contact line back to the solid wall (alleviate the crunch effect)

	displacements=newPoints-previousPoints2;

	for (int iKern=0; iKern<2; ++iKern)
	{
 	  vectorField displacementsTmp = displacements;
      forAll(pNeips, pI)
      {
		label markI = pMarks[pI];
		bool markIEq2 = markI==2;
		//if( markI != 2  )
		{ ///. contact line smoothing Vol preserve
			const labelLoop& neiPoints = pNeips[pI];
			vector avgDisp(vector::zero);
			scalar sumWeights(1e-28);
			forAll(neiPoints, iii)
			{
				const label pJ = neiPoints[iii];
				//if( pMarks[pJ] == markI  )
				if( pMarks[pJ]==markI || !( (pMarks[pJ]==2) ^ markIEq2)  || pMarks[pJ] == 4  )
				{
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
	//Info<< " "<< 0.004+(double(average(max(min((pKc-0.004)*1e8,1e6),-1e6)*max(pMarks-3,0))))/(0.0001+double(average(100000*max(pMarks-3,0))))/1000.; cout.flush();
	Info<< " "<<max(pKc); cout.flush();
	Info<< "],  DispMag  max:  "<<max(mag(displacements)); cout.flush();
	Info<< "  avg:  "<<average(mag(displacements)); cout.flush();
	Info<<endl;
