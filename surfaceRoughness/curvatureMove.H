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



         Info <<iter<<"     Curvature,"; cout.flush();

         ///. point-normal curvature
         scalarField pKc(pNs.size(),0.);
         scalarField pKcw(pNw.size(),0.);
         scalarField pKcCL(pNs.size(),0.);
         scalarField pLCL(pNs.size(),1e-12);
         //~ vectorField Le_ses(newPoints.size(), vector(0,0,0));
         //~ vectorField Nex2s(newPoints.size(), vector(0,0,0));
         //~ vectorField Ces(newPoints.size(), vector(0,0,0));
         //~ vectorField L_Es(newPoints.size(), vector(0,0,0));
         forAll(faces,faceI)
         {///. compute Kc
         label corLabel=fMarks[faceI];
             if( corLabel == 2 )
             {
                 const face & f=faces[faceI];
                 forAll(f, pI)
                 {

                     edge e=f.faceEdge(pI);
					 vector L_E=newPoints[e[1]] - newPoints[e[0]];
                     vector Nex2=(pNs[e[0]] + pNs[e[1]]);
                     vector Ce=0.5*(newPoints[e[0]] + newPoints[e[1]]);

                     /// big impact...!!
                     //~ vector Le_se = 0.5*(0.5*Nex2+fNorms[faceI]) ^ (Ce-Cf[faceI]); ///. wrong
                     //~ vector Le_se = 0.33333*(Nex2+fNorms[faceI]) ^ (Ce-Cf[faceI]);
                     vector Le_se = 0.5*(0.5*(Nex2/*+fNorms[faceI]*/) ^ (Ce-Cf[faceI])) + 0.5*mag(Ce-Cf[faceI])/(mag(L_E)+1e-18)*(L_E);



					 pKc[e[0]] += Le_se & pNs[e[0]];
					 pKc[e[1]] -= Le_se & pNs[e[1]];
					 if (pMarks[e[0]]==4 && pMarks[e[1]]==4)
					 {
						Le_se = (L_E)/(mag(L_E)+1e-18);
						pKcCL[e[0]] += Le_se & pNs[e[0]];
						pKcCL[e[1]] -= Le_se & pNs[e[1]];
						pLCL[e[0]] += (mag(L_E));
					 }
                 }
             //}else if ( corLabel!=4 && corLabel!=2) 
             }else if ( corLabel==1 || corLabel==3) 
             {///. compute Kcw
                 const face & f=faces[faceI];
                 forAll(f, pI)
                 {

                     edge e=f.faceEdge(pI);
					 vector L_E=newPoints[e[1]] - newPoints[e[0]];
                     vector Nex2=(pNs[e[0]] + pNs[e[1]]);
                     vector Ce=0.5*(newPoints[e[0]] + newPoints[e[1]]);

                     //~ vector Le_se = 0.5*(0.5*Nex2+fNorms[faceI]) ^ (Ce-Cf[faceI]); ///. wrong
                     //~ vector Le_se = 0.33333*(Nex2+fNorms[faceI]) ^ (Ce-Cf[faceI]);
                     vector Le_se = 0.5*(0.5*(Nex2/*+fNorms[faceI]*/) ^ (Ce-Cf[faceI])) + 0.5*mag(Ce-Cf[faceI])/(mag(L_E)+1e-18)*(L_E);

					 pKcw[e[0]] += Le_se & pNw[e[0]];
					 pKcw[e[1]] -= Le_se & pNw[e[1]];
                 }
             }
         }
         pKc /= (pAreas);
         pKcw /= (pAreas);
         pKcCL /= pLCL;
		pAreas=pAreastmp;

     //~ forAll(pMarks, pI)    if(pMarks[pI]==4) pKc[pI]*=1.1;













	  scalarField pKcCLSmooth(pKcCL);
      for (int iKern=0; iKern<kernelRadius; ++iKern)
      {
			scalarField pKcCLSmoothTmp(pKcCLSmooth);
			forAll(pointPoints, vertI)
			{
				const labelLoop& neiPoints = pointPoints[vertI];
				label currentLabel = pMarks[vertI];
				if( currentLabel == 4 )
				{ ///. calculating smooth curvature  
					scalar avgAKcp(pKcCLSmoothTmp[vertI]);
					scalar sumWeightsKc(1.);
					forAll(neiPoints, neiPointI)
					{
						const label neiVertI = neiPoints[neiPointI];
						if( pMarks[neiVertI] == 4  )
						{
							scalar weight=1.;
							avgAKcp += weight * pKcCLSmoothTmp[neiVertI];
							sumWeightsKc += weight;
						}
					}
					if (sumWeightsKc>0.5)
						pKcCLSmooth[vertI] =     0.5*avgAKcp/sumWeightsKc + (1.-0.5)*pKcCLSmoothTmp[vertI];
					else
						pKcCLSmooth[vertI] = pKcCLSmoothTmp[vertI];
				}
		   }
	  }


	  scalarField pKcSmooth(pKc);
      for (int iKern=0; iKern<kernelRadius; ++iKern)
      {
			scalarField pKcSmoothTmp(pKcSmooth);
			forAll(pointPoints, vertI)
			{

				const labelLoop& neiPoints = pointPoints[vertI];
				label currentLabel = pMarks[vertI];
				if( currentLabel == 2 || currentLabel == 4  )
				{ ///. calculating smooth curvature  
					scalar avgAKcp(pKcSmoothTmp[vertI]);
					scalar sumWeightsKc(1.);
					forAll(neiPoints, neiPointI)
					{
						const label neiVertI = neiPoints[neiPointI];
						if( pMarks[neiVertI] == 2 || pMarks[neiVertI] == 4  )
						{
							scalar weight=pWeights[neiVertI];  weight*=weight;weight*=weight;
							if( pMarks[neiVertI] == 4) weight*=0.01;
							//~ if( pMarks[neiVertI] == 4 && currentLabel == 4) weight=10.;
							avgAKcp += weight * pKcSmoothTmp[neiVertI];
							sumWeightsKc += weight;
						}
					}
					if (sumWeightsKc>0.5)
						pKcSmooth[vertI] =     0.5*avgAKcp/sumWeightsKc + (1.-0.5)*pKcSmoothTmp[vertI];
					else
						pKcSmooth[vertI] = pKcSmoothTmp[vertI];
				}
		   }
	  }



	  scalarField pKcSmooth2(pKcSmooth);
      for (int iKern=0; iKern<kernelRadius+2; ++iKern)
      {
			scalarField pKcSmooth2Tmp(pKcSmooth2);
			forAll(pointPoints, vertI)
			{
				const labelLoop& neiPoints = pointPoints[vertI];
				label currentLabel = pMarks[vertI];
				if( currentLabel == 2 || currentLabel == 4  )
				{ ///. calculating smooth curvature  
					scalar avgAKcp(pKcSmooth2Tmp[vertI]);
					scalar sumWeightsKc(1.);
					forAll(neiPoints, neiPointI)
					{
						const label neiVertI = neiPoints[neiPointI];
						if( pMarks[neiVertI] == 2 || pMarks[neiVertI] == 4  )
						{
							scalar weight=pWeights[neiVertI];  weight*=weight;weight*=weight;
							if( pMarks[neiVertI] == 4) weight*=0.01;  ///.  TODO change 0.01 : 0-1
							//~ if( pMarks[neiVertI] == 4 && currentLabel == 4) weight=10.;
							avgAKcp += weight * pKcSmooth2Tmp[neiVertI];
							sumWeightsKc += weight;
						}
					}
					if (sumWeightsKc>0.5)
						pKcSmooth2[vertI] =     0.5*avgAKcp/sumWeightsKc + (1.-0.5)*pKcSmooth2Tmp[vertI];
					else
						pKcSmooth2[vertI] = pKcSmooth2Tmp[vertI];
				}

		   }
	  }







        vectorField displacementsRc(newPoints.size(), vector(0,0,0));
        forAll(pointPoints, vertI)
        {


            label currentLabel = pMarks[vertI];

             if( currentLabel == 4  )
             { ///. contact line curvature uniformization
                scalar DispCL = 1.*relaxCL*( -1.*(pKcSmooth[vertI]-pKcSmooth2[vertI]) + 0.9*(pKc[vertI]-pKcSmooth[vertI])+ 1.*(pKcCL[vertI]-pKcCLSmooth[vertI]))
                *Foam::sqrt(pAreas[vertI]/ (pKcSmooth[vertI]*pKcSmooth[vertI] +0.1/pAreasS[vertI]));
                newPoints[vertI] +=  DispCL*(pNs[vertI] - (pNs[vertI] & pNw[vertI])*pNw[vertI]);  ///.  TODO change 0.3 , -2. , etc		/// 0.3 >> bag impact...!!
             }
             else  if(currentLabel == 2)
             {  ///. curvature uniformization
                 newPoints[vertI] +=  pNs[vertI] * relax*(0.5*(pKc[vertI]-pKcSmooth[vertI])+8.*(pKcSmooth[vertI]-pKcSmooth2[vertI]))
												*Foam::sqrt(pAreas[vertI]/ (pKcSmooth[vertI]*pKcSmooth[vertI]  +0.1/pAreasS[vertI]));  ///.  TODO change -2. , etc
             }



            //~ if( currentLabel == 4  )
            //~ { ///. contact line curvature uniformization
               //~ vector DispCL =  relaxCL*(-0.05*(pKcSmooth[pI]-pKc[pI])+ 5.*(pKcCL[pI]-pKcCLSmooth[pI]))* pAreas[pI]*pNs[pI];
               //~ displacementsRc[pI] =  (DispCL - (DispCL & pNw[pI])*pNw[pI]);
            //~ }
            //~ else  if(currentLabel == 2)
            //~ {  ///. curvature uniformization
                //~ displacementsRc[pI] =  relax*(pKc[pI]-pKcSmooth[pI]) * pAreas[pI]*pNs[pI];
            //~ }
		}







/*
		///K_x.txt
		std::ofstream  offKcX("Kc_x.txt");
		offKcX<<"K 1 "<<pointPoints.size()<<" float"<<std::endl;
        forAll(pointPoints, vertI)
        {
			label currentLabel = pMarks[vertI];
            if(iter==iters-1)
            {
				// if (currentLabel==2||currentLabel==4)	
				if (currentLabel==2 && currentLabel!=4)	
				{
					offKcX<<"";
					offKcX<<points[vertI][0]<<" "<<points[vertI][1]<<" "<<points[vertI][2]<<" "<<pKc[vertI]<<" \n";
				}  
			}
       }		
       offKcX.close();
       
       	///pA_x.txt
		std::ofstream  ofpA("pA_x.txt");
		ofpA<<"pArea 1 "<<pointPoints.size()<<" float"<<std::endl;
        forAll(pointPoints, vertI)
        {
			label currentLabel = pMarks[vertI];
            if(iter==iters-1)
            {	
				if (currentLabel==1||currentLabel==3 ||currentLabel==4)		
				{
					ofpA<<"";
					ofpA<<points[vertI][0]<<" "<<points[vertI][1]<<" "<<points[vertI][2]<<" "<<pAreas[vertI]<<" \n";
				}  
			}
       }
       ofpA.close();		
		
		std::ofstream  offRc("Kc.txt");
		offRc<<"POINT_DATA "<<pointPoints.size()<<std::endl;
		offRc<<"FIELD attributes 7 "<<std::endl;
		offRc<<"Kc 1 "<<pointPoints.size()<<" float"<<std::endl;
        forAll(pointPoints, vertI)
        {
            if(iter==iters-1)
            {
				if ((vertI+1)%1==0)	offRc<<"";
				offRc<<pKc[vertI]<<" "<<"  \n";    
			}
       }
       offRc<<"Kcw 1 "<<pointPoints.size()<<" float"<<std::endl;
        forAll(pointPoints, vertI)
        {
            if(iter==iters-1)
            {
				if ((vertI+1)%1==0)	offRc<<"";
				offRc<<pKcw[vertI]<<" "<<"  \n";    
			}
       }
		offRc<<"\n\n RcCL 1 "<<pointPoints.size()<<" float"<<std::endl;
        forAll(pointPoints, vertI)
        {
            if(iter==iters-1)
            {
				if ((vertI+1)%1==0)	offRc<<"";
				offRc<<pKcCL[vertI]<<" "<<"  \n";    
			}
       }
		offRc<<"\n\n InterfNorms 3 "<<pointPoints.size()<<" float"<<std::endl;
        forAll(pointPoints, vertI)
        {
            if(iter==iters-1)
            {
				if ((vertI+1)%1==0)	offRc<<"";
				offRc<<pNs[vertI][0]<<" "<<pNs[vertI][1]<<" "<<pNs[vertI][2]<<" "<<" \n";    
			}
       }
		offRc<<"\n\n wallNorms 3 "<<pointPoints.size()<<" float"<<std::endl;
        forAll(pointPoints, vertI)
        {
            if(iter==iters-1)
            {
				if ((vertI+1)%1==0)	offRc<<"";
				offRc<<pNw[vertI][0]<<" "<<pNw[vertI][1]<<" "<<pNw[vertI][2]<<" "<<"  \n";    
			}
       }
       
       
		//offRc<<"\n\n InterfNormsCL 3 "<<pointPoints.size()<<" float"<<std::endl;
        //forAll(pointPoints, vertI)
        //{
			//label currentLabel = pMarks[vertI];
            //if(iter==iters-1)
            //{
				//if (currentLabel==4)		
				//{
					//if ((vertI+1)%1==0)	offRc<<"";
					//offRc<<pNs[vertI][0]<<" "<<pNs[vertI][1]<<" "<<pNs[vertI][2]<<" "<<" \n";    
				//}
			//}
       //}
		//offRc<<"\n\n wallNormsCL 3 "<<pointPoints.size()<<" float"<<std::endl;
        //forAll(pointPoints, vertI)
        //{
            //label currentLabel = pMarks[vertI];
            //if(iter==iters-1)
            //{
				//if (currentLabel==4)		
				//{
					//if ((vertI+1)%1==0)	offRc<<"";
					//offRc<<pNw[vertI][0]<<" "<<pNw[vertI][1]<<" "<<pNw[vertI][2]<<" "<<"  \n";    
				//}
			//}
       //}
		offRc.close();
		
		//~ newPoints += displacementsRc;
        
*/

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
		if ( pMarks[vertI] == 4 )
		{
			forAll(neiPoints, neiPointI)
			{
				const label neiVertI = neiPoints[neiPointI];
			
				scalar weight=(pAreas[neiVertI]+0.1*pAreas[vertI]);//*max(NormalsOld[neiVertI]&NormalsOld[vertI],0.001);// /mag(newPoints[neiVertI]-newPoints[vertI]);
				weight*=weight*pWeights[neiVertI];
				if( pMarks[neiVertI]!=currentLabel ) weight*=0.1;  ///.  TODO change 0.5,   , affects convergence 
				point neiPos=newPoints[neiVertI];
				if( pMarks[neiVertI]==2 )
				{
					neiPos -= 1.01*((neiPos-newPoints[vertI])&pNw[vertI])*pNw[vertI];    ///.  TODO change 1.05
				}
				else if( pMarks[neiVertI] == 4) weight *= 8.; ///. TODO: Change 3.333, remove ...
				avgPos += weight * neiPos;
				sumWeights += weight;
			}

		}
		else
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


     for (int iKern=0; iKern<2; ++iKern)
     {
 	  vectorField displacementsTmp = displacements;
      forAll(pointPoints, vertI)
      {
		label currentLabel = pMarks[vertI];
		bool curentLabelEq2 = currentLabel==2;
		//~ if( currentLabel != 2  )
		{ ///. contact line smoothing Vol preserve
			const labelLoop& neiPoints = pointPoints[vertI];
			vector avgDisp(vector::zero);
			scalar sumWeights(1e-28); 
			forAll(neiPoints, iii)
			{
				const label neiVertI = neiPoints[iii];
				//~ if( pMarks[neiVertI] == currentLabel  )
				if( pMarks[neiVertI]==currentLabel || !( (pMarks[neiVertI]==2) ^ curentLabelEq2)  || pMarks[neiVertI] == 4  )
				{
					scalar weight=pAreas[neiVertI];
					avgDisp += weight * displacementsTmp[neiVertI];
					sumWeights += weight;
				}
			}

			displacements[vertI] = avgDisp/sumWeights;//myEdges.size();

		}
		//~ else displacements[vertI] = vector::zero;
		//~ if( currentLabel == 4  ) displacements[vertI] = vector::zero;

      }
     }

	  //~ newPoints -= 0.3*displacements+0.79*(displacements& pNw)*pNw;     ///.  TODO change 0.5 0.5,  sum should be 1

	displacements = 0.3*displacements+0.7*(displacements& pNw)*pNw;
	displacements = 0.3*displacements+0.7*(displacements& pNs)*pNs;
	newPoints -= 1.09*displacements;     ///.  TODO change 0.5 0.5,  sum should be 1
	
		Info  <<"CL "; cout.flush();






	surf123.movePoints(newPoints);


	displacements = newPoints-previousPoints;

	Info<< "  A:  ["<<min(pAreas)<<" "<<max(pAreas)<<"] "; cout.flush();
	Info<< "  kw:  ["<<min(pKcw); cout.flush();
	//~ Info<< " "<< 0.004+(double(average(max(min((pKc-0.004)*1e8,1e6),-1e6)*max(pMarks-3,0))))/(0.0001+double(average(100000*max(pMarks-3,0))))/1000.; cout.flush();
	Info<< " "<<max(pKcw); cout.flush();
	Info<< "],  DispMag  max:  "<<max(mag(displacements)); cout.flush();
	Info<< "  avg:  "<<average(mag(displacements)); cout.flush();
	Info<<endl;

