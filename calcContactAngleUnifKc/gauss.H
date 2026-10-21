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
    Volume-preserving Gaussian smoothing

\*---------------------------------------------------------------------------*/




///. point-normal vectors
pointField Cf=surf123.faceCentres();

vectorField pNw(pNeips.size(),vector(0,0,0));
vectorField pNs(pNeips.size(),vector(0,0,0));
	forAll(faces,fI)  {
		const face & f=faces[fI];
		forAll(f, pI)  {
			edge e=f.faceEdge(pI);
			vector Ce=0.5*(newPoints[e[0]] + newPoints[e[1]]);

            vector pNE=0.5*((newPoints[e[0]]-Cf[fI]) ^ (Ce-Cf[fI]));
            if( pMarks[e[0]]==fMarks[fI] ) {    pNw[e[0]] += pNE;    pNs[e[0]] += pNE;    }
            else if( fMarks[fI] != 2 )        pNw[e[0]] += pNE;
            else                                pNs[e[0]] += pNE;

            pNE=0.5*((newPoints[e[1]]-Cf[fI]) ^ (Ce-Cf[fI]));
            if( pMarks[e[1]]==fMarks[fI] ) {    pNw[e[1]] -= pNE;    pNs[e[1]] -= pNE;    }
            else if( fMarks[fI] != 2 )        pNw[e[1]] -= pNE;
            else                                pNs[e[1]] -= pNE;//2.*pNE-mag(pNE)*fNorms[fI];
		}
	}
     scalarField pAreas(0.5*mag(pNs)+1e-18+0.5*mag(pNw));
pNs/=(mag(pNs)+1e-18);
pNw/=(mag(pNw)+1e-18);

const vectorField previousPoints(newPoints);



//===========================================================================
Info<< iter   <<"  "; cout.flush();
for (int iKern=0; iKern<kernelRadius; ++iKern)  {
	  vectorField displacements(newPoints.size(), vector(0,0,0));
	  vectorField displacementsCL(newPoints.size(), vector(0,0,0));
	  forAll(pNeips, pI)  {

		const labelLoop& neiPoints = pNeips[pI];
		{
			vector avgPos(vector::zero);
			scalar sumWeights(1e-16);
			if ( pMarks[pI] == 4 )
			  for(label pJ:neiPoints)  {
				if( pMarks[pJ]!=2 ) {
					//~ scalar weight=pWeights[pJ];//*max(NormalsOld[pJ]&NormalsOld[pI],0.001);// /mag(newPoints[pJ]-newPoints[pI]);
					scalar weight=(pAreas[pJ]+0.3*pAreas[pI]);//*max(NormalsOld[pJ]&NormalsOld[pI],0.001);// /mag(newPoints[pJ]-newPoints[pI]);
					weight*=weight*pWeights[pJ];
					avgPos += weight * newPoints[pJ];
					sumWeights += weight;
				}
			 }
		  else for(label pJ:neiPoints)  {
					//~ scalar weight=pWeights[pJ];//*max(NormalsOld[pJ]&NormalsOld[pI],0.001);// /mag(newPoints[pJ]-newPoints[pI]);
					scalar weight=(pAreas[pJ]+0.3*pAreas[pI]);//*max(NormalsOld[pJ]&NormalsOld[pI],0.001);// /mag(newPoints[pJ]-newPoints[pI]);
					weight*=weight*pWeights[pJ];
					avgPos += weight * newPoints[pJ];
					sumWeights += weight;
			  }
			if (sumWeights>0.1) {
				avgPos /= sumWeights;//myEdges.size();
				displacements[pI] = avgPos-newPoints[pI]; //(0.8*relax)*((avgPos-newPoints[pI])&pNw[pI])*pNw[pI] + (0.2*relax)*(avgPos-newPoints[pI]);
			}
		}



		if( pMarks[pI] == 4 ) { ///. contact line smoothing
			vector avgPos(vector::zero);
			scalar sumWeights(1e-16);
			for(label pJ: neiPoints)
			 if( pMarks[pJ] == 4  ) {
						scalar weight=1.;
						avgPos += weight * newPoints[pJ];
						sumWeights += weight;
			 }

            vector DispCL =  (avgPos/sumWeights-newPoints[pI]); // relax*((avgPos/sumWeights-newPoints[pI])&pNs[pI])*pNs[pI];
            DispCL =  (DispCL & pNs[pI])*pNs[pI];
            displacementsCL[pI] =  (DispCL - (DispCL & pNw[pI])*pNw[pI]);

        }

    }

	Info  <<"Gauss "; cout.flush();
	newPoints += (1.-extraRelaxFactorParIntrf)*relax*(displacements & pNw)*pNw + (extraRelaxFactorParIntrf*relax)*(displacements);
	newPoints += relaxCL*displacementsCL;
}




Info <<"   Preserve Volume, "; cout.flush();

vectorField dispAverage(newPoints-previousPoints);

for (int iKern=0; iKern<kernelRadius+2; ++iKern)  {
	const vectorField dispAverageTmp(dispAverage);

	forAll(pNeips, pI) {
		vector avgDisp(vector::zero);
		scalar sumWeights(1e-16);

		const labelLoop& neiPoints = pNeips[pI];
		label currentLabel = pMarks[pI];
		bool curentLabelEq2 = currentLabel==2;

		for(label pJ:neiPoints)  {

			if( pMarks[pJ]==currentLabel || !( (pMarks[pJ]==2) ^ curentLabelEq2)  || pMarks[pJ] == 4  ) {
				scalar weight=pAreas[pJ];
				if( pMarks[pJ]!=currentLabel ) weight*=0.05;
				avgDisp += weight * dispAverageTmp[pJ];
				sumWeights += weight;
			}
		}
		dispAverage[pI] = avgDisp/sumWeights;//myEdges.size();
	}
}


dispAverage -= 0.3*(dispAverage-(dispAverage& pNw)*pNw);
dispAverage -= 0.3*(dispAverage-(dispAverage& pNs)*pNs);
newPoints -= 1.09*dispAverage;     ///.  TODO change 0.5 0.5,  sum should be 1



surf123.movePoints(newPoints);


dispAverage=(newPoints-previousPoints);

Info<< "  Disp max :  "<<max(mag(dispAverage));
Info<< "  avg:  "<<average(mag(dispAverage));
//Info<< "  avg:  "<<average(dispAverage);
Info<< "  Amin:  "<<min(pAreas);
Info<<endl;
