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
    Normal vectors calculations

\*---------------------------------------------------------------------------*/



     pointField Cf=surf123.faceCentres();


     //~ vectorField fNorms(faces.size(),vector(0,0,0));
     //~ forAll(faces,fI)
     //~ {   const face & f=faces[fI];
         //~ forAll(f, pI)
         //~ {    label pCur=f[pI];
             //~ fNorms[fI] += (newPoints[pCur]-newPoints[f.prevLabel(pI)])^(newPoints[f.nextLabel(pI)]-newPoints[pCur]);
     //~ }    }
     //~ fNorms=fNorms/(mag(fNorms)+1e-12);

///. point-normal vectors
     vectorField pAw(pointPoints.size(),vector(0,0,0));
     vectorField pAs(pointPoints.size(),vector(0,0,0));
     forAll(faces,fI) {
         const face & f=faces[fI];
         forAll(f, pI) {
            edge e=f.faceEdge(pI);
             vector Ce=0.5*(newPoints[e[0]] + newPoints[e[1]]);

             vector AEdgeTri=0.5*((newPoints[e[0]]-Cf[fI]) ^ (Ce-Cf[fI]));

				 if( pMarks[e[0]]==fMarks[fI] ) {
					 pAw[e[0]] += AEdgeTri;
					 pAs[e[0]] += AEdgeTri;
				 }
				 else if( fMarks[fI] != 2 )  {
					 pAw[e[0]] += AEdgeTri;
				 }
				 else
					 pAs[e[0]] += AEdgeTri;

				AEdgeTri=0.5*((newPoints[e[1]]-Cf[fI]) ^ (Ce-Cf[fI]));

				 if( pMarks[e[1]]==fMarks[fI] )
				 {
					 pAw[e[1]] -= AEdgeTri;
					 pAs[e[1]] -= AEdgeTri;
				 }
				 else if( fMarks[fI] != 2 )
					 pAw[e[1]] -= AEdgeTri;
				 else
					 pAs[e[1]] -= AEdgeTri;//2.*AEdgeTri-mag(AEdgeTri)*fNorms[fI];
		  }
	  }
     //~ scalarField pAreas(0.5*mag(pAs)+1e-16);
     scalarField pAreas(0.5*mag(pAs)+1e-18+0.5*mag(pAw));
     scalarField pAreastmp(0.5*mag(pAs)+1e-18+0.5*mag(pAw));
     vectorField pNs=pAs/(mag(pAs)+1e-18);


///. point-normal vectors BC correction
      //{
          //vectorField pNNormB(pointPoints.size(),vector(0,0,0));
          //forAll(faces,fI)
          //{   const face & f=faces[fI];
              //vector pNsS(0.,0.,0.);
              //forAll(f, pI)
              //{    if( fMarks[fI]!=2 && pMarks[f[pI]]==fMarks[fI] )
                      //pNsS+=pNs[f[pI]];
              //}
              //forAll(f, pI)
              //{
                  //if( fMarks[fI]!=2 && pMarks[f[pI]]!=fMarks[fI] )
                      //pNNormB[f[pI]]+= pNsS*f.mag(points);
              //}
          //}
          //pNs-=pNNormB/(3.*mag(pNNormB)+1e-12); ///.  TODO change 7. theoretical is 2 or 4, to check.
          //pNs/=(mag(pNs)+1e-18);
      //}

	///. point-normal vectors correction at contact line (BC), using linear extrapolation
	forAll(pointPoints, pI)
	{
		const labelLoop& neiPoints = pointPoints[pI];
		if( pMarks[pI] == 4  ) ///. only for contact line
		{
			vector  pNsNei(0,0,0);
			for(label pJ: neiPoints) {
				if( pMarks[pJ] == 2 )
					pNsNei += pNs[pJ];
			}
			if (mag(pNsNei)>0.1)
			{
				pNsNei /= mag(pNsNei);
				pNs[pI]=1.333*pNs[pI]-0.333*pNsNei;
				pNs[pI] /= mag(pNs[pI]);
			}
		}
   }
