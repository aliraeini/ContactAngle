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


     //vectorField fNorms(faces.size(),vector(0,0,0));
     //forAll(faces,faceI)
     //{   const face & f=faces[faceI];
         //forAll(f, pI)
         //{    label pCur=f[pI];
             //fNorms[faceI] += (newPoints[pCur]-newPoints[f.prevLabel(pI)])^(newPoints[f.nextLabel(pI)]-newPoints[pCur]);
     //}    }
     //fNorms=fNorms/(mag(fNorms)+1e-12);

	/// Compute point-normal vectors
     vectorField pAw(pNeips.size(),vector(0,0,0));
     vectorField pAs(pNeips.size(),vector(0,0,0));
     forAll(faces,faceI)
     {   const face & f=faces[faceI];
         forAll(f, pI)
         {   edge e=f.faceEdge(pI);
             vector Ce=0.5*(newPoints[e[0]] + newPoints[e[1]]);

             vector AEdgeTri=0.5*((newPoints[e[0]]-Cf[faceI]) ^ (Ce-Cf[faceI]));

				 if( pMarks[e[0]]==fMarks[faceI] )
				 {
					 pAw[e[0]] += AEdgeTri;
					 pAs[e[0]] += AEdgeTri;
				 }
				 else if( fMarks[faceI] != 2 )
				 {
					 pAw[e[0]] += AEdgeTri;
				 }
				 else
					 pAs[e[0]] += AEdgeTri;

				AEdgeTri=0.5*((newPoints[e[1]]-Cf[faceI]) ^ (Ce-Cf[faceI]));

				 if( pMarks[e[1]]==fMarks[faceI] )
				 {
					 pAw[e[1]] -= AEdgeTri;
					 pAs[e[1]] -= AEdgeTri;
				 }
				 else if( fMarks[faceI] != 2 )
					 pAw[e[1]] -= AEdgeTri;
				 else
					 pAs[e[1]] -= AEdgeTri;//2.*AEdgeTri-mag(AEdgeTri)*fNorms[faceI];
		  }
	  }
     //scalarField pAreas(0.5*mag(pAs)+1e-16);
     scalarField pAreas(0.5*mag(pAs)+1e-18+0.5*mag(pAw));
     scalarField pAreastmp(0.5*mag(pAs)+1e-18+0.5*mag(pAw));
     vectorField pNs=pAs/(mag(pAs)+1e-18);



	/// Correct point-normal vectors at contact line (BC), using linear extrapolation
	forAll(pNeips, vertI)
	{
		const labelLoop& neiPoints = pNeips[vertI];
		if( pMarks[vertI] == 4  ) //. only for contact line
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
