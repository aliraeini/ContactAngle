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
    Curvature calculation

\*---------------------------------------------------------------------------*/



/*

POINT_DATA 75466
FIELD attributes 2
norms 3 75466 float


Rc 1 75466 float


=====================================


POINT_DATA 75466
FIELD attributes 1
Rc 1 75466 float


*/





    scalarField curvature(fNorms.size(),1e-12);
     ///. TODO to define normal vectors separately
    forAll(faces,faceI)
    {
        const face & f=faces[faceI];
        forAll(f, pI)
        {
			if( pMarks[f[pI]] != 4 ||  fMarks[faceI] == 1 )
			{
				edge e=f.faceEdge(pI);
			    vector Le_se = 0.5*(pNw[e[0]] + pNw[e[1]]) ^ e.vec(newPoints);
			    curvature[faceI] += Le_se & fNorms[faceI];
		    }
        }
		curvature[faceI] /= f.mag(newPoints);
    }





(Info<< " *  ").flush();

		std::ofstream  offCurvature("curvatures.txt");
        for (int i=0; i<curvature.size(); ++i)
        {
				if ((i+1)%20==0)	offCurvature<<"\n";
				offCurvature<<curvature[i]<<" ";
		} offCurvature.close();


(Info<< " *  ").flush();
