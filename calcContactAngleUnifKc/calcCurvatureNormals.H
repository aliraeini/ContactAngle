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
			{	scalar weight=1.;
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

