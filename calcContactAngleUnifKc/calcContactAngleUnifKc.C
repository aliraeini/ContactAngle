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
    Volume-preserving Gaussian and curvature uniform smoothing and contact angle measurement

\*---------------------------------------------------------------------------*/

//~ #include "triSurface.H"
#include "MeshedSurfaces.H"

#include "argList.H"
#include "OFstream.H"
#include "boundBox.H"

using namespace Foam;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:
int appendUnique(DynamicList<label>& dynList, label value)
{
    bool append=true;
    forAll(dynList,i)
    {
        if(dynList[i]==value) append=false;
    }
    if(append)
    {
        dynList.append(value);
        return 1;
    }
    else
    {
        return 0;
    }
}


int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"

runTime++;
Info<< "Time = " << runTime.timeName() << "\n" << endl;


    IOdictionary meshingDict
    (
        IOobject
        (
            "meshingDict",
            runTime.system(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    dictionary smoothingDict(meshingDict.subDict("surfaceSmoothing"));

    word surfFileName(smoothingDict.lookup("inputSurface"));

    fileName outFileName(smoothingDict.lookup("outputSurface"));






    Info<< "Reading surface from " << surfFileName << " ..." << endl;

    meshedSurface surf123(surfFileName);


    Info<< "nFaces    : " << surf123.size() << endl;
    Info<< "nVertices     : " << surf123.nPoints() << endl;
    Info<< "nZones     : " << surf123.surfZones().size() << endl;
    Info<< "Bounding Box : " << boundBox(surf123.localPoints()) << endl;




    #define labelLoop face

    const pointField & points=surf123.points();
    const List<face> & faces=surf123.faces();
    DynamicList<DynamicList<label> > pointPointsTmp(points.size());
    List<labelLoop> pNeips(points.size()); //< adjacent nodes to each node
    forAll(faces,fI)
    {
        const face & f=faces[fI];
        forAll(f, pI)
        {
           appendUnique(pointPointsTmp[ f[pI] ], f.nextLabel(pI));
           appendUnique(pointPointsTmp[ f[pI] ], f.prevLabel(pI));
        }
    }

    forAll(pNeips,i)
    {
        //pNeips[i].setSize(pointPointsTmp[i].size());
        pNeips[i]=face(pointPointsTmp[i]);
    }

    pointField newPoints(points);

    Info<<"\n........"<<endl<<endl ;
//~ exit(-1);










		/// rock-water interface -> 1
		/// rock-oil   interface -> 3
		/// oil-water interface -> 2
		/// contact line points ->  4

    labelField fMarks(faces.size(),1);  /// 1; GrainWalls, 2: OW interface
    labelField pMarks(pNeips.size(),2); /// 1; GrainWalls and contactline, 2: OW interface
    scalarField pWeights(pNeips.size(),1.); /// 1; GrainWalls and contactline, 2: OW interface

	{
		List< surfZone >   zones = surf123.surfZones(); for(int i=zones.size();i<4;++i) {zones.append(surfZone()); zones[zones.size()-1].start()=faces.size();}
		for (label i=zones[0].start(); i<zones[1].start(); ++i)     fMarks[i]=1;
		for (label i=zones[1].start(); i<zones[2].start(); ++i)     fMarks[i]=2;
		for (label i=zones[2].start(); i<zones[3].start(); ++i)     fMarks[i]=3;


		for (label i=zones[0].start(); i<zones[1].start();  ++i)  forAll(faces[i],pi) pMarks[faces[i][pi]]=1;
		for (label i=zones[2].start(); i<zones[3].start();  ++i)  forAll(faces[i],pi) pMarks[faces[i][pi]]=3;
		for (label i=zones[1].start(); i<zones[2].start();  ++i)  forAll(faces[i],pi) if (pMarks[faces[i][pi]]!=2) pMarks[faces[i][pi]]=4;

		Info<< "faceZones     : " << min(fMarks)<<" - "<< max(fMarks) << endl;
		Info<< "pMarks     : " << min(pMarks)<<" - "<< max(pMarks)<< endl;
		Info<<"........"<<endl<<endl ;





      for (int ii=0; ii<3; ++ii)
      {
		for (label i=zones[1].start(); i<zones[2].start();  ++i)  forAll(faces[i],pi) if (pMarks[faces[i][pi]]!=2) { pWeights[faces[i][pi]]=0.25;   }

		scalarField pWeightsOrig(pWeights);
		forAll(pNeips, vertI)
		{

			const labelLoop& neiPoints = pNeips[vertI];

			scalar avgPWeights(pWeightsOrig[vertI]);
			scalar sumWeightsKc(1.);
			forAll(neiPoints, neiPointI)
			{
				const label neiVertI = neiPoints[neiPointI];

				scalar weight=1.;
				avgPWeights += weight * pWeightsOrig[neiVertI];
				sumWeightsKc += weight;
			}
			if (sumWeightsKc>0.5)
			{
				avgPWeights /= sumWeightsKc;//myEdges.size();
				pWeights[vertI] =     0.99*avgPWeights + (1.-0.99)*pWeightsOrig[vertI]; ///. 0.3 affects convergence,  smaller value makes cl faces bigger
			}
			else
			{
				pWeights[vertI] = pWeightsOrig[vertI];
			}

		}
	  }

	}








	scalar extraRelaxFactorParIntrf(readScalar(smoothingDict.lookup("extraRelaxFactorParIntrf")));




///   ///////////////////////////////////////////////////////


	{
		Info<< "volume preserving Gaussian smoothing"<< endl;
		Info<< "Reading keywords from system/meshingDict -> surfaceSmoothing:"<< endl;
		label nIters(readLabel(smoothingDict.lookup("nIterationsGaussVP")));
		label kernelRadius(readLabel(smoothingDict.lookup("kernelRadiusGaussVP")));
		scalar relax(readScalar(smoothingDict.lookup("relaxFactorGaussVP")));
		scalar relaxCL(readScalar(smoothingDict.lookup("relaxFactorGaussVPCL")));
		if ((relax < 0) || (relax > 1))
				FatalErrorIn(args.executable()) << "Illegal relaxation factor " << relax
				 << endl << "0: no change 1: move vertices to average of neighbours" << exit(FatalError);
		Info<< "   nIterationsGaussVP:  " << nIters << endl;
		Info<< "   kernelRadiusGaussVP: " << kernelRadius << endl;
		Info<< "   relaxFactorGaussVP:  " << relax << endl;
		Info<< "   relaxFactorGaussVPCL:" << relaxCL << endl;


		for(label iter = 0; iter < nIters; iter++)
		{
			#include "./gauss.H"
		}
		Info<<"Gasss smoothed  :/"<<endl<<endl;
	}//===================================================================================



	{
		Info<< "curvature move/smoothing 1"<< endl;
		Info<< "Reading keywords from system/meshingDict -> surfaceSmoothing:"<< endl;
		label nIters(readLabel(smoothingDict.lookup("nIterationsCurvature1")));
		label kernelRadius(readLabel(smoothingDict.lookup("kernelRadiusCurvature1")));
		scalar relax(readScalar(smoothingDict.lookup("relaxFactorCurvature1")));
		scalar relaxCL(readScalar(smoothingDict.lookup("relaxFactorCurvature1CL")));

		if ((relax < 0) || (relax > 1))
				FatalErrorIn(args.executable()) << "Illegal relaxation factor " << relax
				 << endl << "0: no change 1: move vertices to average of neighbours" << exit(FatalError);
		if ((relaxCL < 0) || (relaxCL > 1))
				FatalErrorIn(args.executable()) << "Illegal CL relaxation factor " << relaxCL
				 << endl << "0: no change 1: move vertices to average of neighbours" << exit(FatalError);

		Info<< "  nIterationsCurvature1:   " << nIters << endl;
		Info<< "  kernelRadiusCurvature1:  " << kernelRadius << endl;
		Info<< "  relaxFactorCurvature1:   " << relax << endl;
		Info<< "  relaxFactorCurvature1CL: " << relaxCL << endl;


		for(label iter = 0; iter < nIters; iter++)
		{
			#include "./curvatureMove.H"
		}
		Info<<"smoothed  :/"<<endl<<endl;
	}//===================================================================================

	{
		Info<< "curvature move/smoothing, second round"<< endl;
		label  nIters(      readLabel (smoothingDict.lookup("nIterationsCurvature")));
		label  kernelRadius(readLabel (smoothingDict.lookup("kernelRadiusCurvature")));
		scalar relax(       readScalar(smoothingDict.lookup("relaxFactorCurvature")));
		scalar relaxCL(     readScalar(smoothingDict.lookup("relaxFactorCurvatureCL")));

		if ((relax < 0) || (relax > 1))
				FatalErrorIn(args.executable()) << "Illegal relaxation factor " << relax
				 << endl << "0: no change 1: move vertices to average of neighbours" << exit(FatalError);
		if ((relaxCL < 0) || (relaxCL > 1))
				FatalErrorIn(args.executable()) << "Illegal CL relaxation factor " << relaxCL
				 << endl << "0: no change 1: move vertices to average of neighbours" << exit(FatalError);

		Info<< "  nIterationsCurvature:   " << nIters << endl;
		Info<< "  kernelRadiusCurvature:  " << kernelRadius << endl;
		Info<< "  relaxFactorCurvature:   " << relax << endl;
		Info<< "  relaxFactorCurvatureCL: " << relaxCL << endl;


		for(label iter = 0; iter < nIters; iter++)
		{
			#include "./curvatureMove.H"
		}
		Info<<"smoothed  :/"<<endl<<endl;
	}//===================================================================================




	Info<<average(mag(surf123.points()))<< "  :avgMagPoints    "<<endl;
	Info<<max(mag(surf123.points()))<< "  :maxMagPoints    "; cout.flush();



    Info<< "Writing surface to " << outFileName << " ..." << endl;


	#include "./calcContAngle.H"
	//~ #include "./calcCurvature.H"



    surf123.write(outFileName);


	Info << "\n"<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
	<< " ClockTime = " << runTime.elapsedClockTime() << " s"
	<< "\n" << endl;

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
