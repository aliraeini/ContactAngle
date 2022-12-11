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
    Constrain each vertex of the 3-phase contact line to have a single edge connections with its neighbors.

\*---------------------------------------------------------------------------*/




#include "argList.H"
#include "OFstream.H"
#include "boundBox.H"
#include "DynamicField.H"
#include "MeshedSurfaces.H"

using namespace Foam;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"

    IOdictionary meshingDict
    (  IOobject ( "meshingDict", runTime.system(),  runTime, IOobject::MUST_READ  ) );
    dictionary smoothingDict(meshingDict.subDict("addLayerToCL"));

    word surfFileName(smoothingDict.lookup("inputSurface"));

    fileName outFileName(smoothingDict.lookup("outputSurface"));


    Info<< "Reading surface from " << surfFileName << " ..." << endl;

    meshedSurface surf123(surfFileName);

	const faceList &   faces = surf123.faces();
	const pointField &   points = surf123.points();


    Info<< "nFaces    : " << surf123.size() << endl;
    Info<< "nVertices     : " << surf123.nPoints() << endl;
    Info<< "nZones     : " << surf123.surfZones().size() << endl;
    Info<< "Bounding Box : " << boundBox(surf123.localPoints()) << endl;





	/// rock-water interface -> 1
	/// rock-oil   interface -> 3
	/// oil-water interface -> 2
	/// contact line points ->  4

    labelField fMarks(faces.size(),1);  /// 1; GrainWalls, 2: OW interface
    labelField pMarks(points.size(),2); /// 1; GrainWalls and contactline, 2: OW interface

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
	}





	const edgeList & 	edges = surf123.edges();
	labelList 	edgesNInternalFaces(edges.size(),0);
	labelList 	pNewPsP3(points.size(),-2);
	const labelList& mps = surf123.meshPoints();



	Info<<"  marking new points" <<endl;

	///. generate a list of new points and faces to be added:
	///. 3 points for each contact line point and 3 faces for
	///. each contact line edge
	int iNewP=-1, iNewF=-1;
	forAll(edges, ei)
	  if (pMarks[mps[edges[ei].start()]]==4 && pMarks[mps[edges[ei].end()]]==4)
	  {
		++iNewF;
		edge e=edges[ei];
		if (pNewPsP3[mps[e.end()]  ]<0) pNewPsP3[mps[e.end()]  ]=++iNewP;
		if (pNewPsP3[mps[e.start()]]<0) pNewPsP3[mps[e.start()]]=++iNewP;
	  }///\\\///\\\///\\\///\\\///
  		Info<<endl<<iNewP*3+3<<" new points,  "<<iNewF*3+3<<" new faces"<<endl;

	///. add 3 faces for each edge in contact line
	faceList  newfacs((iNewF+1), face(4));
	iNewF=-1;
	forAll(edges, ei)
	  if (pMarks[mps[edges[ei].start()]]==4 && pMarks[mps[edges[ei].end()]]==4)
	  {
		++iNewF;
		edge e=edges[ei];
		newfacs[iNewF][0]=mps[e.end()];///. store two of the vertex indices
		newfacs[iNewF][1]=mps[e.start()];///. store two of the vertex indices
		newfacs[iNewF][2]=-1;           ///. take care of later
		newfacs[iNewF][3]=-1;           ///. take care of later
	  }///\\\///\\\///\\\///\\\///







    { ///. add newfacs to surf123  faces

		faceList facesSorted_All(faces.size()+3*newfacs.size());
		labelList zoneSizes(3);

		///. detach old faces of zone 1 from the contact line  and connect to new vertex indices
		int indF=-1;
		for (int i=0; i<faces.size(); ++i) /// water-rock interface
		  if(fMarks[i]==1)
		  {
			face f=faces[i];
			forAll(f, fi)	if(pNewPsP3[f[fi]]>=0) f[fi] = points.size()+pNewPsP3[f[fi]]*3+0; ///. detach and connect to new vertex indices
			facesSorted_All[++indF] = f;  ///. store the face into a new list
		  };
		  ///. take care of the vertex indices of the new face to be added to zone 1
		for (int i=0; i<newfacs.size(); ++i)
		{
			face f=newfacs[i];

			f[2] = f[1];  /// contact line edge vertices
			f[3] = f[0];  /// contact line edge vertices
			f[1] = points.size()+pNewPsP3[f[1]]*3+0;  ///. new vertex index
			f[0] = points.size()+pNewPsP3[f[0]]*3+0;  ///. new vertex index
			facesSorted_All[++indF] = f;
		}
		zoneSizes[0]=indF+1;

		for (int i=0; i<faces.size(); ++i) /// water-oil interface		///. same as above, zone 2
		  if(fMarks[i]==2)
		  {
			face f=faces[i];
			forAll(f, fi)	if(pNewPsP3[f[fi]]>=0) f[fi] = points.size()+pNewPsP3[f[fi]]*3+1;
			facesSorted_All[++indF] = f;
		  };
		for (int i=0; i<newfacs.size(); ++i)
		{
			face f=newfacs[i];
			f[2] = points.size()+pNewPsP3[f[1]]*3+1;
			f[3] = points.size()+pNewPsP3[f[0]]*3+1;
			facesSorted_All[++indF] = f;
		}
		zoneSizes[1]=indF+1-zoneSizes[0];


		for (int i=0; i<faces.size(); ++i)  /// oil-rock interface		///. same as above, zone 3
		  if(fMarks[i]==3)
		  {
			face f=faces[i];
			forAll(f, fi)	if(pNewPsP3[f[fi]]>=0) f[fi] = points.size()+pNewPsP3[f[fi]]*3+2;
			facesSorted_All[++indF] = f;
		  };
		for (int i=0; i<newfacs.size(); ++i)
		{
			face f=newfacs[i];
			f[2] = points.size()+pNewPsP3[f[1]]*3+2;
			f[3] = points.size()+pNewPsP3[f[0]]*3+2;
			facesSorted_All[++indF] = f;
		}
		zoneSizes[2]=indF+1-zoneSizes[1]-zoneSizes[0];

  		Info<<faces.size()<<" fs?<= "<<facesSorted_All.size()<<endl;
  		Info<<indF+1-faces.size()<<" fs?= "<<facesSorted_All.size()-faces.size()<<endl;
  		Info<<iNewP*3+3<<"new points,  "<<indF+1-faces.size()<<" new faces"<<endl;

		///. create new point positions, 3 for each conact line vertex
		Field<point> pointsNew((iNewP+1)*3, point(100.,100.,100.));
		forAll(points, pi)
		  if(pNewPsP3[pi]>=0)
		  {
			  pointsNew[pNewPsP3[pi]*3]=points[pi];
			  pointsNew[pNewPsP3[pi]*3+1]=points[pi];
			  pointsNew[pNewPsP3[pi]*3+2]=points[pi];
		  }

		DynamicField<point> pointsAll;
		pointsAll.reserve(points.size()+pointsNew.size());
		pointsAll.append(points);
		pointsAll.append(pointsNew);
  		Info<<points.size()<<" ps?<= "<<pointsAll.size()<<"      "<<faces.size()<<" fs?<= "<<facesSorted_All.size()<<endl;

		/// convert to openfoam weired data pointers
        Xfer<List<face> > facesFer(facesSorted_All,true);
        Field<point> & pointsSF=pointsAll;
		Xfer<Field<point> > pointsFer(pointsSF,true);


        meshedSurface surfOut(pointsFer, facesFer, zoneSizes);

        surfOut.write(outFileName);

    }


    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
