/*-------------------------------------------------------------------------*\
This code is part of poreFOAM, a suite of codes written using OpenFOAM
for direct simulation of flow at the pore scale.
You can redistribute this code and/or modify this code under the
terms of the GNU General Public License (GPL) as published by the
Free Software Foundation, either version 3 of the License, or (at
your option) any later version. see <http://www.gnu.org/licenses/>.


The code has been developed by Ali Qaseminejad Raeini as a part his PhD
at Imperial College London, under the supervision of Branko Bijeljic
and Martin Blunt.

Please see our website for relavant literature:
https://www.imperial.ac.uk/earth-science/research/research-groups/pore-scale-modelling/

Developed by:
- Ali Q Raeini (2010-2025)
- Ahmed AlRarout (2015?-2018)

 Description:
	creates a surface between two fluid phases and solid phase in a 3D image
\*-------------------------------------------------------------------------*/

using namespace Foam; // !!!

void correct( faceList & faces, labelList& fMarks, DynamicField<point> & points, bool handlemultipliConnectedEdges  );
void correctbioti( faceList & faces, labelList& fMarks, DynamicField<point> & points, int stage );
void correctbioti2( faceList & faces, labelList& fMarks, DynamicField<point> & points, int stage );



void writeSTLBINARY( const voxelImage & vxlImg, std::string outputSurface)
{


	Info<<"writeSTLBINARY: "<<endl;
	int3 n=vxlImg.size3();
	int nx_2=n[0]-2, ny_2=n[1]-2, nz_2=n[2]-2;
	dbl3 dx=vxlImg.dx();
	dbl3 X1=vxlImg.X0() + dx;



//=======================================================================




    int nInternalFaces=0;
    int nAllFaces=0;
    int nLeftFaces=0;
    int nRightFaces=0;
    int nBottomFaces=0;
    int nTopFaces=0;
    int nBackFaces=0;
    int nFrontFaces=0;

#define  count3(type) ++n##type ##Faces;
#define ucount3(type) count3(type)
#define  count2(type) count3(type)
#define ucount2(type) count3(type)
#define  count1(type) count3(type)
#define ucount1(type) count3(type)


    for(int kk=1;kk<=nz_2;kk++)
     for(int jj=1;jj<=ny_2;jj++)
      for(int ii=1;ii<=nx_2;ii++)
       {
			unsigned char vv=vxlImg(ii,jj,kk);
            if (vv==1 || vv==2)
            {
				if (ii!=1)
				{
					if (vv==vxlImg(ii-1,jj,kk))    {/*count3( Internal);internal[iface]=true;*/}
					else  if(1!=vxlImg(ii-1,jj,kk)) count3( All);
				}//else {count3( Back);}
				if (ii!=nx_2)
				{
					if (vv==vxlImg(ii+1,jj,kk))    {/*ucount3( Internal);internal[iface]=true;*/}
					else if(1!=vxlImg(ii+1,jj,kk)) ucount3( All);
				}//else {ucount3( Front); }


				if (jj!=1)
				{
					if (vv==vxlImg(ii,jj-1,kk))    {/*count2( Internal);internal[iface]=true;*/}
					else if(1!=vxlImg(ii,jj-1,kk)) count2( All);
				}//else {count2( Bottom)}
				if (jj!=ny_2)
				{
					if (vv==vxlImg(ii,jj+1,kk))    {/*ucount2( Internal);internal[iface]=true;*/}
					else if(1!=vxlImg(ii,jj+1,kk)) ucount2( All);
				}//else {ucount2( Top)}


				if (kk!=1)
				{
					if (vv==vxlImg(ii,jj,kk-1))    {/*count1( Internal);internal[iface]=true;*/}
					else  if(1!=vxlImg(ii,jj,kk-1)) count1( All);
				}//else {count1( Left)}
				if (kk!=nz_2)
				{
					if (vv==vxlImg(ii,jj,kk+1))    {/*ucount1( Internal);internal[iface]=true;*/}
					else  if(1!=vxlImg(ii,jj,kk+1)) ucount1( All);
				}//else {ucount1( Right)}
            }
        }

    Info<<"nInternalFaces: "<<nInternalFaces<<endl;
    Info<<"nAllFaces: "<<nAllFaces<<endl;
    Info<<"nLeftFaces: "<<nLeftFaces<<endl;
    Info<<"nRightFaces: "<<nRightFaces<<endl;
    Info<<"nBottomFaces: "<<nBottomFaces<<endl;
    Info<<"nTopFaces: "<<nTopFaces<<endl;
    Info<<"nBackFaces: "<<nBackFaces<<endl;
    Info<<"nFrontFaces: "<<nFrontFaces<<endl;

    std::cout.flush();



     int iPoints=-1;



	voxelField<int> point_mapper(n[0]+1,n[1]+1,n[2]+1,-1);

	DynamicField<point> points;

	faceList faces_All(nAllFaces);
	labelList fMarks(nAllFaces);
	faceList faces_Left(nLeftFaces);
	faceList faces_Right(nRightFaces );
	faceList faces_Bottom(nBottomFaces);
	faceList faces_Top(nTopFaces);
	faceList faces_Back(nBackFaces);
	faceList faces_Front(nFrontFaces);



	#define addPointToFace_m(pointIndex, kkk,jjj,iii, type)                    \
		if (point_mapper(iii,jjj,kkk)<0) {                                        \
			points.append(point(dx.x*(iii)+X1.x, dx.y*(jjj)+X1.y, dx.z*(kkk)+X1.z)); \
			point_mapper(iii,jjj,kkk)=++iPoints;                                      \
		}                                                                            \
		faces_##type[i##type ##Faces][pointIndex]=point_mapper(iii,jjj,kkk);


	#define recordFaces_m( l10,l11,l20,l21,l30,l31, ii,jj,kk,type )                      \
	  {                                                                                  \
		      ++i##type##Faces ;                                                        \
		      faces_##type[i##type ##Faces].setSize(4);                                  \
		      addPointToFace_m(0, ii, jj, kk, type)                                      \
		      addPointToFace_m(1, ii+l10, jj+l20, kk+l30, type)                          \
		      addPointToFace_m(2, ii+l10+l11, jj+l20+l21, kk+l30+l31, type)              \
		      addPointToFace_m(3, ii+l11, jj+l21, kk+l31, type)                          \
	  }


	#define  clockwiserecordFaces3(type) recordFaces_m( 0,1,1,0,0,0,  kk-1,jj-1,ii-1, type)
	#define uclockwiserecordFaces3(type) recordFaces_m( 1,0,0,1,0,0,  kk-1,jj-1,ii  , type)
	#define  clockwiserecordFaces2(type) recordFaces_m( 1,0,0,0,0,1,  kk-1,jj-1,ii-1, type)
	#define uclockwiserecordFaces2(type) recordFaces_m( 0,1,0,0,1,0,  kk-1,jj  ,ii-1, type)
	#define  clockwiserecordFaces1(type) recordFaces_m( 0,0,0,1,1,0,  kk-1,jj-1,ii-1, type)
	#define uclockwiserecordFaces1(type) recordFaces_m( 0,0,1,0,0,1,  kk  ,jj-1,ii-1, type)








    int iAllFaces=-1;

    for (int kk=1;kk<=nz_2;kk++)
      for (int jj=1;jj<=ny_2;jj++)
        for (int ii=1;ii<=nx_2;ii++)
        {
			unsigned char vv=vxlImg(ii,jj,kk);
            if (vv==1 || vv==2)
            {

				if (ii!=1)
				{
					if (vv==vxlImg(ii-1,jj,kk))    {/*clockwiserecordFaces3( Internal);internal[iface]=true;*/}
					else  if(1!=vxlImg(ii-1,jj,kk)) {clockwiserecordFaces3( All); fMarks[iAllFaces]=vv+vxlImg(ii-1,jj,kk); }
				}//else {clockwiserecordFaces3( Back);}
				if (ii!=nx_2)
				{
					if (vv==vxlImg(ii+1,jj,kk))    {/*uclockwiserecordFaces3( Internal);internal[iface]=true;*/}
					else if(1!=vxlImg(ii+1,jj,kk)){uclockwiserecordFaces3( All); fMarks[iAllFaces]=vv+vxlImg(ii+1,jj,kk); }
				}//else {uclockwiserecordFaces3( Front); }


				if (jj!=1)
				{
					if (vv==vxlImg(ii,jj-1,kk))    {/*clockwiserecordFaces2( Internal);internal[iface]=true;*/}
					else if(1!=vxlImg(ii,jj-1,kk)) {clockwiserecordFaces2( All);  fMarks[iAllFaces]=vv+vxlImg(ii,jj-1,kk); }
				}//else {clockwiserecordFaces2( Bottom)}
				if (jj!=ny_2)
				{
					if (vv==vxlImg(ii,jj+1,kk))    {/*uclockwiserecordFaces2( Internal);internal[iface]=true;*/}
					else if(1!=vxlImg(ii,jj+1,kk)) {uclockwiserecordFaces2( All);  fMarks[iAllFaces]=vv+vxlImg(ii,jj+1,kk); }
				}//else {uclockwiserecordFaces2( Top)}


				if (kk!=1)
				{
					if (vv==vxlImg(ii,jj,kk-1))    {/*clockwiserecordFaces1( Internal);internal[iface]=true;*/}
					else  if(1!=vxlImg(ii,jj,kk-1)) {clockwiserecordFaces1( All); fMarks[iAllFaces]=vv+vxlImg(ii,jj,kk-1);}
				}//else {clockwiserecordFaces1( Left)}
				if (kk!=nz_2)
				{
					if (vv==vxlImg(ii,jj,kk+1))    {/*uclockwiserecordFaces1( Internal);internal[iface]=true;*/}
					else  if(1!=vxlImg(ii,jj,kk+1)) {uclockwiserecordFaces1( All); fMarks[iAllFaces]=vv+vxlImg(ii,jj,kk+1);}
				}//else {uclockwiserecordFaces1( Right)}

            }
        }



//______________________________________________

	point_mapper.reset({0,0,0});


    Info<<"nPoints: "<<points.size()<<"      "<<endl;        /*Info.flush()*/;


	correct(faces_All,fMarks,points,true);
	correct(faces_All,fMarks,points,true);
	correct(faces_All,fMarks,points,true);
	correct(faces_All,fMarks,points,true);
	correct(faces_All,fMarks,points,false);
	correct(faces_All,fMarks,points,false);
	correct(faces_All,fMarks,points,false);
	correct(faces_All,fMarks,points,false);


Info<<1<<endl;
	correctbioti(faces_All,fMarks,points,0);
	correctbioti(faces_All,fMarks,points,0);
	correctbioti(faces_All,fMarks,points,0);
	correctbioti(faces_All,fMarks,points,0);
	correctbioti(faces_All,fMarks,points,0);
	correctbioti(faces_All,fMarks,points,0);
	correctbioti(faces_All,fMarks,points,0);
	correctbioti(faces_All,fMarks,points,1);
	correctbioti(faces_All,fMarks,points,1);
	correctbioti(faces_All,fMarks,points,0);
	correctbioti(faces_All,fMarks,points,1);
	correctbioti(faces_All,fMarks,points,0);
	correctbioti(faces_All,fMarks,points,1);
	correctbioti(faces_All,fMarks,points,1);
	correctbioti(faces_All,fMarks,points,1);
	correctbioti(faces_All,fMarks,points,1);
	correctbioti(faces_All,fMarks,points,1);

Info<<2<<endl;


Info<<endl;




		faceList facesSorted_All(faces_All.size());
		labelList zoneSizes(3);

		int indF=-1;
		for (int i=0; i<faces_All.size(); ++i)
			if(fMarks[i]==1) facesSorted_All[++indF]=faces_All[i];
		zoneSizes[0]=indF+1;
		for (int i=0; i<faces_All.size(); ++i)
			if(fMarks[i]==2) facesSorted_All[++indF]=faces_All[i];
		zoneSizes[1]=indF+1-zoneSizes[0];
		for (int i=0; i<faces_All.size(); ++i)
			if(fMarks[i]==3) facesSorted_All[++indF]=faces_All[i];
			else if(fMarks[i]>3) Info<<" Error wrong  label in the image"<<endl;
		zoneSizes[2]=indF+1-zoneSizes[1]-zoneSizes[0];

		Info<<faces_All.size()<<" ?= "<<facesSorted_All.size()<<endl;
		Info<<zoneSizes[0]<<" + "<<zoneSizes[1]<<" + "<<zoneSizes[2]<<" = "<<indF+1<<endl;

    {
        Xfer<List<face> > facesFer(facesSorted_All,true);
        Field<point> & pointsSF=points;
        Xfer<Field<point> > pointsFer(pointsSF,true);
        meshedSurface surf1(pointsFer, facesFer, zoneSizes);
		//~ surf1.addZones();
        //~ surf.triangulate ();

        surf1.write(outputSurface);


        //std::ofstream  offile("fMarks.txt");
        //for (int i=0; i<fMarks.size(); ++i)		offile<<fMarks[i]<<"  ";
		//offile.close();
    }

}
