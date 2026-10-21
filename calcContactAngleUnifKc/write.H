





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
        if ((pI+1)%1==0)    offRc<<"";
        offRc<<pKc[pI]<<" "<<"  \n";
    }
}
offRc<<"Kcw 1 "<<pNeips.size()<<" float"<<std::endl;
forAll(pNeips, pI)
{
    if(iter==nIters-1)
    {
        if ((pI+1)%1==0)    offRc<<"";
        offRc<<pKcw[pI]<<" "<<"  \n";
    }
}
offRc<<"\n\n RcCL 1 "<<pNeips.size()<<" float"<<std::endl;
forAll(pNeips, pI)
{
    if(iter==nIters-1)
    {
        if ((pI+1)%1==0)    offRc<<"";
        offRc<<pKcCL[pI]<<" "<<"  \n";
    }
}
offRc<<"\n\n InterfNorms 3 "<<pNeips.size()<<" float"<<std::endl;
forAll(pNeips, pI)
{
    if(iter==nIters-1)
    {
        if ((pI+1)%1==0)    offRc<<"";
        offRc<<pNs[pI][0]<<" "<<pNs[pI][1]<<" "<<pNs[pI][2]<<" "<<" \n";
    }
}
offRc<<"\n\n wallNorms 3 "<<pNeips.size()<<" float"<<std::endl;
forAll(pNeips, pI)
{
    if(iter==nIters-1)
    {
        if ((pI+1)%1==0)    offRc<<"";
        offRc<<pNw[pI][0]<<" "<<pNw[pI][1]<<" "<<pNw[pI][2]<<" "<<"  \n";
    }
}

offRc.close();

//newPoints += displacementsRc;
