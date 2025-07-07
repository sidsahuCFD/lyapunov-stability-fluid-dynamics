 volVectorField Uvec2Uvvf(Eigen::VectorXd U_vec)
 {

   #include "createMesh.H"
   
   volVectorField U_vvf
   (  
   IOobject
   (
    "U_vvf",
    0,
    mesh,
    IOobject::NO_READ,
    IOobject::AUTO_WRITE
   ),
   mesh,
   dimensionedVector("U_vvf", dimensionSet(0,1,-1,0,0,0,0), Foam::vector(0,0,0))
   );
   
   int k(0);
   forAll(U_vvf, id)
   {	
     U_vvf[id].component(0) = U_vec(k);k++;
     U_vvf[id].component(1) = U_vec(k);k++;
     U_vvf[id].component(2) = 0;
   }
   return(U_vvf);
}
