volVectorField U_d
(
   IOobject
   (
    "U_d",
    "0",
    mesh,
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE
   ),
   mesh//,
   //dimensionedVector("U_d", dimensionSet(0,1,-1,0,0,0,0), Foam::vector(0,0,0))
);
		
Info<< "Reading field p_d\n" << endl;
volScalarField p_d
(
    IOobject
    (
        "p_d",
        "0",
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
   mesh
   //dimensionedScalar("p_d", dimensionSet(0,2,-2,0,0,0,0), Foam::scalar(0))
);	

U_vec = U_D.col(i);

if (Pstream::master()){
exam << "In createDelftaFields, after accepting the U_d from U_D the norm of U_d is" << U_vec.norm() <<endl;
}

k = 0;				
forAll(U_d, id)
{	
    U_d[id].component(0) = U_vec(k);k++; 
    U_d[id].component(1)= U_vec(k);k++; 	 	  
}

   		



