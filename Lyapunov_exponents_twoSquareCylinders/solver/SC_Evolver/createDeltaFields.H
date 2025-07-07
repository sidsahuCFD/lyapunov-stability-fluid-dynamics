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
   mesh//,
   //dimensionedScalar("p_d", dimensionSet(0,2,-2,0,0,0,0), Foam::scalar(0))
);	



//////////////////////////////////////////////////////

volVectorField U_d2
(
   IOobject
   (
    "U_d2",
    "0",   
    mesh,
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE
   ),
   mesh//,
   //dimensionedVector("U_d2", dimensionSet(0,1,-1,0,0,0,0), Foam::vector(0,0,0))
);		
	
Info<< "Reading field p_d2\n" << endl;
volScalarField p_d2
(
    IOobject
    (
        "p_d2",
    	"0",   
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
   mesh//,
   //dimensionedScalar("p_d2", dimensionSet(0,2,-2,0,0,0,0), Foam::scalar(0))
);	


////////////////////////////////////////////////////////

volVectorField U_d3
(
   IOobject
   (
    "U_d3",
    "0",   
    mesh,
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE
   ),
   mesh//,
   //dimensionedVector("U_d3", dimensionSet(0,1,-1,0,0,0,0), Foam::vector(0,0,0))
);		
	
Info<< "Reading field p_d3\n" << endl;
volScalarField p_d3
(
    IOobject
    (
        "p_d3",
    	"0",   
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
   mesh//,
   //dimensionedScalar("p_d3", dimensionSet(0,2,-2,0,0,0,0), Foam::scalar(0))
);	


////////////////////////////////////////////////////////

volVectorField U_d4
(
   IOobject
   (
    "U_d4",
    "0",
    mesh,
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE
   ),
   mesh//,
   //dimensionedVector("U_d4", dimensionSet(0,1,-1,0,0,0,0), Foam::vector(0,0,0))
);		
	
Info<< "Reading field p_d4\n" << endl;
volScalarField p_d4
(
    IOobject
    (
        "p_d4",
    	"0",        
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
   mesh//,
   //dimensionedScalar("p_d4", dimensionSet(0,2,-2,0,0,0,0), Foam::scalar(0))
);	


////////////////////////////////////////////////////////

volVectorField U_d5
(
   IOobject
   (
    "U_d5",
    "0",
    mesh,
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE
   ),
   mesh
);		
	
Info<< "Reading field p_d5\n" << endl;
volScalarField p_d5
(
    IOobject
    (
        "p_d5",
    	"0",        
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
   mesh
);	


////////////////////////////////////////////////////////

volVectorField U_d6
(
   IOobject
   (
    "U_d6",
    "0",
    mesh,
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE
   ),
   mesh
);		
	
Info<< "Reading field p_d6\n" << endl;
volScalarField p_d6
(
    IOobject
    (
        "p_d6",
    	"0",        
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
   mesh
);	



