
surfaceScalarField phi_d
(
    IOobject
    (
        "phi_d",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    fvc::flux(U_d)
);
mesh.schemes().setFluxRequired(p_d.name()); 
autoPtr<incompressible::momentumTransportModel> turbulence_d 
(
    incompressible::momentumTransportModel::New(U_d, phi_d, viscosity)
);


//////////////////////////////////////////////////////


surfaceScalarField phi_d2
(
    IOobject
    (
        "phi_d2",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    fvc::flux(U_d2)
);

mesh.schemes().setFluxRequired(p_d2.name()); 
autoPtr<incompressible::momentumTransportModel> turbulence_d2 
(
    incompressible::momentumTransportModel::New(U_d2, phi_d2, viscosity)
);

////////////////////////////////////////////////////////


surfaceScalarField phi_d3
(
    IOobject
    (
        "phi_d3",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    fvc::flux(U_d3)
);
mesh.schemes().setFluxRequired(p_d3.name()); 
autoPtr<incompressible::momentumTransportModel> turbulence_d3 
(
    incompressible::momentumTransportModel::New(U_d3, phi_d3, viscosity)
);

////////////////////////////////////////////////////////

surfaceScalarField phi_d4
(
    IOobject
    (
        "phi_d4",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    fvc::flux(U_d4)
);
mesh.schemes().setFluxRequired(p_d4.name()); 
autoPtr<incompressible::momentumTransportModel> turbulence_d4 
(
    incompressible::momentumTransportModel::New(U_d4, phi_d4, viscosity)
);


////////////////////////////////////////////////////////

surfaceScalarField phi_d5
(
    IOobject
    (
        "phi_d5",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    fvc::flux(U_d5)
);
mesh.schemes().setFluxRequired(p_d5.name()); 
autoPtr<incompressible::momentumTransportModel> turbulence_d5 
(
    incompressible::momentumTransportModel::New(U_d5, phi_d5, viscosity)
);

////////////////////////////////////////////////////////

surfaceScalarField phi_d6
(
    IOobject
    (
        "phi_d6",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    fvc::flux(U_d6)
);
mesh.schemes().setFluxRequired(p_d6.name()); 
autoPtr<incompressible::momentumTransportModel> turbulence_d6 
(
    incompressible::momentumTransportModel::New(U_d6, phi_d6, viscosity)
);
