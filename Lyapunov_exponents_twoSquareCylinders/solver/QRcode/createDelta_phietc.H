
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





