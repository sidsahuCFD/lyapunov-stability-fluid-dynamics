volScalarField rAU_d(1.0/U_dEqn.A());
volVectorField HbyA_d(constrainHbyA(rAU_d*U_dEqn.H(), U_d, p_d));
surfaceScalarField phiHbyA_d
(
    "phiHbyA_d",
    fvc::flux(HbyA_d)
  + MRF.zeroFilter(fvc::interpolate(rAU_d)*fvc::ddtCorr(U_d, phi_d, U_df))
);

MRF.makeRelative(phiHbyA_d);

if (p_d.needReference())
{
    fvc::makeRelative(phiHbyA_d, U_d);
    adjustPhi(phiHbyA_d, U_d, p_d);
    fvc::makeAbsolute(phiHbyA_d, U_d);
}

tmp<volScalarField> rAtU_d(rAU_d);

if (pimple.consistent())
{
    rAtU_d = 1.0/max(1.0/rAU_d - U_dEqn.H1(), 0.1/rAU_d);
    phiHbyA_d +=
        fvc::interpolate(rAtU_d() - rAU_d)*fvc::snGrad(p_d)*mesh.magSf();
    HbyA_d -= (rAU_d - rAtU_d())*fvc::grad(p_d);
}

if (pimple.nCorrPiso() <= 1)
{
    tU_dEqn.clear();
}

// Update the pressure BCs to ensure flux consistency
constrainPressure(p_d, U_d, phiHbyA_d, rAtU_d(), MRF);

// Non-orthogonal pressure corrector loop
while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix p_dEqn
    (
        fvm::laplacian(rAtU_d(), p_d) == fvc::div(phiHbyA_d)
    );

    p_dEqn.setReference
    (
        pressureReference.refCell(), // Not fighting this, since I'm not using the non-orthogonal correction
        pressureReference.refValue()
    );

    p_dEqn.solve();

    if (pimple.finalNonOrthogonalIter())
    {
        phi_d = phiHbyA_d - p_dEqn.flux();
    }
}

#include "deltaContinuityErrs.H"

// Explicitly relax pressure for momentum corrector
p_d.relax();

U_d = HbyA_d - rAtU_d*fvc::grad(p_d);
U_d.correctBoundaryConditions();
fvConstraints.constrain(U_d);

// Correct Uf if the mesh is moving
fvc::correctUf(U_df, U_d, phi_d, MRF);

// Make the fluxes relative to the mesh motion
fvc::makeRelative(phi_d, U_d);
