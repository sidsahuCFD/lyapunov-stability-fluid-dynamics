volScalarField rAU_d4(1.0/U_d4Eqn.A());
volVectorField HbyA_d4(constrainHbyA(rAU_d4*U_d4Eqn.H(), U_d4, p_d4));
surfaceScalarField phiHbyA_d4
(
    "phiHbyA_d4",
    fvc::flux(HbyA_d4)
  + MRF.zeroFilter(fvc::interpolate(rAU_d4)*fvc::ddtCorr(U_d4, phi_d4, U_d4f))
);

MRF.makeRelative(phiHbyA_d4);

if (p_d4.needReference())
{
    fvc::makeRelative(phiHbyA_d4, U_d4);
    adjustPhi(phiHbyA_d4, U_d4, p_d4);
    fvc::makeAbsolute(phiHbyA_d4, U_d4);
}

tmp<volScalarField> rAtU_d4(rAU_d4);

if (pimple.consistent())
{
    rAtU_d4 = 1.0/max(1.0/rAU_d4 - U_d4Eqn.H1(), 0.1/rAU_d4);
    phiHbyA_d4 +=
        fvc::interpolate(rAtU_d4() - rAU_d4)*fvc::snGrad(p_d4)*mesh.magSf();
    HbyA_d4 -= (rAU_d4 - rAtU_d4())*fvc::grad(p_d4);
}

if (pimple.nCorrPiso() <= 1)
{
    tU_d4Eqn.clear();
}

// Update the pressure BCs to ensure flux consistency
constrainPressure(p_d4, U_d4, phiHbyA_d4, rAtU_d4(), MRF);

// Non-orthogonal pressure corrector loop
while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix p_d4Eqn
    (
        fvm::laplacian(rAtU_d4(), p_d4) == fvc::div(phiHbyA_d4)
    );

    p_d4Eqn.setReference
    (
        pressureReference.refCell(), // Not fighting this, since I'm not using the non-orthogonal correction
        pressureReference.refValue()
    );

    p_d4Eqn.solve();

    if (pimple.finalNonOrthogonalIter())
    {
        phi_d4 = phiHbyA_d4 - p_d4Eqn.flux();
    }
}

#include "deltaContinuityErrs.H"

// Explicitly relax pressure for momentum corrector
p_d4.relax();

U_d4 = HbyA_d4 - rAtU_d4*fvc::grad(p_d4);
U_d4.correctBoundaryConditions();
fvConstraints.constrain(U_d4);

// Correct Uf if the mesh is moving
fvc::correctUf(U_d4f, U_d4, phi_d4, MRF);

// Make the fluxes relative to the mesh motion
fvc::makeRelative(phi_d4, U_d4);
