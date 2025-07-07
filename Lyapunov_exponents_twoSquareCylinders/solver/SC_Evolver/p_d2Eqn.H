volScalarField rAU_d2(1.0/U_d2Eqn.A());
volVectorField HbyA_d2(constrainHbyA(rAU_d2*U_d2Eqn.H(), U_d2, p_d2));
surfaceScalarField phiHbyA_d2
(
    "phiHbyA_d2",
    fvc::flux(HbyA_d2)
  + MRF.zeroFilter(fvc::interpolate(rAU_d2)*fvc::ddtCorr(U_d2, phi_d2, U_d2f))
);

MRF.makeRelative(phiHbyA_d2);

if (p_d2.needReference())
{
    fvc::makeRelative(phiHbyA_d2, U_d2);
    adjustPhi(phiHbyA_d2, U_d2, p_d2);
    fvc::makeAbsolute(phiHbyA_d2, U_d2);
}

tmp<volScalarField> rAtU_d2(rAU_d2);

if (pimple.consistent())
{
    rAtU_d2 = 1.0/max(1.0/rAU_d2 - U_d2Eqn.H1(), 0.1/rAU_d2);
    phiHbyA_d2 +=
        fvc::interpolate(rAtU_d2() - rAU_d2)*fvc::snGrad(p_d2)*mesh.magSf();
    HbyA_d2 -= (rAU_d2 - rAtU_d2())*fvc::grad(p_d2);
}

if (pimple.nCorrPiso() <= 1)
{
    tU_d2Eqn.clear();
}

// Update the pressure BCs to ensure flux consistency
constrainPressure(p_d2, U_d2, phiHbyA_d2, rAtU_d2(), MRF);

// Non-orthogonal pressure corrector loop
while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix p_d2Eqn
    (
        fvm::laplacian(rAtU_d2(), p_d2) == fvc::div(phiHbyA_d2)
    );

    p_d2Eqn.setReference
    (
        pressureReference.refCell(), // Not fighting this, since I'm not using the non-orthogonal correction
        pressureReference.refValue()
    );

    p_d2Eqn.solve();

    if (pimple.finalNonOrthogonalIter())
    {
        phi_d2 = phiHbyA_d2 - p_d2Eqn.flux();
    }
}

#include "deltaContinuityErrs.H"

// Explicitly relax pressure for momentum corrector
p_d2.relax();

U_d2 = HbyA_d2 - rAtU_d2*fvc::grad(p_d2);
U_d2.correctBoundaryConditions();
fvConstraints.constrain(U_d2);

// Correct Uf if the mesh is moving
fvc::correctUf(U_d2f, U_d2, phi_d2, MRF);

// Make the fluxes relative to the mesh motion
fvc::makeRelative(phi_d2, U_d2);
