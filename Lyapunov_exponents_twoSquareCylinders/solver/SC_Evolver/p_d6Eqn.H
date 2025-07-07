volScalarField rAU_d6(1.0/U_d6Eqn.A());
volVectorField HbyA_d6(constrainHbyA(rAU_d6*U_d6Eqn.H(), U_d6, p_d6));
surfaceScalarField phiHbyA_d6
(
    "phiHbyA_d6",
    fvc::flux(HbyA_d6)
  + MRF.zeroFilter(fvc::interpolate(rAU_d6)*fvc::ddtCorr(U_d6, phi_d6, U_d6f))
);

MRF.makeRelative(phiHbyA_d6);

if (p_d6.needReference())
{
    fvc::makeRelative(phiHbyA_d6, U_d6);
    adjustPhi(phiHbyA_d6, U_d6, p_d6);
    fvc::makeAbsolute(phiHbyA_d6, U_d6);
}

tmp<volScalarField> rAtU_d6(rAU_d6);

if (pimple.consistent())
{
    rAtU_d6 = 1.0/max(1.0/rAU_d6 - U_d6Eqn.H1(), 0.1/rAU_d6);
    phiHbyA_d6 +=
        fvc::interpolate(rAtU_d6() - rAU_d6)*fvc::snGrad(p_d6)*mesh.magSf();
    HbyA_d6 -= (rAU_d6 - rAtU_d6())*fvc::grad(p_d6);
}

if (pimple.nCorrPiso() <= 1)
{
    tU_d6Eqn.clear();
}

// Update the pressure BCs to ensure flux consistency
constrainPressure(p_d6, U_d6, phiHbyA_d6, rAtU_d6(), MRF);

// Non-orthogonal pressure corrector loop
while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix p_d6Eqn
    (
        fvm::laplacian(rAtU_d6(), p_d6) == fvc::div(phiHbyA_d6)
    );

    p_d6Eqn.setReference
    (
        pressureReference.refCell(), // Not fighting this, since I'm not using the non-orthogonal correction
        pressureReference.refValue()
    );

    p_d6Eqn.solve();

    if (pimple.finalNonOrthogonalIter())
    {
        phi_d6 = phiHbyA_d6 - p_d6Eqn.flux();
    }
}

#include "deltaContinuityErrs.H"

// Explicitly relax pressure for momentum corrector
p_d6.relax();

U_d6 = HbyA_d6 - rAtU_d6*fvc::grad(p_d6);
U_d6.correctBoundaryConditions();
fvConstraints.constrain(U_d6);

// Correct Uf if the mesh is moving
fvc::correctUf(U_d6f, U_d6, phi_d6, MRF);

// Make the fluxes relative to the mesh motion
fvc::makeRelative(phi_d6, U_d6);
