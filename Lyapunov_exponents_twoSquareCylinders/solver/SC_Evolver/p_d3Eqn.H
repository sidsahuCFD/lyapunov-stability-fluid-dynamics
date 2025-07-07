volScalarField rAU_d3(1.0/U_d3Eqn.A());
volVectorField HbyA_d3(constrainHbyA(rAU_d3*U_d3Eqn.H(), U_d3, p_d3));
surfaceScalarField phiHbyA_d3
(
    "phiHbyA_d3",
    fvc::flux(HbyA_d3)
  + MRF.zeroFilter(fvc::interpolate(rAU_d3)*fvc::ddtCorr(U_d3, phi_d3, U_d3f))
);

MRF.makeRelative(phiHbyA_d3);

if (p_d3.needReference())
{
    fvc::makeRelative(phiHbyA_d3, U_d3);
    adjustPhi(phiHbyA_d3, U_d3, p_d3);
    fvc::makeAbsolute(phiHbyA_d3, U_d3);
}

tmp<volScalarField> rAtU_d3(rAU_d3);

if (pimple.consistent())
{
    rAtU_d3 = 1.0/max(1.0/rAU_d3 - U_d3Eqn.H1(), 0.1/rAU_d3);
    phiHbyA_d3 +=
        fvc::interpolate(rAtU_d3() - rAU_d3)*fvc::snGrad(p_d3)*mesh.magSf();
    HbyA_d3 -= (rAU_d3 - rAtU_d3())*fvc::grad(p_d3);
}

if (pimple.nCorrPiso() <= 1)
{
    tU_d3Eqn.clear();
}

// Update the pressure BCs to ensure flux consistency
constrainPressure(p_d3, U_d3, phiHbyA_d3, rAtU_d3(), MRF);

// Non-orthogonal pressure corrector loop
while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix p_d3Eqn
    (
        fvm::laplacian(rAtU_d3(), p_d3) == fvc::div(phiHbyA_d3)
    );

    p_d3Eqn.setReference
    (
        pressureReference.refCell(), // Not fighting this, since I'm not using the non-orthogonal correction
        pressureReference.refValue()
    );

    p_d3Eqn.solve();

    if (pimple.finalNonOrthogonalIter())
    {
        phi_d3 = phiHbyA_d3 - p_d3Eqn.flux();
    }
}

#include "deltaContinuityErrs.H"

// Explicitly relax pressure for momentum corrector
p_d3.relax();

U_d3 = HbyA_d3 - rAtU_d3*fvc::grad(p_d3);
U_d3.correctBoundaryConditions();
fvConstraints.constrain(U_d3);

// Correct Uf if the mesh is moving
fvc::correctUf(U_d3f, U_d3, phi_d3, MRF);

// Make the fluxes relative to the mesh motion
fvc::makeRelative(phi_d3, U_d3);
