// Solve the Momentum equation

MRF.correctBoundaryVelocity(U_d2);

tmp<fvVectorMatrix> tU_d2Eqn
(
    fvm::ddt(U_d2) + fvc::div(phi_d2,U) + fvm::div(phi, U_d2)
  + MRF.DDt(U_d2)
  + (turbulence_d2->divDevSigma(U_d2))/Re
);

fvVectorMatrix& U_d2Eqn = tU_d2Eqn.ref();

U_d2Eqn.relax();

fvConstraints.constrain(U_d2Eqn);

if (pimple.momentumPredictor())
{
    solve(U_d2Eqn == -fvc::grad(p_d2));

    fvConstraints.constrain(U_d2);
}
