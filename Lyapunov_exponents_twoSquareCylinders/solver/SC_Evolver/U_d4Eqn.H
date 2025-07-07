// Solve the Momentum equation

MRF.correctBoundaryVelocity(U_d4);

tmp<fvVectorMatrix> tU_d4Eqn
(
    fvm::ddt(U_d4) + fvc::div(phi_d4,U) + fvm::div(phi, U_d4)
  + MRF.DDt(U_d4)
  + (turbulence_d4->divDevSigma(U_d4))/Re
);

fvVectorMatrix& U_d4Eqn = tU_d4Eqn.ref();

U_d4Eqn.relax();

fvConstraints.constrain(U_d4Eqn);

if (pimple.momentumPredictor())
{
    solve(U_d4Eqn == -fvc::grad(p_d4));

    fvConstraints.constrain(U_d4);
}
