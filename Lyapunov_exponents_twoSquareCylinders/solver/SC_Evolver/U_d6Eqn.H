// Solve the Momentum equation

MRF.correctBoundaryVelocity(U_d6);

tmp<fvVectorMatrix> tU_d6Eqn
(
    fvm::ddt(U_d6) + fvc::div(phi_d6,U) + fvm::div(phi, U_d6)
  + MRF.DDt(U_d6)
  + (turbulence_d6->divDevSigma(U_d6))/Re
);

fvVectorMatrix& U_d6Eqn = tU_d6Eqn.ref();

U_d6Eqn.relax();

fvConstraints.constrain(U_d6Eqn);

if (pimple.momentumPredictor())
{
    solve(U_d6Eqn == -fvc::grad(p_d6));

    fvConstraints.constrain(U_d6);
}
