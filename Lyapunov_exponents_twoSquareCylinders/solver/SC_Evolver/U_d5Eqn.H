// Solve the Momentum equation

MRF.correctBoundaryVelocity(U_d5);

tmp<fvVectorMatrix> tU_d5Eqn
(
    fvm::ddt(U_d5) + fvc::div(phi_d5,U) + fvm::div(phi, U_d5)
  + MRF.DDt(U_d5)
  + (turbulence_d5->divDevSigma(U_d5))/Re
);

fvVectorMatrix& U_d5Eqn = tU_d5Eqn.ref();

U_d5Eqn.relax();

fvConstraints.constrain(U_d5Eqn);

if (pimple.momentumPredictor())
{
    solve(U_d5Eqn == -fvc::grad(p_d5));

    fvConstraints.constrain(U_d5);
}
