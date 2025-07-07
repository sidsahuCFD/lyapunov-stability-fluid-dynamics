// Solve the Momentum equation

MRF.correctBoundaryVelocity(U_d3);

tmp<fvVectorMatrix> tU_d3Eqn
(
    fvm::ddt(U_d3) + fvc::div(phi_d3,U) + fvm::div(phi, U_d3)
  + MRF.DDt(U_d3)
  + (turbulence_d3->divDevSigma(U_d3))/Re
);

fvVectorMatrix& U_d3Eqn = tU_d3Eqn.ref();

U_d3Eqn.relax();

fvConstraints.constrain(U_d3Eqn);

if (pimple.momentumPredictor())
{
    solve(U_d3Eqn == -fvc::grad(p_d3));

    fvConstraints.constrain(U_d3);
}
