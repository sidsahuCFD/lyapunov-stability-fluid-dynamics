        // Momentum predictor
  
        fvVectorMatrix U_dEqn
        (
            fvm::ddt(U_d)
          + fvc::div(phi_d,U)
          + fvm::div(phi,U_d)         
          - fvm::laplacian(nu,U_d)/Re
        );
        
        if (piso.momentumPredictor())
        {
            solve(U_dEqn == -fvc::grad(p_d));
        }

        // --- PISO loop
        while (piso.correct())
        {
            volScalarField rAU_d(1.0/U_dEqn.A());
            volVectorField HbyA_d(constrainHbyA(rAU_d*U_dEqn.H(), U_d, p_d));
            surfaceScalarField phiHbyA_d
            (
                "phiHbyA_d",
                fvc::flux(HbyA_d)
              + fvc::interpolate(rAU_d)*fvc::ddtCorr(U_d, phi_d)
            );

            adjustPhi(phiHbyA_d, U_d, p_d);

            // Update the pressure BCs to ensure flux consistency
            constrainPressure(p_d, U_d, phiHbyA_d, rAU_d);

            // Non-orthogonal pressure corrector loop
            while (piso.correctNonOrthogonal())
            {
                // Pressure corrector

                fvScalarMatrix pEqn_d
                (
                    fvm::laplacian(rAU_d, p_d) == fvc::div(phiHbyA_d)
                );

                pEqn_d.setReference(p_dRefCell, p_dRefValue);

                pEqn_d.solve();

                if (piso.finalNonOrthogonalIter())
                {
                    phi_d = phiHbyA_d - pEqn_d.flux();
                }
            }

            #include "deltaContinuityErrs.H"

            U_d = HbyA_d - rAU_d*fvc::grad(p_d);
            U_d.correctBoundaryConditions();
        }
