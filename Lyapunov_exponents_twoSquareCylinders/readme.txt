Run Directory: Lyapunov Exponents of Chaotic Flow over Two Square Cylinders

This folder contains the core setup and execution logic for computing the Lyapunov exponents (LEs) of the chaotic flow past two side-by-side square cylinders.
🧰 Overview

Two solvers are used in sequence:

1. SC_Evolver – Nonlinear Flow + Linearised Perturbations

A custom parallelisable solver, modified from pimpleFoam, used to:

    Evolve the nonlinear Navier–Stokes equations in time

    Simultaneously evolve perturbations using the Linearised Navier-Stokes equations

Purpose: To linearise the flow over its unsteady system trajectory.

🔁 2. QR_code – Orthonormalisation via Gram-Schmidt-based QR Decomposition

A serial solver that:

    Applies the Gram–Schmidt process (QR decomposition) on the tangent vectors

    Is executed every t_step = 0.3 time units

    Resets the perturbations using the new orthonormalised basis (Q)

    Stores the R matrices for post-processing of Lyapunov exponents

Note:

QR_code is not parallelisable.

📄 Output Files

NonL.txt: Snapshot of the nonlinear flow solution at each orthonormalisation point
U_D_collection.txt: Contains the six Gram–Schmidt vectors (Q matrix columns)
R_full.txt: Full upper-triangular R matrices from QR decomposition, saved at each step

📈 Post-Processing: LE Calculation

Once the simulation has evolved for a sufficient time (recommended: T > 2000 time units), you can compute the Lyapunov exponents by running LE_calc.m

This MATLAB script processes the R_full.txt file to extract the spectrum of Lyapunov exponents using the Benettin's LE algorithm;

✅ Suggested Workflow

    Run SC_Evolver in parallel to evolve the nonlinear flow + perturbations (Gram-Schmidt vectors)

    At each t = 0.3, stop and run QR_code to orthonormalise perturbations and save R

    Repeat until desired final time

    Run LE_calc.m in MATLAB to compute the Lyapunov exponents
    
📄 Bash Script 

SC_run.sh is a sample bash script that performs the above workflow on a High-Performance-Computing cluster by decomposing the domain over 32 processors.     
 
