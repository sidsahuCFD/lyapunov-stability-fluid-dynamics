// Calculate absolute flux
// from the mapped surface velocity

phi_d = mesh.Sf() & U_df();

correctUphiBCs(U_d, phi_d, true);

CorrectPhi
(
    phi_d,
    U_d,
    p_d,
    dimensionedScalar("rAUf", dimTime, 1),
    geometricZeroField(),
    pressureReference,
    pimple
);

//#include "continuityErrs_d.H"

// Make the flux relative to the mesh motion
fvc::makeRelative(phi_d, U_d);

// 09_02 This is not necessary
