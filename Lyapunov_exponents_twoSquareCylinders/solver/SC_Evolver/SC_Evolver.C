/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    pimpleFoam to study Lyapunov spectrums

Description
    Transient solver for incompressible, laminar flow of Newtonian fluids.

\*---------------------------------------------------------------------------*/
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "linearalgebra.H"
#include "fvCFD.H"
#include "viscosityModel.H"
#include "OFstream.H"
#include <fstream>
#include "incompressibleMomentumTransportModels.H"
#include "pimpleControl.H"
#include "pressureReference.H"
#include "CorrectPhi.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include <Eigen/Dense>
#include "volVectorField2VectorXd.H"
#include "gramSchmidt.H"
#include "std2Eigen.H"
#include "Eigen2std.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
   
   // Read current time
   scalar T = 0; 
   ifstream myfile;
   myfile.open("currentTime.txt");
   myfile >> T;
   myfile.close();
   
   const int N(84253*2); 
   const int N_p = 6;
   std::ofstream currentTime("currentTime.txt");
   
   #include "postProcess.H"
   #include "setRootCaseLists.H"
   #include "createTime.H"
   #include "createMesh.H"
   #include "initContinuityErrs.H"
   #include "createDyMControls.H"
   #include "createFields.H"
   #include "createUfIfPresent.H"
   #include "createU_dfIfPresent.H" //added
    
   turbulence->validate();
   turbulence_d->validate();
   turbulence_d2->validate();
   turbulence_d3->validate();
   turbulence_d4->validate();
   turbulence_d5->validate(); 
   turbulence_d6->validate();    	    	    	    	    	    	    	    	    	    	    	    	   	
    
    if (!LTS)
    {
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"
    }
        
    while (pimple.run(runTime))
    {
    
	#include "readDyMControls.H"
        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "CourantNo.H"
            #include "setDeltaT.H"
        }
   
        fvModels.preUpdateMesh();

        // Update the mesh for topology change, mesh to mesh mapping
        mesh.update();

        runTime++; 
     
        Info<< "Time = " << runTime.userTimeName() << nl << endl;
	
	#include "primalPIMPLE.H"
  	#include "deltaPIMPLE.H"
  	runTime.write();	    																								
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;	
        
     } // end of time loop
     
    Info<< "End\n" << endl;
    
    if (Pstream::master())
    {
	currentTime << T + runTime.userTimeValue() << endl;
    }
    
    return 0;
}


// ************************************************************************* //
