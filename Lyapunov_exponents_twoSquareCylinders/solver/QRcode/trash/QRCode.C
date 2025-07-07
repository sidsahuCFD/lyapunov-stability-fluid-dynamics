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
   
   scalar T = 0; 
   ifstream myfile;
   myfile.open("$HOME/OpenFOAM/sidsahu-10/SC/run/currentTime.txt");
   myfile >> T;
   myfile.close();
   
   scalar stept = 10;
   const int N(18483*2); 
   const int N_p = 4;
   
   #include "acceptFields.H"
   
   Eigen::MatrixXd U_D(Eigen::MatrixXd::Random(N,N_p)), U_primal(Eigen::MatrixXd::Random(N,1));   
   
   U_primal.col(0) = volVectorField2VectorXd(U,N);
   U_D.col(0) = volVectorField2VectorXd(U_d,N);
   U_D.col(1) = volVectorField2VectorXd(U_d2,N);
   U_D.col(2) = volVectorField2VectorXd(U_d3,N);
   U_D.col(3) = volVectorField2VectorXd(U_d4,N);
   
   Eigen::VectorXd cum(Eigen::VectorXd::Zero(N_p));
   Eigen::VectorXd lp(Eigen::VectorXd::Zero(N_p));
   
   //Accepting cum array
   ifstream myfile;
   myfile.open("$HOME/OpenFOAM/sidsahu-10/SC/run/cum.txt");
   for (int i = 0; i < N_p; i++) 
   {
       myfile >> cum(i);
   }
   myfile.close();
   
   OFstream U_D_store("$HOME/OpenFOAM/sidsahu-10/SC/run/U_D_collection.txt",std::ios::app);     
   std::ofstream R_full("$HOME/OpenFOAM/sidsahu-10/SC/run/R_full.txt",std::ios::app);     
   std::ofstream NonL("$HOME/OpenFOAM/sidsahu-10/SC/run/NonL.txt",std::ios::app);
   NonL << U_primal << endl;
     
   int i, j, n(N_p), m(N), q_n, r_m, k(0);
   bool full(0);
   if(full==0)
   {
   q_n = n;
   r_m = n;
   }
   else if(full==1)
   {
   q_n = m;
   r_m = m;
   Info <<"This is a full Gram-Schmidt decomposition" <<endl;
   }
   /* allocate memory for the matrices A and R */
   double ** a = new double*[q_n];
   double ** r = new double*[n];
   for(i = 0; i < n; i++) 
   {
        a[i] = new double[m];
        r[i] = new double[r_m];
   }
   for(; i < q_n; i++) 
   {
        a[i] = new double[m];
   }
   
   a = Eigen2std(U_D,N,N_p,q_n,r_m);
   
   Eigen::MatrixXd q(std2Eigen(a,N,N_p));
   Eigen::MatrixXd R(std2Eigen(r,N_p,N_p));
   
   gramSchmidt(a, r, m, n, full);
   q = std2Eigen(a,N,N_p);
   U_D = q;			
   Info << "Finished Reorthonormalisation at" << T << endl;

   // Saving U_d in U_D at the end of the time segment
   for(int i=0;i<N;i++)
   {
   	for(int j=0;j<N_p;j++)
   	{
   		U_D_store << U_D(i,j) << "\t";
   	}
   	U_D_store << endl;
   }
   		
   R = std2Eigen(r,N,N_p);
   for(label i=0; i<N_p; i++)
   {
	R_full << R.row(i) << endl;		
   }
   
   std::ofstream LE_t;
   LE_t.open("$HOME/OpenFOAM/sidsahu-10/SC/run/LE.txt", std::ios_base::app); // append instead of overwrite
   LE_t << T << "\t";     		
   for(label ii=0; ii<N_p; ii++)
   { 
     	cum(ii) = cum(ii) + Foam::log(Foam::mag(r[ii][ii]));
      	lp(ii)  = cum(ii)/T;
      	LE_t << lp(ii) << "\t";
   }  
   LE_t << endl;
   LE_t.close();
   
   //Updating cum
   std::ofstream cum1("$HOME/OpenFOAM/sidsahu-10/SC/run/cum.txt");
   for (int i = 0; i < N_p; i++) 
   {
       cum1 << cum(i) << "\t";
   }
   
   #include "postProcess.H"
   #include "setRootCaseLists.H"

   #include "createTime.H"
   #include "createMesh.H"
   #include "initContinuityErrs.H"
   #include "createDyMControls.H"
   #include "createFields.H"
   #include "createUfIfPresent.H"
  
   k = 0;
   forAll(U_d, id)
   {	
    	U_d[id].component(0) = U_D(k,0);k++; 
    	U_d[id].component(1)= U_D(k,0);k++; 	 	  
   }
   k = 0;
   forAll(U_d2, id)
   {	
    	U_d2[id].component(0) = U_D(k,1);k++; 
    	U_d2[id].component(1)= U_D(k,1);k++; 	 	  
   }
   k = 0;
   forAll(U_d3, id)
   {	
    	U_d3[id].component(0) = U_D(k,2);k++; 
    	U_d3[id].component(1)= U_D(k,2);k++; 	 	  
   }   
   k = 0;
   forAll(U_d4, id)
   {	
    	U_d4[id].component(0) = U_D(k,3);k++; 
    	U_d4[id].component(1)= U_D(k,3);k++; 	 	  
   }    	
   runTime.write();
	
    
    
    /* free memory */
    for(i = 0; i < n; i++) {
        delete[] a[i];
        delete[] r[i];
    }
    for(; i < q_n; i++) {
        delete[] a[i];
    }
    delete[] a;  
    delete[] r;
    

    return 0;
}


// ************************************************************************* //
