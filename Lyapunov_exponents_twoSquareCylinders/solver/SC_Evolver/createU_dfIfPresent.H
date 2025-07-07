/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
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

Global
    createUfIfPresent

Description
    Creates and initialises the velocity field Uf if required.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

autoPtr<surfaceVectorField> U_df;

if (mesh.dynamic())
{
    Info<< "Constructing face velocity U_df\n" << endl;

    U_df = new surfaceVectorField
    (
        IOobject
        (
            "U_df",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::interpolate(U_d)
    );
}

autoPtr<surfaceVectorField> U_d2f;

if (mesh.dynamic())
{
    Info<< "Constructing face velocity U_d2f\n" << endl;

    U_d2f = new surfaceVectorField
    (
        IOobject
        (
            "U_d2f",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::interpolate(U_d2)
    );
}

// ************************************************************************* //

autoPtr<surfaceVectorField> U_d3f;

if (mesh.dynamic())
{
    Info<< "Constructing face velocity U_d3f\n" << endl;

    U_d3f = new surfaceVectorField
    (
        IOobject
        (
            "U_d3f",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::interpolate(U_d3)
    );
}

// ************************************************************************* //

autoPtr<surfaceVectorField> U_d4f;

if (mesh.dynamic())
{
    Info<< "Constructing face velocity U_d4f\n" << endl;

    U_d4f = new surfaceVectorField
    (
        IOobject
        (
            "U_d4f",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::interpolate(U_d4)
    );
}

// ************************************************************************* //

autoPtr<surfaceVectorField> U_d5f;

if (mesh.dynamic())
{
    Info<< "Constructing face velocity U_d5f\n" << endl;

    U_d5f = new surfaceVectorField
    (
        IOobject
        (
            "U_d5f",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::interpolate(U_d5)
    );
}

// ************************************************************************* //

autoPtr<surfaceVectorField> U_d6f;

if (mesh.dynamic())
{
    Info<< "Constructing face velocity U_d6f\n" << endl;

    U_d6f = new surfaceVectorField
    (
        IOobject
        (
            "U_d6f",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::interpolate(U_d6)
    );
}

