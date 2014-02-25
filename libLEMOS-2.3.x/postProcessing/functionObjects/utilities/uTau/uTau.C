/*---------------------------------------------------------------------------*\ 
| File modified by LEMOS (University of Rostock) 2013                         |
\*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "uTau.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "incompressible/turbulenceModel/turbulenceModel.H"
#include "compressible/turbulenceModel/turbulenceModel.H"
#include "wallPolyPatch.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


namespace Foam
{
defineTypeNameAndDebug(uTau, 0);
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::uTau::writeFileHeader(const label i)
{
    // Add headers to output data
    file() <<  token::TAB << "# Friction velocity based on sqrt(mag(snGrad(Reff)))" 
           <<  token::TAB << token::TAB 
	   << "# Friction velocity based on sqrt(nu*mag(snGrad(U)))" << nl
           <<  "# time " 
           << token::TAB << token::TAB << "uTauMin" 
           << token::TAB << token::TAB << "uTauMax" 
           << token::TAB << token::TAB << "uTauAvg" 
	   <<  token::TAB << token::TAB 
           << token::TAB << token::TAB << "uTauMin"
           << token::TAB << token::TAB << "uTauMax" 
           << token::TAB << token::TAB << "uTauAvg"
	   << endl;
}


void Foam::uTau::calcFrictionVelocity
(
    const fvMesh& mesh,
    const volSymmTensorField& Reff
)
{
    uTauMin = 0.0;
    uTauMax = 0.0;
    uTauAvg = 0.0;

    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        label patchI = iter.key();
        const polyPatch& pp = mesh.boundaryMesh()[patchI];

        const vectorField& Sfp = mesh.Sf().boundaryField()[patchI];
        const scalarField& magSfp = mesh.magSf().boundaryField()[patchI];
        const symmTensorField& Reffp = Reff.boundaryField()[patchI];

        vectorField ssp = (-Sfp/magSfp) & Reffp;

        scalar uTauMinp = gMin(mag(ssp));
        uTauMin += uTauMinp;

        scalar uTauMaxp = gMax(mag(ssp));
        uTauMax += uTauMaxp;

        scalar uTauAvgp = gAverage(mag(ssp));
        uTauAvg += uTauAvgp;

        if (log_ > 2)
        {
            Info<< "  Reff:  min/max/avg(" << pp.name() << ") = " 
                << uTauMinp << ", " << uTauMaxp << ", " << uTauAvgp << endl;
        }
    }


    scalar numWallPatches = patchSet_.size();

    uTauMin = sqrt(uTauMin/numWallPatches);
    uTauMax = sqrt(uTauMax/numWallPatches);
    uTauAvg = sqrt(uTauAvg/numWallPatches);
        
    if (Pstream::master())
    {
            file() << mesh.time().timeName() << token::TAB
                << token::TAB << uTauMin << token::TAB << uTauMax
                << token::TAB << uTauAvg  << token::TAB;
    }


     if (log_ == 1)
     {
	Info<< "   Reff:  avg = " << uTauAvg << endl;
     }
     else  if (log_ == 2)
     {
	Info<< "   Reff:  min/max/avg = " 
                << uTauMin << ", " << uTauMax 
		<< ", " << uTauAvg << endl;
     }

}


void Foam::uTau::calcFrictionVelocity
(
    const fvMesh& mesh,
    const volVectorField& U,
    const volScalarField& nu
)
{
    
    uTauMin = 0.0;
    uTauMax = 0.0;
    uTauAvg = 0.0;

    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        label patchI = iter.key();
        const polyPatch& pp = mesh.boundaryMesh()[patchI];


        const vectorField wallGradU = -U.boundaryField()[patchI].snGrad();

        scalar uTauMinp = gMin(nu*mag(wallGradU));
        uTauMin += uTauMinp;

        scalar uTauMaxp = gMax(nu*mag(wallGradU));
        uTauMax += uTauMaxp;

        scalar uTauAvgp = gAverage(nu*mag(wallGradU));
        uTauAvg += uTauAvgp;

        if (log_ > 2)
        {
            Info<< " gradU:  min/max/avg(" << pp.name() << ") = "
                << uTauMinp << ", " << uTauMaxp << ", " << uTauAvgp << endl;
        }
    }

    scalar numWallPatches = patchSet_.size();

    uTauMin = sqrt(uTauMin/numWallPatches);
    uTauMax = sqrt(uTauMax/numWallPatches);
    uTauAvg = sqrt(uTauAvg/numWallPatches);

    if (Pstream::master())
    {
            file() << token::TAB 
                << token::TAB << uTauMin << token::TAB << uTauMax
                << token::TAB << uTauAvg  << token::TAB << endl;
    }

   if (log_ == 1)
   {
   	Info<< "  gradU:  avg = " << uTauAvg << endl;
   }
   else  if (log_ == 2)
   {
        Info<< "  gradU:  min/max/avg = "
                << uTauMin << ", " << uTauMax
                << ", " << uTauAvg << endl;
   }

}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uTau::uTau
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    functionObjectFile(obr, name, typeName),
    name_(name),
    obr_(obr),
    active_(true),
    log_(0),
    patchSet_()
{
    // Check if the available mesh is an fvMesh, otherwise deactivate
    if (!isA<fvMesh>(obr_))
    {
        active_ = false;
        WarningIn
        (
            "uTau::uTau"
            "("
                "const word&, "
                "const objectRegistry&, "
                "const dictionary&, "
                "const bool"
            ")"
        )   << "No fvMesh available, deactivating." << nl
            << endl;
    }

    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::uTau::~uTau()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::uTau::read(const dictionary& dict)
{
    if (active_)
    {
        log_ = dict.lookupOrDefault<label>("log", 0);


        const fvMesh& mesh = refCast<const fvMesh>(obr_);
        const polyBoundaryMesh& pbm = mesh.boundaryMesh();

        patchSet_ =
            mesh.boundaryMesh().patchSet
            (
                wordReList(dict.lookupOrDefault("patches", wordReList()))
            );

        Info<< type() << " output:" << nl;

        if (patchSet_.empty())
        {
            forAll(pbm, patchI)
            {
                if (isA<wallPolyPatch>(pbm[patchI]))
                {
                    patchSet_.insert(patchI);
                }
            }

            Info<< "    processing all wall patches" << nl << endl;
        }
        else
        {
            Info<< "    processing wall patches: " << nl;
            labelHashSet filteredPatchSet;
            forAllConstIter(labelHashSet, patchSet_, iter)
            {
                label patchI = iter.key();
                if (isA<wallPolyPatch>(pbm[patchI]))
                {
                    filteredPatchSet.insert(patchI);
                    Info<< "        " << pbm[patchI].name() << endl;
                }
                else
                {
                    WarningIn("void uTau::read(const dictionary&)")
                        << "Requested friction velocity on non-wall boundary "
                        << "type patch: " << pbm[patchI].name() << endl;
                }
            }

            Info<< endl;

            patchSet_ = filteredPatchSet;
        }

    }
}


void Foam::uTau::execute()
{
    // Do nothing - only valid on write
}


void Foam::uTau::end()
{
    // Do nothing - only valid on write
}

void Foam::uTau::timeSet()
{
    // Do nothing
}

void Foam::uTau::write()
{
    typedef compressible::turbulenceModel cmpModel;
    typedef incompressible::turbulenceModel icoModel;

    if (active_)
    {
        functionObjectFile::write();

        const fvMesh& mesh = refCast<const fvMesh>(obr_);

        if (log_ > 0)
        {
            Info<< type() << " output:" << nl;
        }


        tmp<volVectorField> U = mesh.lookupObject<volVectorField>("U");


        if (mesh.foundObject<cmpModel>("turbulenceModel"))
        {
            const cmpModel& model =
                mesh.lookupObject<cmpModel>("turbulenceModel");
        
	    calcFrictionVelocity(mesh, model.devRhoReff());
            calcFrictionVelocity(mesh, U(), model.mu()/model.rho());
        }
        else if (mesh.foundObject<icoModel>("turbulenceModel"))
        {
            const icoModel& model =
                mesh.lookupObject<icoModel>("turbulenceModel");

	    calcFrictionVelocity(mesh, model.devReff());
            calcFrictionVelocity(mesh, U(), model.nu());


        }
        else
        {
            FatalErrorIn("void Foam::uTau::write()")
                << "Unable to find turbulence model in the "
                << "database" << exit(FatalError);
        }

        if (log_ > 0)
        {
            Info<< endl;
        }

    }
}


// ************************************************************************* //
