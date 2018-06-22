/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "FENE_PA.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(FENE_PA, 0);
    addToRunTimeSelectionTable(constitutiveEq, FENE_PA, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FENE_PA::FENE_PA
(
    const word& name,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const dictionary& dict
)
:
    constitutiveEq(name, U, phi),
    tau_
    (
        IOobject
        (
            "tau" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
	dimensionedSymmTensor
        (
            "zero",
            dimensionSet(1, -1, -2, 0, 0, 0, 0),
            symmTensor::zero
        )
    ),
    A_
    (
        IOobject
        (
            "A" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()//,
        // dimensionedSymmTensor("I", dimless, symmTensor::I),
    	// tau_.boundaryField().types()
    ),
    I_
    (
        dimensionedSymmTensor
        (
            "I",
            dimensionSet(0, 0, 0, 0, 0, 0, 0),
            symmTensor
            (
                1, 0, 0,
                   1, 0,
                      1
            )
        )
    ),
    rho_(dict.lookup("rho")),
    etaS_(dict.lookup("etaS")),
    etaP_(dict.lookup("etaP")),
    b_(dict.lookup("b")),
    lambda_(dict.lookup("lambda"))
{
    checkForStab(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::FENE_PA::correct()
{
    dimensionedScalar eta0 = etaP_/(etaS_+etaP_);

    // Velocity gradient tensor
    volTensorField L = fvc::grad(U());

    // Convected derivate term
    volTensorField C = A_ & L;

    // Twice the rate of deformation tensor
    volSymmTensorField twoD = twoSymm(L);

    dimensionedScalar a = b_/(b_-3);

     // Stress transport equation
    fvSymmTensorMatrix AEqn
    (
        fvm::ddt(A_)
      + fvm::div(phi(),A_)
     ==
     // ((b_ / lambda_ + tr(tau_)/etaP_)/(b_ - 3.0))*etaP_*twoD
        twoSymm(C)
     - fvm::Sp(1/lambda_*(b_/(b_-tr(A_))), A_)
     + 1/lambda_*a*I_
    );

    AEqn.relax();
    AEqn.solve();

    tau_ = etaP_/lambda_*((b_/(b_-tr(A_)))*A_-a*I_);

    // NOTE: Did not compile correctly with line under!
    // tau_.correctBoundaryCondition();
}


// ************************************************************************* //
