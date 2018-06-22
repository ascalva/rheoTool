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

#include "PomPom.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(PomPom, 0);
    addToRunTimeSelectionTable(constitutiveEq, PomPom, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PomPom::PomPom
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
  IOobject::NO_READ,
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
lambda_
(
 IOobject
 (
  "lambda" + name,
  U.time().timeName(),
  U.mesh(),
  IOobject::MUST_READ,
  IOobject::AUTO_WRITE
  ),
 U.mesh()
 ),
rho_(dict.lookup("rho")),
etaS_(dict.lookup("etaS")),
etaP_(dict.lookup("etaP")),
lambdaM_(dict.lookup("lambdaM")),
taub_(dict.lookup("taub")),
taus_(dict.lookup("taus"))
{
    checkForStab(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Foam::tmp<Foam::fvVectorMatrix>
// Foam::PomPom::divTau(volVectorField& U) const
// {
//     dimensionedScalar etaPEff = etaP_;
//
//     return
//     (
//      fvc::div(tau_/rho_, "div(tau)")
//      - fvc::div((etaPEff/rho_)*fvc::grad(U))
//      //      - fvc::laplacian(etaPEff/rho_, U, "laplacian(etaPEff,U)")
//      + fvm::laplacian( (etaPEff + etaS_)/rho_, U, "laplacian(etaPEff+etaS,U)")
//      //        + fvm::laplacian( etaS_/rho_, U, "laplacian(etaS,U)")
//
//      );
// }


void Foam::PomPom::correct()
{
    // Velocity gradient tensor
    volTensorField L = fvc::grad(U());

    // Convected derivate term
    volTensorField C = A_ & L;//a

    volScalarField G = A_ && L;

    // Twice the rate of deformation tensor
    volSymmTensorField twoD = twoSymm(L);

    volScalarField G2 = A_ && twoD;

    // Stress transport equation
    fvSymmTensorMatrix AEqn
    (
     fvm::ddt(A_)//a
     + fvm::div(phi(), A_)//a-i
     ==
     twoSymm(C)
     - fvm::Sp(2*G2/2,A_)
     - fvm::Sp(1/Foam::sqr(lambda_)*(1/taub_),A_)//a-i
     + 1/Foam::sqr(lambda_)*(1/taub_)*I_/3

     );

    AEqn.relax();
    AEqn.solve();

    //Equation for lambda
    fvScalarMatrix lambdaEqn
    (
     fvm::ddt(lambda_)
     + fvm::div(phi(), lambda_)
     ==
     fvm::Sp(G2/tr(A_/2),lambda_)
     +1/taus_*Foam::exp(2*(lambda_ - 1)/lambdaM_)
     -fvm::Sp(1/taus_*Foam::exp(2*(lambda_ - 1)/lambdaM_),lambda_)
     //- (lambda_ - 1)/taus_*Foam::exp(2*(lambda_ - 1)/lambdaM_)
     );
    lambdaEqn.relax();
    lambdaEqn.solve();



    tau_ = etaP_/taub_*Foam::sqr(lambda_)*(A_ - (1/3)*I_);

    // tau_.correctBoundaryConditions();
}


// ************************************************************************* //
