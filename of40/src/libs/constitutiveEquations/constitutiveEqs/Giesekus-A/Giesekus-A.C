/*
 * License: GPL3+
 */

#include "Giesekus-A.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam {
    defineTypeNameAndDebug(Giesekus_A, 0);
    addToRunTimeSelectionTable(constitutiveEq, Giesekus_A, dictionary);
}

Foam::Giesekus_A::Giesekus_A (
    const word& name,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const dictionary dict
):
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
/*    I_
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
    ),*/
    rho_(dict.lookup("rho")),
    etaS_(dict.lookup("etaS")),
    etaP_(dict.lookup("etaP")),
    lambda_(dict.lookup("lambda")),
    alpha_(dict.lookup("alpha"))
{
    checkForStab(dict);
}

// Foam::tmp<Foam::fvVectorMatrix>
// Foam::Giesekus_A::divTau(volVectorField& U) const {
//     dimensionedScalar etaPeff = etaP_;
//
//     return (
// 	fvc::div(tau_ / rho_, "div(tau)")
//       - fvc::div((etaPeff / rho_) * fvc::grad(U))
//       + fvm::laplacian((etaPeff + etaS_) / rho_, U, "laplacian(etaPEff+etaS,U)")
//     );
// }

void Foam::Giesekus_A::correct() {

    // Velocity gradient tensor
    volTensorField L = fvc::grad(U());

    dimensionedSymmTensor Ist
    (
     "Identity",
     A_.dimensions(),
     symmTensor::I
     );

    // Convected derivate term
    volTensorField C = A_ & L;

    dimensionedScalar eta0 = etaP_ / (etaS_ + etaP_);

    // Stress transport equation
    fvSymmTensorMatrix aEqn
    (
        fvm::ddt(A_)
      + fvm::div(phi(), A_)
     ==
        twoSymm(C)
     - fvm::Sp(1 / lambda_,A_)
      + Ist/lambda_
      - (alpha_ / lambda_) * symm((A_ - Ist) & (A_ - Ist))
    );

    aEqn.relax();
    aEqn.solve();

    tau_ = etaP_ / lambda_ * (A_ - Ist);

//    tau_.correctBoundaryConditions();
}
