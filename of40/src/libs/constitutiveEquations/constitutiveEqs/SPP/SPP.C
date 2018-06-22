/*
License
    GPL3+
*/

#include "SPP.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam {
    defineTypeNameAndDebug(SPP, 0);
    addToRunTimeSelectionTable(constitutiveEq, SPP, dictionary);
}

Foam::SPP::SPP (
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
    Y_
    (
        IOobject
        (
            "Y" + name,
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
    K_
    (
         dimensionedSymmTensor
         (
          "K",
          dimensionSet(0, 0, 0, 0, 0, 0, 0),
          symmTensor
          (
               -1, 0, 0,
                   1, 0,
                      1
           )
          )
     ),
    rho_(dict.lookup("rho")),
    etaS_(dict.lookup("etaS")),
    etaP_(dict.lookup("etaP")),
    c1(dict.lookup("c1")),
    c2(dict.lookup("c2")),
    c3(dict.lookup("c3")),
    c4(dict.lookup("c4")),
    c5(dict.lookup("c5")),
    c6(dict.lookup("c6")),
    c7(dict.lookup("c7")),
    c8(dict.lookup("c8")),
    c9(dict.lookup("c9")),
    k1(dict.lookup("k1")),
    k2(dict.lookup("k2")),
    k3(dict.lookup("k3")),
    k4(dict.lookup("k4")),
    k5(dict.lookup("k5")),
    k6(dict.lookup("k6")),
    k7(dict.lookup("k7")),
    phi0(dict.lookup("phi0")),
    phim0(dict.lookup("phim0"))

{
    checkForStab(dict);
}

// Foam::tmp<Foam::fvVectorMatrix>
// Foam::SPP::divTau(volVectorField& U) const {
//     dimensionedScalar etaPEff = etaP_;
//
//     return
//     (
//         fvc::div(tau_/rho_, "div(tau)")
// //      + k7*twoSymm(fvc::grad(U))*fvc::div(Y_/rho_)+fvm::laplacian(k7/rho_*Y_,U)
//       - fvc::div((etaPEff/rho_)*fvc::grad(U))
// //      - fvc::laplacian(etaPEff/rho_, U, "laplacian(etaPEff,U)")
//       + fvm::laplacian( (etaPEff + etaS_)/rho_, U, "laplacian(etaPEff+etaS,U)")
// //      + fvm::laplacian( (etaS_)/rho_, U, "laplacian(etaPEff+etaS,U)")
//     );
// }

void Foam::SPP::correct() {

    dimensionedScalar eta0 = etaP_/(etaS_+etaP_);

    // Velocity gradient tensor
    volTensorField L = fvc::grad(U());

    volSymmTensorField Gdot = symm(L);

    volTensorField LT = 2*Gdot-L;

    volTensorField Omega = 1/2*(L-LT);

    volTensorField YL = Y_ & L;

    volTensorField LY = L & Y_;

    volSymmTensorField YY = symm(K_ & Y_ & Y_);

    volSymmTensorField OmegaY = symm(LY)-symm(YL);

    volTensorField YGdot1 = Y_ & Gdot;

    volSymmTensorField YGdot = twoSymm(YL);//(YGdot1);

    dimensionedScalar fphi = Foam::pow((1/phi0-1/phim0),-1);

    volScalarField gdot = Foam::sqrt(2*Gdot&&Gdot);//sqrt(2*L&&L);//


    //If L is 0, gdot is 0, leads to a division by zero,
    volScalarField h = Foam::pow(gdot,0);//-1/8);

    volScalarField gdoth = Foam::pow(gdot,1);//7/8);

    // Stress transport equation
    fvSymmTensorMatrix YEqn
    (
        fvm::ddt(Y_)
      + fvm::div(phi(), Y_,"div(phi,Y)")
      ==
      - 1/2*twoSymm(LY)+1/2*twoSymm(YL)
      + gdoth*(c2/3*tr(Y_) + c1*fphi)*I_ + c3*tr(YGdot1)*I_
     + fvm::Sp(c4*gdoth,Y_) -c4*gdoth*1/3*tr(Y_)*I_ + (c5*fphi + c6*1/3*tr(Y_))*Gdot
      + c7*(YGdot - 2/3*tr(YGdot1)*I_)

    );

    YEqn.relax();
    YEqn.solve();
//k2*gdot*(1/3*tr(Y_)-3*fphi)*I_
    tau_ = etaP_*(k2*gdot*(1/3*tr(Y_)-3*fphi)*I_
    + h*k3*tr(YGdot1)*I_
    + k4*gdot*Y_
    + h*(k5*fphi + k6/3*tr(Y_))*Gdot
    + h*k7*gdot*YY);//(YGdot - 1/3*tr(YGdot)*I_));

    // tau_.correctBoundaryConditions();
}
