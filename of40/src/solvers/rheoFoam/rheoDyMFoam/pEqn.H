/**
 * author: Alberto Serrano - ascalva@gmail.com
 */

{
    volScalarField rAU(1.0/UEqn.A());
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
      + fvc::interpolate(rAU)*fvc::ddtCorr(U, Uf)
    );

    // MRF.makeRelative(phiHbyA);

    if (p.needReference())
    {
        fvc::makeRelative(phiHbyA, U);
        adjustPhi(phiHbyA, U, p);
        fvc::makeAbsolute(phiHbyA, U);
    }

    tmp<volScalarField> rAtU(rAU);

    rAtU = 1.0/(1.0/rAU - UEqn.H1());
    phiHbyA += fvc::interpolate(rAtU() - rAU)*fvc::snGrad(p)*mesh.magSf();
    HbyA -= (rAU - rAtU())*fvc::grad(p);

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p, U, phiHbyA, rAtU());

    tUEqn.clear();

    // Non-orthogonal pressure corrector loop
    while( simple.lookup("correctNonOrthogonal") )
    {
        fvScalarMatrix pEqn
        (
           fvm::laplacian(rAtU(), p, "laplacian(p|(ap-H1))") == fvc::div(phiHbyA)
        );

        pEqn.setReference(pRefCell, pRefValue);

        pEqn.solve();

        if( simple.lookup("finalNonOrthogonalIter") )
        {
            phi = phiHbyA - pEqn.flux();
        }
    }

    #include "continuityErrs.H"

    // Explicitly relax pressure for momentum corrector
    p.relax();

    U = HbyA - rAtU*fvc::grad(p);
    U.correctBoundaryConditions();
    fvOptions.correct(U);

    {
        Uf = fvc::interpolate(U);
        surfaceVectorField n(mesh.Sf()/mesh.magSf());
        Uf += n*(phi/mesh.magSf() - (n & Uf));
    }

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi, U);
}
