/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Application
    rheoDyMFoam.C

Description
    Transient solver for incompressible, laminar flow. Any GNF or viscoelastic
    model of library lconstitutiveEquations can be selected. Pressure-velocity
    coupling is using the SIMPLEC algorithm.

    Dynamic mesh capabilities have been added to original rheoFoam solver.

Changes made by: Alberto Serrano - ascalva@gmail.com

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "IFstream.H"
#include "OFstream.H"
#include "dynamicFvMesh.H"
#include "simpleControl.H"
#include "CorrectPhi.H"
#include "fvOptions.H"

#include "extrapolatedCalculatedFvPatchField.H"
#include "ppUtilInterface.H"
#include "constitutiveModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createControls.H"
    #include "createFields.H"
    #include "createUf.H"
    #include "createFvOptions.H"
    #include "createPPutil.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    // Read extra-controls

    auto   simple = mesh.solutionDict().subDict("SIMPLE");
    int    nInIter = simple.lookupOrDefault("nInIter", 1);
    bool   sPS = cttProperties.subDict("passiveScalarProperties").lookupOrDefault<Switch>("solvePassiveScalar", false);
    if (sPS) C.writeOpt() = IOobject::AUTO_WRITE;

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    // --- Time loop ---

    while( simple.lookup("loop") )
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "readControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        // --- Inner loop iterations ---

        /// Adaptive/Dynamic mesh section
        /******** START ********/

        mesh.update();

        // Calculate absolute flux from the mapped surface velocity
        phi = mesh.Sf() & Uf;

        if( mesh.changing() && correctPhi )
        {
            #include "CorrectPhi.H"
        }

        fvc::makeRelative(phi, U);

        if( mesh.changing() && checkMeshCourantNo )
        {
            #include "meshCourantNo.H"
        }

        /******** END ********/

        for (int i=0; i<nInIter; i++)
	  {

            Info<< "Inner iteration:  " << i << nl << endl;

            // --- Pressure-velocity SIMPLEC corrector
            {
               // ---- Solve constitutive equation ----
               constEq.correct();

               // ---- Solve U and p ----
               #include "UEqn.H"
               #include "pEqn.H"
            }

            // --- Passive Scalar transport
            if (sPS)
             {
               #include "CEqn.H"
             }
         }

        postProc.update();
        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
