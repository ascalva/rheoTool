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

Class
    FENE_CRA

Description
    Finite Extensible Nonlinear Elastic viscoelastic fluid model.
    Type FENE-CR(Chilcott and Rallison, 1988).

Author
    Jovani L. Favero.  All rights reserved

SourceFiles
    FENE_CRA.C

\*---------------------------------------------------------------------------*/

#ifndef FENE_CRA_H
#define FENE_CRA_H

#include "constitutiveEq.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class FENE_CRA Declaration
\*---------------------------------------------------------------------------*/

class FENE_CRA
:
    public constitutiveEq
{
    // Private data

        //- Transported viscoelastic stress
        volSymmTensorField tau_;

	volSymmTensorField A_;

	dimensionedSymmTensor I_;


        // Model constants

            //- Density
            dimensionedScalar rho_;

            //- Solvent viscosity
            dimensionedScalar etaS_;

            //- Zero shear rate polymer viscosity
            dimensionedScalar etaP_;

            //- Molecular adimensional extensibility
            dimensionedScalar b_;

            //- Relaxation time
            dimensionedScalar lambda_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        FENE_CRA(const FENE_CRA&);

        //- Disallow default bitwise assignment
        void operator=(const FENE_CRA&);

    protected:

    //- Return the solvent viscosity
    virtual const dimensionedScalar etaS() const
    {
        return etaS_;
    }

    //- Return the polymeric viscosity
    virtual const dimensionedScalar etaP() const
    {
        return etaP_;
    }

public:

    //- Runtime type information
    TypeName("FENE_CRA");


    // Constructors

        //- Construct from components
        FENE_CRA
        (
            const word& name,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const dictionary& dict
        );


    // Destructor

        virtual ~FENE_CRA()
        {}


    // Member Functions

    //- Return the viscoelastic stress tensor
    virtual tmp<volSymmTensorField> tau() const
    {
        return tau_;
    }

    //- Return the coupling term for the momentum equation
    // virtual tmp<fvVectorMatrix> divTau(volVectorField& U) const;

    //- Return the density
    virtual const dimensionedScalar rho() const
    {
        return rho_;
    }

    //- Return true if GNF (non-elastic)
    virtual bool isGNF() const
    {
        return false;
    };

    //- Correct the viscoelastic stress
    virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
