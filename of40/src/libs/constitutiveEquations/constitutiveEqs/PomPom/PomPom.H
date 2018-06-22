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
    PomPom

Description
    Oldroyd-BAn non linear viscoelastic fluid model (Oldroyd, J.G., 1950).

Author
    Jovani L. Favero.  All rights reserved

SourceFiles
    PomPom.C

\*---------------------------------------------------------------------------*/

#ifndef PomPom_H
#define PomPom_H

#include "constitutiveEq.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class Oldroyd_BAn Declaration
\*---------------------------------------------------------------------------*/

class PomPom
:
    public constitutiveEq
{
    // Private data

	//- Transported viscoelastic stress
	volSymmTensorField tau_;

	//- conformation tensor A
	volSymmTensorField A_;

	//- Identity tensor
	dimensionedSymmTensor I_;

	//- number density
	volScalarField lambda_;



	// Model constants

	//- Density
	dimensionedScalar rho_;

	//- Solvent viscosity
	dimensionedScalar etaS_;

	//- Zero shear rate polymer viscosity
	dimensionedScalar etaP_;

	//- Relaxation time
	dimensionedScalar lambdaM_;

	//- Relaxation time
	dimensionedScalar taub_;

	//- Relaxation time
	dimensionedScalar taus_;

    // Private Member Functions

	//- Disallow default bitwise copy construct
	PomPom(const PomPom&);

	//- Disallow default bitwise assignment
	void operator=(const PomPom&);


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
    TypeName("PomPom");

    // Constructors

	//- Construct from components
	PomPom
	(
		const word& name,
		const volVectorField& U,
		const surfaceScalarField& phi,
		const dictionary& dict
	);


    // Destructor

	virtual ~PomPom()
	{}


    // Member Functions

	//- Return the viscoelastic stress tensor
	virtual tmp<volSymmTensorField> tau() const
	{
		return tau_;
	}

	//- Return the coupling term for the momentum equation
	// virtual tmp<fvVectorMatrix> divTau(volVectorField& U) const;

	//- Correct the viscoelastic stress
	virtual void correct();

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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
