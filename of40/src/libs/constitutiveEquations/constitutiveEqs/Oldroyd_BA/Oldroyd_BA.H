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
    Oldroyd_BA

Description
    Oldroyd-B non linear viscoelastic fluid model (Oldroyd, J.G., 1950).

Author
    Jovani L. Favero.  All rights reserved

SourceFiles
    Oldroyd_BA.C

\*---------------------------------------------------------------------------*/

#ifndef Oldroyd_BA_H
#define Oldroyd_BA_H

#include "constitutiveEq.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class Oldroyd_BA Declaration
\*---------------------------------------------------------------------------*/

class Oldroyd_BA
:
    public constitutiveEq
{
    // Private data

    //- Tensor known as A
	volSymmTensorField A_;
    
	//- Transported viscoelastic stress
	volSymmTensorField tau_;

	//- Identity tensor
	dimensionedSymmTensor I_;

	//Gradient of U
//	volTensorField L;


	// Model constants

	//- Density
	dimensionedScalar rho_;

	//- Solvent viscosity
	dimensionedScalar etaS_;

	//- Zero shear rate polymer viscosity
	dimensionedScalar etaP_;

	//- Relaxation time
	dimensionedScalar lambda_;


    // Private Member Functions

	//- Disallow default bitwise copy construct
	Oldroyd_BA(const Oldroyd_BA&);

	//- Disallow default bitwise assignment
	void operator=(const Oldroyd_BA&);


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
    TypeName("Oldroyd_BA");

    // Constructors

	//- Construct from components
	Oldroyd_BA
	(
		const word& name,
		const volVectorField& U,
		const surfaceScalarField& phi,
		const dictionary& dict
	);


    // Destructor

	virtual ~Oldroyd_BA()
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
