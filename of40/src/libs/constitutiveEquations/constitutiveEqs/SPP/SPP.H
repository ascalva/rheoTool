/*
License
    GPL3+

Class
    SPP

Description
    ROuse LInear Elastic POLYmer viscoelastic fluid model (need citation).

Author
    Daniel James Ryckman

SourceFiles
    SPP.C
*/

#ifndef SPP_H
#define SPP_H

#include "constitutiveEq.H"

namespace Foam {

class SPP:
    public constitutiveEq {

        // Tensors

        // Viscoelastic stress tensor
        volSymmTensorField tau_;

        // The Y tensor
        volSymmTensorField Y_;

        // Identity tensor
        dimensionedSymmTensor I_;
        dimensionedSymmTensor K_;

        // Constants

        // Density
        dimensionedScalar rho_;

        //- Solvent viscosity
        dimensionedScalar etaS_;

        //- Zero shear rate polymer viscosity
        dimensionedScalar etaP_;

	//Other Parameters
	dimensionedScalar c1;

	dimensionedScalar c2;

	dimensionedScalar c3;

	dimensionedScalar c4;

	dimensionedScalar c5;

	dimensionedScalar c6;

	dimensionedScalar c7;

	dimensionedScalar c8;

	dimensionedScalar c9;

	dimensionedScalar k1;

	dimensionedScalar k2;

	dimensionedScalar k3;

	dimensionedScalar k4;

	dimensionedScalar k5;

	dimensionedScalar k6;

	dimensionedScalar k7;

	dimensionedScalar phi0;

        dimensionedScalar phim0;

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
        TypeName("SPP");

        // Constructor
        SPP (
            const word& name,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const dictionary dict
        );

        // Destructor
        virtual ~SPP() {}

        // Private member functions

        // Return the viscoelastic stress tensor
        virtual tmp<volSymmTensorField> tau() const { return tau_; };

        // Return the coupling term in the momentum equation
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
        // Correct the viscoelastic stress term
        virtual void correct();
    };
}

#endif
