/*
 * License: GPL3+
 * Class: Giesekus-A
 * Description: Updated version of Giesekus viscoelastic fluid model
 * Author: Hector A. Escobedo
 * SourceFiles: Giesekus-A.C
 */

#ifndef GIESEKUS_A_H
#define GIESEKUS_A_H

#include "constitutiveEq.H"

namespace Foam {

class Giesekus_A:
    public constitutiveEq {

        // Tensors

        // Viscoelastic stress tensor
        volSymmTensorField tau_;

        // The conformation tensor
        volSymmTensorField A_;

        // Identity tensor
        dimensionedSymmTensor I_;

        // Constants

        // Density
        dimensionedScalar rho_;

        // Solvent viscosity
        dimensionedScalar etaS_;

        // Zero shear rate polymer viscosity
        dimensionedScalar etaP_;

        // Relaxation time
        dimensionedScalar lambda_;

        dimensionedScalar alpha_;

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
        TypeName("Giesekus_A");

        // Constructor
        Giesekus_A (
            const word& name,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const dictionary dict
        );

        // Destructor
        virtual ~Giesekus_A() {}

        // Private member functions

        // Return the viscoelastic stress tensor
        virtual tmp<volSymmTensorField> tau() const { return tau_; };

        // virtual tmp<fvVectorMatrix> divTau(volVectorField& U) const;

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
