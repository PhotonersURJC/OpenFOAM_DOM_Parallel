/*---------------------------------------------------------------------------*\
             Discrete Ordinate Method Radiation Model for OpenFOAM.
	    V 2.0. Including optimization and shared memory parallelization


Code corresponding to the article entitled:


"Optimization and parallelization of the Discrete Ordinate Method
     for radiation transport simulation in OpenFOAM: Hierarchical 
          combination of shared and distributed memory approaches"


by:


Jose Moreno-SanSegundo, Cintia Casado, Javier Marugán from

    Department of Chemical and Environmental Technology,

    ESCET, Universidad Rey Juan Carlos,

and David Concha, Antonio Sanz from

	Department of Computer Science and Software Architecture,

C/Tulipán s/n, 28933 Móstoles (Madrid), Spain

Tel. +34 91 664 7466; E-mail: javier.marugan@urjc.es

License
	This file is is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
	You should have received a copy of the GNU General Public License
    along with the code.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::radiation::sourceGreyRadiationMixedFvPatchScalarField

Group
    grpThermoBoundaryConditions

Description
    This boundary condition provides a grey-diffuse condition for radiation
    intensity, \c I, for use with the finite-volume discrete-ordinates model
    (fvDOM), in which the radiation temperature is retrieved from the
    temperature field boundary condition.

    \heading Patch usage

    \table
        Property     | Description					|Required|Default value
        refCoeff     | reflection Coefficient		| no     | 0.0
		refDifF      | reflection Diffuse fraction	| no     | 0.0
		transCoeff   | transmission Coefficient 	| no     | 0.0
		transDifF    | transmission Diffuse fraction| no     | 0.0
		sampleRegion | neighbour region,transmission| no     | empty
		samplePatchI | patch in neighbour region	| no     | empty
		order        | scheme order in reflection	| no     | first
    \endtable
	
	
	
	
    Example of the boundary condition specification:
	(full specular reflection)
    \verbatim
    myPatch
    {
        type            sourceGreyRadiation;
        refCoeff		1.0;
        value           uniform 0;
    }
    \endverbatim

SeeAlso
    Foam::radiation::radiationModel
    Foam::radiation::DORT
    Foam::mixedFvPatchField

SourceFiles
    sourceGreyRadiationMixedFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef sourceGreyRadiationMixedFvPatchScalarField_H
#define sourceGreyRadiationMixedFvPatchScalarField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{
/*---------------------------------------------------------------------------*\
            Class sourceGreyRadiationMixedFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class sourceGreyRadiationMixedFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private data
        //- Underlying patch
        const fvPatch& patch_;
		const DORT& dort;
		quadrature& quad;
		label rayId;
		label lambdaId;
		label patchI;
		vector direction_;

	scalar refCoeff_;	//Reflexion Coefficient
	scalar refDifF_;
	scalar transCoeff_;
	scalar transDifF_;
	word sampleMesh_;
	word samplePatchI_;
	bool secondOrder_;

	void initialise();
	

public:

    //- Runtime type information
    TypeName("sourceGreyRadiation");


    // Constructors

        //- Construct from patch and internal field
        sourceGreyRadiationMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        sourceGreyRadiationMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given a
        //  sourceGreyRadiationMixedFvPatchScalarField onto a new patch
        sourceGreyRadiationMixedFvPatchScalarField
        (
            const sourceGreyRadiationMixedFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        sourceGreyRadiationMixedFvPatchScalarField
        (
            const sourceGreyRadiationMixedFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new sourceGreyRadiationMixedFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        sourceGreyRadiationMixedFvPatchScalarField
        (
            const sourceGreyRadiationMixedFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new sourceGreyRadiationMixedFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access

	    const scalar& refCoeff() const
	    {
		return refCoeff_;
	    }
	    const scalar& refDifF() const
	    {
		return refDifF_;
	    }
	    const scalar& transCoeff() const
	    {
		return transCoeff_;
	    }
	    const scalar& transDifF() const
	    {
		return transDifF_;
	    }
	    const word& sampleMesh() const
	    {
		return sampleMesh_;
	    }
	    const word& samplePatchI() const
	    {
		return samplePatchI_;
	    }


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
