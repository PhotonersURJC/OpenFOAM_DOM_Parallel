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
    Foam::radiation::discreteOrdinate

Description
    Abstract class for radiative intensity storage in a given solid angle.
	May be:
		- cleanDiscOrd: Transparent media
		- homoDiscOrd: Absorptive and scattering media, homogeneous coeffs
		- fullDiscOrd: same as homoDiscOrd, but heterogeneous coeffs

SourceFiles
    discreteOrdinate.C

\*---------------------------------------------------------------------------*/

#ifndef discreteOrdinate_H
#define discreteOrdinate_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

// Forward declaration of classes
class DORT;
class quadrature;


/*---------------------------------------------------------------------------*\
                    Class discreteOrdinate Declaration
\*---------------------------------------------------------------------------*/

class discreteOrdinate
{
public:

    static const word intensityPrefix;


protected:

    // Private data
	
        //- Refence to the owner DORT object
        const DORT& dort_;

	//- Reference to the owner quadrature
	const quadrature& quad_;

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Direction
        vector d_;

        //- Average direction vector inside the solid angle
        vector dAve_;

        //- Theta angle
        scalar theta_;

        //- Phi angle
        scalar phi_;

        //- Solid angle
        scalar omega_;

        //- Number of wavelengths/bands
        label nLambda_;

        //- List of pointers to radiative intensity fields for given wavelengths
        PtrList<volScalarField> ILambda_;
		
        //- Global ray id - incremented in constructor
        static label rayId;

        //- My ray Id
        label myRayId_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        discreteOrdinate(const discreteOrdinate&);

        //- Disallow default bitwise assignment
        void operator=(const discreteOrdinate&);


public:

    // Constructors

        //- Construct form components
        discreteOrdinate
        (
            const DORT& dort,
	    const quadrature& quad,
            const fvMesh& mesh,
            const scalar phi,
            const scalar theta,
            const scalar deltaPhi,
            const scalar deltaTheta,
            const label lambda,
            const label rayId,
	    const List<vector> mainAxis
        );


    //- Destructor
    virtual ~discreteOrdinate();


    // Member functions

        // Edit

        //- Update radiative intensity on i direction
        virtual scalar correct() = 0;

	    virtual void setBandScatter (const label lambdaI, const volScalarField value) = 0;
		
		virtual void changeBand (const label lambdaI) = 0;

        // Access

            //- Return direction
            inline const vector& d() const;

            //- Return the average vector inside the solid angle
            inline const vector& dAve() const;

            //- Return the number of bands
            inline scalar nLambda() const;

            //- Return the phi angle
            inline scalar phi() const;

            //- Return the theta angle
            inline scalar theta() const;

            //- Return the solid angle
            inline scalar omega() const;

            //- Return the radiative intensity for a given wavelength
            inline const volScalarField& ILambda(const label lambdaI) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "discreteOrdinateI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
