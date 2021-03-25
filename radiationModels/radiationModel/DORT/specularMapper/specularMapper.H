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
    Radiation intensity for a ray in a given direction

SourceFiles
    discreteOrdinate.C

\*---------------------------------------------------------------------------*/

#ifndef radiationSpecularIndexer_H
#define radiationSpecularIndexer_H

//#include "fvPatch.H"
//#include "fvMatrices.H"
#include "specularMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{


/*---------------------------------------------------------------------------*\
                    Class specularIndexer Declaration
\*---------------------------------------------------------------------------*/

class specularIndexer
{

        //- Reference to the patch
        const fvPatch& patch_;
		
		List<short> indexes_;
		
        //- Intensity in [band][direction][face]
        List<PtrList<scalarField> > storedIntensities_;

		PtrList<specularMatrix> interpolators_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        specularIndexer(const specularIndexer& /*toCopy*/);

        //- Disallow default bitwise assignment
        void operator=(const specularIndexer&);

public:

    // Constructors
	
		//- Null constructor
		/*specularIndexer();*/

        //- Construct from components
        specularIndexer
        (
            const fvPatch& patch,
            const List<vector> mainDirs,
			const List<vector> vertex,
			const List<face> faces,
			const label nLambda
        );


    //- Destructor
    ~specularIndexer();


    // Member functions
		//- Assignment
			void updateSpecular(const PtrList<scalarField> patchIntensity, const label lambdaI);

	    //- Return the specular index for a given face
            scalarField getSpec(const label rayId, const label lambdaId) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
