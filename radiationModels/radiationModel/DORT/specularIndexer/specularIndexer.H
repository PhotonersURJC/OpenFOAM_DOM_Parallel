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
    Foam::radiation::specularIndexer

Description
    Storage of specular reflection directions for each discrete Ordinate in
	a given face.

SourceFiles
    specularIndexer.C

\*---------------------------------------------------------------------------*/

#ifndef radiationSpecularIndexer_H
#define radiationSpecularIndexer_H

#include "volFields.H"

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
		//- Determines where to look for spotters, if any (else -1)
		List<short> outers_;
		
		//- Reflected incoming DOs looking at current outgoing DO
		List<List<short > > spotters_;
		
		//- Factors to add intensities of spotters
		List<List<scalar> > viewFactors_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        specularIndexer(const specularIndexer&);

        //- Disallow default bitwise assignment
        void operator=(const specularIndexer&);


public:

    // Constructors
		
        //- Construct from components
        specularIndexer
        (
			const List<List< label > > facesNeighbours,
			const List<vector> mainDirs,
			const List<scalar> faceMagnitudes,
			const vector faceN,
			const List<vector> quadVertex,
			const List<face> quadFaces,
			const bool secondOrder
        );

    //- Destructor
    ~specularIndexer();


    // Member functions

	    //- Return the reflected Field for a given face
            scalarField getSpec(const scalarField originalI) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
