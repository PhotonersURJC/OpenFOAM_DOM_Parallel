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
    Foam::radiation::specularGuide

Description
    Class to storage specular indexes for a given face in a region.
	Stores a specularIndexer per face direction (a flat face stores only
	a single specularIndexer).

SourceFiles
    specularGuide.C

\*---------------------------------------------------------------------------*/

#ifndef radiationSpecularGuide_H
#define radiationSpecularGuide_H

#include "specularIndexer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{


/*---------------------------------------------------------------------------*\
                    Class specularGuide Declaration
\*---------------------------------------------------------------------------*/

class specularGuide
{

        //- Reference to the patch
        const fvPatch& patch_;
		
		List<short> indexes_;

		PtrList<specularIndexer> interpolators_;
		
		//- Intensity in [band][direction][face]
        List<PtrList<scalarField> > storedIntensities_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        specularGuide(const specularGuide&);

        //- Disallow default bitwise assignment
        void operator=(const specularGuide&);

public:

    // Constructors

        //- Construct from components
        specularGuide
        (
            const Foam::fvPatch& patch,
			const Foam::List<Foam::vector> mainDirs,
			const List<List< label > > facesNeighbours,
			const List<scalar> faceMagnitudes,
			const List<vector> quadVertex,
			const List<face> quadFaces,
			const bool secondOrder,
			const label nLambda
        );


    //- Destructor
    ~specularGuide();


    // Member functions
		//- Assignment
			void updateSpecular(const PtrList<scalarField> patchIntensity, const label lambdaI);

	    //- Return the reflecteds Intensities for a given face
            scalarField getSpec(const label rayId, const label lambdaId) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
