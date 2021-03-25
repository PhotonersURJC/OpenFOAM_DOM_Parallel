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

#ifndef radiationSpecularMatrix_H
#define radiationSpecularMatrix_H

#include "volFields.H"
#include "PatchToPatchInterpolation.H"
//#include "MeshedSurface.H"
#include "PrimitivePatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{


/*---------------------------------------------------------------------------*\
                    Class specularIndexer Declaration
\*---------------------------------------------------------------------------*/

class specularMatrix
{
				
		PrimitivePatch<face,List,vectorField,point> fromPatch_;
		
		PrimitivePatch<face,List,vectorField,point> toPatch_;
		
		PatchToPatchInterpolation<PrimitivePatch<face,List,vectorField,point>,PrimitivePatch<face,List,vectorField,point> > interpolator_;
		
		List<bool> dirOut_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        specularMatrix(const specularMatrix&);

        //- Disallow default bitwise assignment
        void operator=(const specularMatrix&);


public:

    // Constructors
		
        //- Construct from components
        specularMatrix
        (
			const PrimitivePatch<face,List,pointField,point> quadraturePatch,
			const PrimitivePatch<face,List,pointField,point> reflectedPatch,
			const List<vector> mainDirs,
			const vector faceN
        );

    //- Destructor
    ~specularMatrix();


    // Member functions

	    //- Return the specular index for a given direction
            scalarField getSpec(const scalarField originalI) const;
			
		//- Determines if a direction is leaving the face, or not
			bool dirOut(const label dirI) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
