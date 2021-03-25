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
    Foam::radiation::parallelRadSource

Description

    Parallel emmision source. Light is emmited in a single discrete ordinate.
	Quadrature is rotated to place a discrete ordinate near the quadrature
	equator in the light direction.
	
	Discrete ordinates near the equator are the biggest, and have a squared
	shape, been considered as "better" ordinates.

SourceFiles
    parallelRadSource.C

\*---------------------------------------------------------------------------*/

#ifndef radSourceparallelRadSource_H
#define radSourceparallelRadSource_H

#include "radSource.H"
#include "fvMatrices.H"
#include "booleanSourceStorage.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                    Class parallelRadSource Declaration
\*---------------------------------------------------------------------------*/

class parallelRadSource
:
    public radSource
{
	// Private data
	PtrList<booleanSourceStorage> sourceValues_;
    List<vector> patchNormals_;
    label quadIndex_;
    scalar omega_;

    // Private Member Functions

        //- Initialise
        void initialise();

        //- Disallow default bitwise copy construct
        parallelRadSource(const parallelRadSource&);

        //- Disallow default bitwise assignment
        void operator=(const parallelRadSource&);


public:

    //- Runtime type information
    TypeName("parallelRadSource");


    // Constructors

        //- Construct from components
        parallelRadSource(const fvMesh& mesh, const label sourceIndex, const label nBands);

        //- Construct from components
        parallelRadSource(const dictionary& dict, const fvMesh& mesh, const label sourceIndex, const label nBands);


    //- Destructor
    virtual ~parallelRadSource();


    // Member functions
	
	    //- Return main axis to draw the new quadrature
	    const List<vector> getAxis(label nTheta, label nPhi);

		    //- Store values of emission
            virtual void generate(const List<vector> directions);

            //- Read radiation properties dictionary
            bool read();

        // Access

			//- Acces to specific stored value
            virtual scalar calculate(const word patchID, const label face, const label direction, const label band) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
