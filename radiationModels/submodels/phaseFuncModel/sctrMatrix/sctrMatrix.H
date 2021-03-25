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
    Foam::radiation::sctrMatrix

Description
    Struct to store scattering factors between a pair of directions.

SourceFiles
    sctrMatrix.C

\*---------------------------------------------------------------------------*/

#ifndef sctrMatrix_H
#define sctrMatrix_H

#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{


/*---------------------------------------------------------------------------*\
                    Class sctrMatrix Declaration
\*---------------------------------------------------------------------------*/

class sctrMatrix
{

        //- List of specular index in patch faces
        List<scalar> sctrFactors_;
		
		label directions_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        specularIndexer(const sctrMatrix&);

        //- Disallow default bitwise assignment
        void operator=(const sctrMatrix&);


public:

    // Constructors

        //- Construct from components
        sctrMatrix
        (
            const label directions
        );
		


    //- Destructor
    ~sctrMatrix();


    // Member functions

		//- Return non-const access to sctrFactor from dir i to dir j
			scalar& Item(const label i, const label j);
			
		//- Return const access to sctrFactor from dir i to dir j
		scalar constItem(const label i, const label j) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
