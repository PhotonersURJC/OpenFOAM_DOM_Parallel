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
    Foam::radiation::isotropicPhaseFunc

Description
    Isotropic phase function

SourceFiles
    isotropicPhaseFunc.C

\*---------------------------------------------------------------------------*/

#ifndef radiationIsotropicPhaseFunc_H
#define radiationIsotropicPhaseFunc_H

#include "phaseFuncModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                  Class isotropicPhaseFunc Declaration
\*---------------------------------------------------------------------------*/

class isotropicPhaseFunc
:
    public phaseFuncModel
{
    // Private data

        //- Scattering isotropic factor
		scalar factor_;

public:

    //- Runtime type information
    TypeName("isotropicPhaseFunction");


    // Constructors

        //- Construct from components
        isotropicPhaseFunc(const dictionary& dict, const fvMesh& mesh);


    //- Destructor
    virtual ~isotropicPhaseFunc();

	    //Scattering Matrix Generation
	    void MatrixGen(const List<vector> directions, const label model);
		
		scalar sctrFactor(const label model, const label band, const label directionFrom, const label directionTo) const;
				
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
