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
    Foam::radiation::absorptionDO

Description
    Discrete Ordinate in a region with homogeneous optical properties
	Scattering matrix includes scattering coefficent, to optimize iterations

SourceFiles
    absorptionDO.C

\*---------------------------------------------------------------------------*/

#ifndef absorptionDO_H
#define absorptionDO_H

#include "discreteOrdinate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                    Class absorptionDO Declaration
\*---------------------------------------------------------------------------*/

class absorptionDO
:
    public discreteOrdinate
{

	//- absorption coefficient
    dimensionedScalar k_;

public:

    // Constructors

        //- Construct form components
        absorptionDO
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
	    const List<vector> mainAxis,
            const scalar k
        );


    //- Destructor
    ~absorptionDO();


    // Member functions

        // Edit

            //- Update radiative intensity on i direction
            scalar correct();

	    void setBandScatter (const label lambdaI, const volScalarField value);
		
		void changeBand (const label lambdaI);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
