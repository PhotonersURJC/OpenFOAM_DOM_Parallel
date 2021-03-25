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
    Foam::radiation::HenyeyGreenteinsPhaseFunc

Description
    Henyey-Greenteins phase function

SourceFiles
    HenyeyGreenteinsPhaseFunc.C

\*---------------------------------------------------------------------------*/

#ifndef radiationHenyeyGreenteinsPhaseFunc_H
#define radiationHenyeyGreenteinsPhaseFunc_H

#include "phaseFuncModel.H"
#include "sctrMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{
	
	class sctrMatrix;

/*---------------------------------------------------------------------------*\
                  Class HenyeyGreenteinsPhaseFunc Declaration
\*---------------------------------------------------------------------------*/

class HenyeyGreenteinsPhaseFunc
:
    public phaseFuncModel
{
    // Private data


        //- PhaseFunc model dictionary
        dictionary coeffsDict_;

		dimensionedScalar g_;
		
		//3D matrix to store multiband scattering factors;
		List< List < List < scalar > > > sctrMatrix_;
		
		label matrixSize_;

public:

    //- Runtime type information
    TypeName("HenyeyGreenteinsPhaseFunction");


    // Constructors

        //- Construct from components
        HenyeyGreenteinsPhaseFunc(const dictionary& dict, const fvMesh& mesh);


    //- Destructor
    virtual ~HenyeyGreenteinsPhaseFunc();

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
