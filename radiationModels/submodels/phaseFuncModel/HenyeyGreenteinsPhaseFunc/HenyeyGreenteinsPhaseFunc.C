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

\*---------------------------------------------------------------------------*/

#include "HenyeyGreenteinsPhaseFunc.H"
#include "addToRunTimeSelectionTable.H"
#include "constants.H"

using namespace Foam::constant;
using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(HenyeyGreenteinsPhaseFunc, 0);

        addToRunTimeSelectionTable
        (
            phaseFuncModel,
            HenyeyGreenteinsPhaseFunc,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::HenyeyGreenteinsPhaseFunc::HenyeyGreenteinsPhaseFunc
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    phaseFuncModel(dict, mesh),
    coeffsDict_(dict.subDict(typeName + "Coeffs")),
    g_(coeffsDict_.lookup("asymmetry")),
	sctrMatrix_(0),
	matrixSize_(0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::HenyeyGreenteinsPhaseFunc::~HenyeyGreenteinsPhaseFunc()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::radiation::HenyeyGreenteinsPhaseFunc::MatrixGen
(const List<vector> directions, const label model)
{
	if(sctrMatrix_.size()<1+model)
	{
		sctrMatrix_.resize(model+1);
	}
	sctrMatrix_[model].resize(directions.size());
	forAll(directions, dir)
	{
		sctrMatrix_[model][dir].resize(directions.size());
	}
    scalar g = g_.value();
    if (g >= 1) // No scattering, all radiation scattered forward
    {
	forAll (directions, Index1)
	{
	    sctrMatrix_[model][Index1][Index1]=1.0;
	}
	return;
    }
    if (g <= -1) //All light scattered backward
    {
	forAll (directions, Index1)
	{
	    scalar minCos = 1;
	    label Index3 = Index1;
	    forAll (directions, Index2)
	    {
		    scalar cos = (directions[Index1]&directions[Index2]);
		    if (cos < minCos)
		    { minCos = cos; Index3 = Index2; }
	    }
	    sctrMatrix_[model][Index3][Index1]=1.0;
	}
	return;
    }
    forAll (directions, Index1)
    {
	forAll (directions, Index2)
	{
	    scalar cos = (directions[Index1]&directions[Index2]);
	    sctrMatrix_[model][Index1][Index2]=(1.0-g*g)/Foam::sqrt(pow(1.0+g*g-2.0*g*cos,3));
	}
    }
    List<scalar> directionWeights(directions.size(),0.0);
    forAll (directions, Index1)
    {
	forAll (directions, Index2)
	{
		directionWeights[Index2]+=sctrMatrix_[model][Index1][Index2];
	}
    }
    forAll (directions, Index1)
    {
	forAll (directions, Index2)
	{
		sctrMatrix_[model][Index1][Index2]=sctrMatrix_[model][Index1][Index2]/directionWeights[Index2];
	}
    }
}

Foam::scalar Foam::radiation::HenyeyGreenteinsPhaseFunc::sctrFactor
(const label model, const label band, const label directionFrom, const label directionTo)
const
{
	return sctrMatrix_[model][directionFrom][directionTo];
}


// ************************************************************************* //
