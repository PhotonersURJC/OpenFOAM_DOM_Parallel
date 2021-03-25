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

#include "specularMatrix.H"
#include "fvm.H"
#include "constants.H"

using namespace Foam::constant;
using namespace Foam::constant::mathematical;
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::specularMatrix::specularMatrix
(
	const PrimitivePatch<face,List,pointField,point> quadraturePatch,
	const PrimitivePatch<face,List,pointField,point> reflectedPatch,
	const List<vector> mainDirs,
	const vector faceN
)
:
fromPatch_(reflectedPatch),
toPatch_(quadraturePatch),
interpolator_(fromPatch_,toPatch_),
dirOut_(mainDirs.size(),false)
{
	forAll(mainDirs, dirI)
	{
		if((-faceN & mainDirs[dirI]) > 0.0)
		{
			dirOut_[dirI]=true;
		}
	}
}

Foam::scalarField Foam::radiation::specularMatrix::getSpec
(
	const Foam::scalarField originalI
) const
{
	return interpolator_.faceInterpolate(originalI);
}

bool Foam::radiation::specularMatrix::dirOut
(
	const label dirI
) const
{
	return dirOut_[dirI];
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::specularMatrix::~specularMatrix()
{}

// ************************************************************************* //
