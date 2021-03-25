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

#include "greyConstExt.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(greyConstExt, 0);

        addToRunTimeSelectionTable
        (
            radExtintionModel,
            greyConstExt,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::greyConstExt::greyConstExt
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    radExtintionModel(dict, mesh),
    coeffsDict_(dict.subDict(typeName + "Coeffs")),
    kappa_(coeffsDict_.lookup("absorptivity")),
	sigma_(coeffsDict_.lookup("scatter"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::greyConstExt::~greyConstExt()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiation::greyConstExt::kappaCont(const label bandI) const
{
    tmp<volScalarField> ta
    (
        new volScalarField
        (
            IOobject
            (
                "kappa",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            kappa_
        )
    );

    return ta;
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::greyConstExt::sigmaCont(const label bandI) const
{
    tmp<volScalarField> ts
    (
        new volScalarField
        (
            IOobject
            (
                "sigma",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            sigma_
        )
    );

    return ts;
}


// ************************************************************************* //
