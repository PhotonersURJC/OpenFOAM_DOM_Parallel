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

#include "newRadiationModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::autoPtr<Foam::radiation::newRadiationModel>
Foam::radiation::newRadiationModel::New
(
    const volScalarField& T
)
{
    IOobject radIO
    (
        "radiationProperties",
        T.time().constant(),
        T.mesh(),
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE,
        false
    );
    word modelType("none");
    if (radIO.good())
    {
        IOdictionary(radIO).lookup("radiationModel") >> modelType;
    }
    else
    {
        Info<< "Radiation model not active: radiationProperties not found"
            << endl;
    }
    Info<< "Selecting radiationModel " << modelType << endl;
    TConstructorTable::iterator cstrIter =
        TConstructorTablePtr_->find(modelType);

    if (cstrIter == TConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "noTradiationModel::New(const fvMesh&)"
        )   << "Unknown noTradiationModel type "
            << modelType << nl << nl
            << "Valid noTradiationModel types are:" << nl
            << TConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }
    return autoPtr<newRadiationModel>(cstrIter()(T));
}

Foam::autoPtr<Foam::radiation::newRadiationModel>
Foam::radiation::newRadiationModel::New
(
    const dictionary& dict,
    const volScalarField& T
)
{
    const word modelType(dict.lookup("radiationModel"));

    Info<< "Selecting radiationModel " << modelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "noTradiationModel::New(const dictionary&, const fvMesh&)"
        )   << "Unknown radiationModel type "
            << modelType << nl << nl
            << "Valid radiationModel types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<newRadiationModel>(cstrIter()(dict, T));
}


// ************************************************************************* //
