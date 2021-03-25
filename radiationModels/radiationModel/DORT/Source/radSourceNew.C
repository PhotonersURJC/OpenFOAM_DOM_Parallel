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

#include "radSource.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::radiation::radSource>
Foam::radiation::radSource::New
(
    const fvMesh& mesh,
    const label sourceIndex,
	const label nBands
)
{
    IOobject radIO
    (
        "radSource" + Foam::name(sourceIndex),
        mesh.time().constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );

    word modelType("none");
    if (radIO.good())
    {
        IOdictionary(radIO).lookup("sourceType") >> modelType;
    }

    meshConstructorTable::iterator cstrIter =
        meshConstructorTablePtr_->find(modelType);

    if (cstrIter == meshConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "radSource::New(const fvMesh&)"
        )   << "Unknown radSource type "
            << modelType << nl << nl
            << "Valid radSource types are:" << nl
            << meshConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<radSource>(cstrIter()(mesh, sourceIndex, nBands));
}


Foam::autoPtr<Foam::radiation::radSource>
Foam::radiation::radSource::New
(
    const dictionary& dict,
    const fvMesh& mesh,
    const label sourceIndex,
	const label nBands
)
{
    const word modelType(dict.lookup("sourceType"));

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "radSource::New(const dictionary&, const fvMesh&)"
        )   << "Unknown radiSource type "
            << modelType << nl << nl
            << "Valid radSource types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<radSource>(cstrIter()(dict, mesh, sourceIndex, nBands));
}


// ************************************************************************* //
