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
#include "fvmSup.H"
#include "constants.H"
#include <cstdlib>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(radSource, 0);
        defineRunTimeSelectionTable(radSource, mesh);
        defineRunTimeSelectionTable(radSource, dictionary);
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::IOobject Foam::radiation::radSource::createIOobject
(
    const fvMesh& mesh, const label sourceIndex
) const
{
    IOobject io
    (
        "radSource" + Foam::name(sourceIndex),
        mesh.time().constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );
    return io;
}


void Foam::radiation::radSource::initialise()
{
	dictionary patches(subOrEmptyDict("patchesList"));
	if(nBands_<2)
	{
		forAll(mesh_.boundary(), patchI)
		{
			scalar currentValue = 0.0;
			if (patches.readIfPresent(mesh_.boundary()[patchI].name(),currentValue))//Searching face and flux
			{
				patchesID_.resize(patchesID_.size()+1);
				patchesQ_.resize(patchesID_.size());
				patchesID_[patchesID_.size()-1]=mesh_.boundary()[patchI].name();	//Storing face name
				List<scalar> currentList (1);
				currentList[0] = currentValue;
				patchesQ_[patchesID_.size()-1]=currentList;		//Saving face flux
			}
		}
	}
	else
	{
		forAll(mesh_.boundary(), patchI)
		{
			scalar currentValue = 0.0;
			if(patches.found(mesh_.boundary()[patchI].name()))
			{
				patchesID_.resize(patchesID_.size()+1);
				patchesQ_.resize(patchesID_.size());
				patchesID_[patchesID_.size()-1]=mesh_.boundary()[patchI].name();
				List<scalar> currentList (nBands_);
				if(patches.isDict(mesh_.boundary()[patchI].name()))
				{
					dictionary currentPatch = patches.subOrEmptyDict(mesh_.boundary()[patchI].name());
					wordList bandNames = currentPatch.toc();
					for(label lambdaI = 0; lambdaI < nBands_; lambdaI++)
					{
						if(currentPatch.readIfPresent(bandNames[lambdaI],currentValue))
						{
							currentList[lambdaI]=currentValue;
						}
					}
				}
				patchesQ_[patchesID_.size()-1]=currentList;
			}
		}
	}

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::radSource::radSource(const fvMesh& mesh, const label sourceIndex, const label nBands)
:
    IOdictionary
    (
        IOobject
        (
            "radSource" + Foam::name(sourceIndex),
            mesh.time().constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    time_(mesh.time()),
    index_(sourceIndex),
	nBands_(nBands),
    patchesID_(0),
    patchesQ_(0),
    direction_(),
    opening_(0.0)
	
{
	initialise();
}


Foam::radiation::radSource::radSource
(
    const word& type,
    const fvMesh& mesh,
    const label sourceIndex,
	const label nBands
)
:
    IOdictionary(createIOobject(mesh, sourceIndex)),
    mesh_(mesh),
    time_(mesh.time()),
    index_(sourceIndex),
	nBands_(nBands),
    patchesID_(0),
    patchesQ_(0),
    direction_(),
    opening_(0.0)
	
{
    initialise();
}


Foam::radiation::radSource::radSource
(
    const word& type,
    const dictionary& dict,
    const fvMesh& mesh,
    const label sourceIndex,
	const label nBands
)
:
    IOdictionary
    (
        IOobject
        (
            "radSource" + Foam::name(sourceIndex),
            mesh.time().constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dict
    ),
    mesh_(mesh),
    time_(mesh.time()),
    index_(sourceIndex),
	nBands_(nBands),
    patchesID_(0),
    patchesQ_(0),
    direction_(),
    opening_(0.0)
	
{
    initialise();
}


// * * * * * * * * * * * * * * * * Destructor    * * * * * * * * * * * * * * //

Foam::radiation::radSource::~radSource()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::radiation::radSource::read()
{
    return true;
}




// ************************************************************************* //
