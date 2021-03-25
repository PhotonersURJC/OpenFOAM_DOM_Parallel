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

#include "constConeRadSource.H"
#include "phaseFuncModel.H"
#include "radExtintionModel.H"
#include "constants.H"
#include "fvm.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam::constant;
using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(constConeRadSource, 0);
        addToSourceRunTimeSelectionTables(constConeRadSource);
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::radiation::constConeRadSource::initialise()
{
    readIfPresent("Direction", direction_);
    direction_=direction_/mag(direction_);
    readIfPresent("Opening", opening_);
	opening_ = (opening_*pi)/360.0;
	sourceValues_.setSize(patchesID_.size());
	patchNormals_.setSize(patchesID_.size());
	label patchIndex = 0;
    while (patchIndex < patchesID_.size())
    {
        forAll(mesh_.boundary(), patchI)
        {
			const polyPatch currentPatch(mesh_.boundary()[patchI].patch());
			if(currentPatch.name()==patchesID_[patchIndex])
            {
				vector patchDir(0,0,0);
				forAll(currentPatch, faceI)
				{
					patchDir += currentPatch.faceNormals()[faceI];
				}
				patchNormals_[patchIndex]=patchDir/mag(patchDir);
			}
		}
		patchIndex++;
	}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::constConeRadSource::constConeRadSource(const fvMesh& mesh, const label sourceIndex, const label nBands)
:
    radSource(typeName, mesh, sourceIndex, nBands),
	sourceValues_(0),
	patchNormals_(0)
{
    initialise();
}


Foam::radiation::constConeRadSource::constConeRadSource
(
    const dictionary& dict,
    const fvMesh& mesh,
    const label sourceIndex,
	const label nBands
)
:
    radSource(typeName, dict, mesh, sourceIndex, nBands),
	sourceValues_(0),
	patchNormals_(0)
{
    initialise();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::constConeRadSource::~constConeRadSource()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::radiation::constConeRadSource::read()
{
    if (radSource::read())
    {


        return true;
    }
    else
    {
        return false;
    }
}

const Foam::List<Foam::vector> Foam::radiation::constConeRadSource::getAxis
(Foam::label nTheta, Foam::label nPhi)
{
    List<vector> mainAxis(3,vector());
    vector Vx(1,0,0); vector Vy(0,1,0); vector Vz(0,0,1);
    if((Vz & direction_) == -1.0)
    { Vx = -Vx; Vy = -Vy; Vz = -Vz; }
    else if ((Vz & direction_) < 1.0)
    {
	vector VCross = Vz ^ direction_;
	VCross = VCross / mag(VCross);
	vector VSum = Vz + direction_;
	VSum = VSum / mag(VSum);
	vector VNormal = VCross ^ VSum;
	Vx = (VCross*VCross.x()) + (VSum*VSum.x()) - (VNormal*VNormal.x());
	Vy = (VCross*VCross.y()) + (VSum*VSum.y()) - (VNormal*VNormal.y());
	Vz = direction_;
    }
    mainAxis[0] = Vx; mainAxis[1] = Vy; mainAxis[2] = Vz;
    return mainAxis;
}

Foam::scalar Foam::radiation::constConeRadSource::calculate
(
	const word patchID, const label face, const label direction, const label band
) const
{
	if(patchesID_.size()>0)
	{
		forAll(patchesID_, patch)
		{
			if(patchID == patchesID_[patch])
			{
				return sourceValues_[patch].getValue(face, direction, band);
			}
		}
	}
	return 0.0;
}

void Foam::radiation::constConeRadSource::generate(const List<vector> directions)
{
	if(patchesID_.size()==0){return;}
    label patchIndex = 0;
	scalar dirCorrection = Foam::constant::mathematical::pi*Foam::sin(opening_)*Foam::sin(opening_)*(direction_&patchNormals_[patchIndex]);
    while (patchIndex < patchesID_.size())
    {
        forAll(mesh_.boundary(), patchI)
        {
			const polyPatch currentPatch(mesh_.boundary()[patchI].patch());
            if(currentPatch.name()==patchesID_[patchIndex])
            {
				sourceValues_.set(patchIndex,new booleanSourceStorage(nBands_,currentPatch.size(),directions.size()));
				List<scalar> modQs = patchesQ_[patchIndex];
				forAll(modQs, currentQ)
				{
					modQs[currentQ]=-modQs[currentQ]/dirCorrection;
				}
				sourceValues_[patchIndex].setQs(modQs);
				forAll(currentPatch, faceI)
				{
					label dirIndex = 0;
					while(dirIndex<directions.size())
					{
						if((directions[dirIndex]&direction_)>Foam::cos(opening_))
						{
							sourceValues_[patchIndex].setBools(faceI,dirIndex,true);
							dirIndex++;
						}
						else
						{
							dirIndex=directions.size();
						}
					}
				}
			}
        }
        patchIndex++;
    }
}




// ************************************************************************* //
