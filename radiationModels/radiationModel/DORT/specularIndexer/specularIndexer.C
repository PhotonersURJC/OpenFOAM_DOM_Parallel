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

#include "specularIndexer.H"
#include "fvm.H"
#include "constants.H"

using namespace Foam::constant;
using namespace Foam::constant::mathematical;
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::specularIndexer::specularIndexer
(
	const List<List< label > > facesNeighbours,
	const List<vector> mainDirs,
	const List<scalar> faceMagnitudes,
	const vector faceN,
	const List<vector> quadVertex,
	const List<face> quadFaces,
	const bool secondOrder
)
:
outers_(mainDirs.size(),(short)-1),
spotters_(0),
viewFactors_(0)
{
	List<List< label > > fNs (facesNeighbours);
	List<short> outIndex(0);
	List<vector> outers(0);
	List<short> refIndex(0);
	List<vector> reflecteds(0);
	
	forAll(quadFaces, dirI)
	{
		bool allOutgoing = true;
		bool allIncoming = true;
		forAll(quadFaces[dirI], vertexI)
		{
			if((-faceN & quadVertex[quadFaces[dirI][vertexI]]) < 0.0) 
			{
				allOutgoing = false;
			}
			if((-faceN & quadVertex[quadFaces[dirI][vertexI]]) > 0.0) 
			{
				allIncoming = false;
			}
		}
		if(allOutgoing)
		{
			outers.setSize(outers.size()+1);
			outers[outers.size()-1]=mainDirs[dirI];
			outIndex.setSize(outIndex.size()+1);
			outIndex[outIndex.size()-1]=dirI;
		}
		else if(allIncoming)
		{
			reflecteds.setSize(reflecteds.size()+1);
			reflecteds[reflecteds.size()-1]=mainDirs[dirI]-2*faceN*(faceN&mainDirs[dirI]);
			refIndex.setSize(refIndex.size()+1);
			refIndex[refIndex.size()-1]=dirI;
			//incoming ordinates are eliminated from neighbours
			forAll(fNs[dirI], neiI) //goes over incoming DO neighbours
			{
				List<label> newNeighbours(0);
				forAll(fNs[fNs[dirI][neiI]], secNeiI) //eliminated from each neighbour
				{
					if(fNs[fNs[dirI][neiI]][secNeiI] != dirI)
					{
						newNeighbours.setSize(newNeighbours.size()+1);
						newNeighbours[newNeighbours.size()-1] = fNs[fNs[dirI][neiI]][secNeiI];
					}
				}
				fNs[fNs[dirI][neiI]] = newNeighbours;
			}
		}
		else
		{
			forAll(fNs[dirI], neiI) //goes over in plane DO neighbours
			{
				List<label> newNeighbours(0);
				forAll(fNs[fNs[dirI][neiI]], secNeiI) //eliminated from each neighbour
				{
					if(fNs[fNs[dirI][neiI]][secNeiI] != dirI)
					{
						newNeighbours.setSize(newNeighbours.size()+1);
						newNeighbours[newNeighbours.size()-1] = fNs[fNs[dirI][neiI]][secNeiI];
					}
				}
				fNs[fNs[dirI][neiI]] = newNeighbours;
			}
		}
	}
	List<short> firstSpot(reflecteds.size(), (short)(-1));
	List<short> secondSpot(reflecteds.size(), (short)(-1));
	//secondSpot is a neighbour ordinate that allows equal both energy and flux
	//Boolean secondOrder determines if a second spotter is used or not.
	//User is chosing between:
	//  -First Order: energy and direction, not flux
	//  -Second Order: energy and flux, not direction
	forAll(reflecteds, refI)
	{
		scalar maxCos = 0.0;
		label maxIndex = 0;
		forAll(outers, outI)
		{
			scalar currentCos = reflecteds[refI] & outers[outI];
			if(currentCos > maxCos)
			{
				maxCos=currentCos;
				maxIndex=outI;
			}
		}
		maxCos = -1.0;
		firstSpot[refI] = outIndex[maxIndex];
		if(secondOrder)
		{
			List<label> neighbours = fNs[outIndex[maxIndex]];
			vector chosenNormal = outers[maxIndex] ^ reflecteds[refI];
			forAll(neighbours, neiI)
			{
				vector currentNormal = reflecteds[refI] ^ mainDirs[neighbours[neiI]];
				currentNormal = currentNormal / mag(currentNormal);
				scalar currentCos = (currentNormal) & chosenNormal;
				if(currentCos > maxCos)
				{
					maxCos=currentCos;
					maxIndex=neiI;
				}
			}
			secondSpot[refI] = neighbours[maxIndex];
		}
	}
	//up to now, two list of index has been generated, 
	//to pair each incoming DO with two outgoings DO
	forAll(firstSpot, spotI)
	{
		scalar omega1 = faceMagnitudes[firstSpot[spotI]];
		scalar omegaR = faceMagnitudes[refIndex[spotI]];
		scalar f1 = omegaR/omega1;
		scalar f2 = 0.0;
		if(secondOrder)
		{
			scalar omega2 = faceMagnitudes[secondSpot[spotI]];
			scalar cos1 = faceN & mainDirs[firstSpot[spotI]];
			scalar cosR = faceN & mainDirs[refIndex[spotI]];
			f1 = omegaR * (1-cosR) / (omega1 * (1-cos1));
			f2 = (omegaR - f1*omega1)/omega2;
		}
		//Intensity in direction i due to this reflected DO is fi*IR
		//filling lists with first DO values
		if(outers_[firstSpot[spotI]]==-1)
		{
			outers_[firstSpot[spotI]]=spotters_.size();
			spotters_.setSize(spotters_.size()+1);
			List<short> newSpot(1, (short)refIndex[spotI]);
			spotters_[spotters_.size()-1]=newSpot;
			viewFactors_.setSize(viewFactors_.size()+1);
			List<scalar> newFactor(1, f1);
			viewFactors_[viewFactors_.size()-1]=newFactor;
		}
		else
		{
			spotters_[outers_[firstSpot[spotI]]].setSize(spotters_[outers_[firstSpot[spotI]]].size()+1);
			spotters_[outers_[firstSpot[spotI]]][spotters_[outers_[firstSpot[spotI]]].size()-1]=refIndex[spotI];
			viewFactors_[outers_[firstSpot[spotI]]].setSize(viewFactors_[outers_[firstSpot[spotI]]].size()+1);
			viewFactors_[outers_[firstSpot[spotI]]][viewFactors_[outers_[firstSpot[spotI]]].size()-1]=f1;
		}
		//filling lists with second DO values
		if(secondOrder)
		{
			if(outers_[secondSpot[spotI]]==-1)
			{
				outers_[secondSpot[spotI]]=spotters_.size();
				spotters_.setSize(spotters_.size()+1);
				List<short> newSpot(1, (short)refIndex[spotI]);
				spotters_[spotters_.size()-1]=newSpot;
				viewFactors_.setSize(viewFactors_.size()+1);
				List<scalar> newFactor(1, f2);
				viewFactors_[viewFactors_.size()-1]=newFactor;
			}
			else
			{
				spotters_[outers_[secondSpot[spotI]]].setSize(spotters_[outers_[secondSpot[spotI]]].size()+1);
				spotters_[outers_[secondSpot[spotI]]][spotters_[outers_[secondSpot[spotI]]].size()-1]=refIndex[spotI];
				viewFactors_[outers_[secondSpot[spotI]]].setSize(viewFactors_[outers_[secondSpot[spotI]]].size()+1);
				viewFactors_[outers_[secondSpot[spotI]]][viewFactors_[outers_[secondSpot[spotI]]].size()-1]=f2;
			}
		}
	}
}

Foam::scalarField Foam::radiation::specularIndexer::getSpec
(
	const Foam::scalarField originalI
) const
{
	scalarField reflectedI(originalI.size(),0.0);
	forAll(reflectedI, dirI)
	{
		if(outers_[dirI]>-1)
		{
			forAll(spotters_[outers_[dirI]], spotI)
			{
				scalar currentI = originalI[spotters_[outers_[dirI]][spotI]];
				currentI = currentI*viewFactors_[outers_[dirI]][spotI];
				reflectedI[dirI] += currentI;
			}
		}
	}
	return reflectedI;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::specularIndexer::~specularIndexer()
{}

// ************************************************************************* //
