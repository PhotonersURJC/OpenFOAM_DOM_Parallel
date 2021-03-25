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
//#include "MeshedSurface.H"
#include "PrimitivePatch.H"

using namespace Foam::constant;
using namespace Foam::constant::mathematical;
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::specularIndexer::specularIndexer
(
    const Foam::fvPatch& patch,
    const Foam::List<Foam::vector> mainDirs,
	const Foam::List<Foam::vector> vertexList,
	const Foam::List<Foam::face> faces,
	const label nLambda
)
:
    patch_(patch),
    indexes_(patch_.size(),(short)(-1)),
	storedIntensities_(0),
	interpolators_(0)
{
	//almacena una specular matrix por cada direccion
	const vectorField n(patch_.nf());
	for(label lambdaI = 0; lambdaI < nLambda; lambdaI++)
	{
		PtrList<scalarField> bandIntensities(mainDirs.size());
		forAll(bandIntensities, dirI)
		{
			bandIntensities.set(dirI, new scalarField(patch_.size(),0.0));
		}
		storedIntensities_.append(bandIntensities);
	}
	pointField vertex = (pointField)(vertexList);
	List<face> reflectedFaces(faces.size());
	forAll(faces, faceI)
	{
		face currentFace(faces[faceI].size());
		forAll(faces[faceI],vertexI)
		{
			currentFace[currentFace.size()-1-vertexI]=faces[faceI][vertexI];
		}
		reflectedFaces[faceI]=currentFace;
	}
	/*****/
	ofstream facesFile;
	facesFile.open ("casoQuad/constant/polyMesh/faces");
	
	ofstream refFacesFile;
	refFacesFile.open ("casoRef/constant/polyMesh/faces");
	
	facesFile<<"FoamFile \n { \n version     2.0; \n format      ascii; \n ";
	facesFile<<"class       faceList; \n object      faces; \n } \n  \n "<<faces.size()<<" \n ( \n ";
	
	refFacesFile<<"FoamFile \n { \n version     2.0; \n format      ascii; \n ";
	refFacesFile<<"class       faceList; \n object      faces; \n } \n  \n "<<faces.size()<<" \n ( \n ";
	
	forAll(faces, faceI)
	{
		int faceVertex = faces[faceI].size();
		facesFile << faceVertex << "(" << faces[faceI][0];
		refFacesFile << faceVertex << "(" << reflectedFaces[faceI][0];
		for(int i = 1; i < faceVertex; i++)
		{
			facesFile<<" "<<faces[faceI][i];
			refFacesFile<<" "<<reflectedFaces[faceI][i];
		}
		facesFile << ") \n ";
		refFacesFile << ") \n ";
	}
	facesFile<<")";
	refFacesFile<<")";
	facesFile.close();
	refFacesFile.close();
	
	ofstream pointsFile;
	pointsFile.open ("casoQuad/constant/polyMesh/points");
	
	pointsFile<<"FoamFile \n { \n version     2.0; \n format      ascii; \n ";
	pointsFile<<"class       vectorField; \n object      points; \n } \n  \n "<<vertex.size()<<" \n ( \n ";
	
	forAll(vertex, vertI)
	{
		pointsFile<<"(";
		pointsFile<<vertex[vertI][0];
		pointsFile<<"  ";
		pointsFile<<vertex[vertI][1];
		pointsFile<<"  ";
		pointsFile<<vertex[vertI][2];
		pointsFile<<") \n ";
	}
	pointsFile<<")";
	pointsFile.close();
	/*********/
	forAll(patch_, faceI)
	{
		if(indexes_[faceI]==-1)
		{
			vector currentN = n[faceI];
			interpolators_.setSize(interpolators_.size()+1);
			pointField reflectedVertex(vertex.size());
			forAll(vertex, vertI)
			{
				reflectedVertex[vertI]=vertexList[vertI]-2*(currentN & vertexList[vertI])*currentN;
			}
			/*******/
			if(currentN==n[2403])
			{
				ofstream refPointsFile;
				refPointsFile.open ("casoRef/constant/polyMesh/points");
	
				refPointsFile<<"FoamFile \n { \n version     2.0; \n format      ascii; \n ";
				refPointsFile<<"class       vectorField; \n object      points; \n } \n  \n "<<reflectedVertex.size()<<" \n ( \n ";
	
				forAll(reflectedVertex, vertI)
				{
					refPointsFile<<"(";
					refPointsFile<<reflectedVertex[vertI][0];
					refPointsFile<<"  ";
					refPointsFile<<reflectedVertex[vertI][1];
					refPointsFile<<"  ";
					refPointsFile<<reflectedVertex[vertI][2];
					refPointsFile<<") \n ";
				}
				refPointsFile<<")";
				refPointsFile.close();
			}
			/*******/
			PrimitivePatch<face,List,pointField,point> quadraturePatch(faces, vertex);
			PrimitivePatch<face,List,pointField,point> reflectedPatch(reflectedFaces, reflectedVertex);		
			interpolators_.set(interpolators_.size()-1,new specularMatrix(quadraturePatch,reflectedPatch,mainDirs, currentN));
			indexes_[faceI]=interpolators_.size()-1;
			for(label i=faceI+1;i<patch_.size();i++)
			{
				if(n[i]==currentN) //may be introduced a relaxed criteria
				{
					indexes_[i]=interpolators_.size()-1;
				}
			}
		}
	}
}

// * * * * * * * * * * * * * * Member functions * * * * * * * * * * * * * * * //

void Foam::radiation::specularIndexer::updateSpecular
(
	const Foam::PtrList<Foam::scalarField> patchIntensity,
	const Foam::label lambdaI
)
{
	forAll(patch_, faceId)
	{
		//generates a scalarField for every face
		scalarField faceIntensities(patchIntensity.size());
		forAll(faceIntensities, dirI)
		{
			if(!interpolators_[indexes_[faceId]].dirOut(dirI)) 
			{ //only if direction is incoming
				faceIntensities[dirI] = patchIntensity[dirI][faceId];
			}
		}
		scalarField reflectedIntensities = 
			interpolators_[indexes_[faceId]].getSpec(faceIntensities);
		/******/
		if(faceId==2403)
		{
			mkDir("casoQuad/" + patch_.boundaryMesh().mesh().time().timeName());
			ofstream intensitiesFile;
			intensitiesFile.open ("casoRef/" + patch_.boundaryMesh().mesh().time().timeName() + "/I");
	
			intensitiesFile<<"FoamFile \n { \n version     2.0; \n format      ascii; \n ";
			intensitiesFile<<"class       volScalarField; \n object      I; \n } \n  \n ";
			intensitiesFile<<"dimensions [0 0 0 0 0 0 0]; \n internalField uniform 0.0; \n";
			intensitiesFile<<"boundaryField \n { \n     sphere \n { \n type   calculated; \n";
			intensitiesFile<<"value   nonuniform List<scalar> \n"<<faceIntensities.size()<<"\n ( \n";
			forAll(faceIntensities, faceI)
			{
				intensitiesFile<<faceIntensities[faceI];
				intensitiesFile<<" \n ";
			}
			intensitiesFile<<"); \n } \n }";
			intensitiesFile.close();
			
			mkDir("casoRef/" + patch_.boundaryMesh().mesh().time().timeName());
			ofstream reflexionFile;
			reflexionFile.open ("casoQuad/" + patch_.boundaryMesh().mesh().time().timeName() + "/I");
	
			reflexionFile<<"FoamFile \n { \n version     2.0; \n format      ascii; \n ";
			reflexionFile<<"class       volScalarField; \n object      I; \n } \n  \n ";
			reflexionFile<<"dimensions [0 0 0 0 0 0 0]; \n internalField uniform 0.0; \n";
			reflexionFile<<"boundaryField \n { \n     sphere \n { \n type   calculated; \n";
			reflexionFile<<"value   nonuniform List<scalar> \n"<<reflectedIntensities.size()<<"\n ( \n";
			forAll(reflectedIntensities, faceI)
			{
				reflexionFile<<reflectedIntensities[faceI];
				reflexionFile<<" \n ";
			}
			reflexionFile<<"); \n } \n }";
			reflexionFile.close();
		}
		/********/
		forAll(reflectedIntensities, dirI)
		{
			if(interpolators_[indexes_[faceId]].dirOut(dirI))
			{ //only if outgoing direction
				storedIntensities_[lambdaI][dirI][faceId] = reflectedIntensities[dirI];
			}
		}
	}
}

Foam::scalarField Foam::radiation::specularIndexer::getSpec
(
	const Foam::label rayId, const Foam::label lambdaId
) const
{
	return storedIntensities_[lambdaId][rayId];
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::specularIndexer::~specularIndexer()
{}

// ************************************************************************* //
