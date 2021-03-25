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

#include "discreteOrdinate.H"
#include "fvm.H"
#include "DORT.H"
#include "quadrature.H"
#include "constants.H"

using namespace Foam::constant;


const Foam::word
Foam::radiation::discreteOrdinate::intensityPrefix("ILambda");


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::discreteOrdinate::discreteOrdinate
(
    const DORT& dort,
    const quadrature& quad,
    const fvMesh& mesh,
    const scalar phi,
    const scalar theta,
    const scalar deltaPhi,
    const scalar deltaTheta,
    const label nLambda,
    const label rayId,
    const List<vector> mainAxis
)
:
    dort_(dort),
    quad_(quad),
    mesh_(mesh),
    d_(vector::zero),
    dAve_(vector::zero),
    theta_(theta),
    phi_(phi),
    omega_(0.0),
    nLambda_(nLambda),
    ILambda_(nLambda),
    myRayId_(rayId)
{

    scalar sinTheta = Foam::sin(theta);
    scalar cosTheta = Foam::cos(theta);
    scalar sinPhi = Foam::sin(phi);
    scalar cosPhi = Foam::cos(phi);

    omega_ = 2.0*sinTheta*Foam::sin(deltaTheta/2.0)*deltaPhi;
    if(omega_<0) { omega_ = -omega_;}
    d_ = vector(sinTheta*sinPhi*mainAxis[0])+ (sinTheta*cosPhi*mainAxis[1]) + (cosTheta*mainAxis[2]);
    dAve_ = omega_ * d_;
    autoPtr<volScalarField> IDefaultPtr;
    forAll(ILambda_, lambdaI)
    {
        IOobject IHeader
        (
            intensityPrefix + "_" + name(quad.modelNumber()) + "_" + name(rayId) + "_" +  name(lambdaI),
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );
		
        // check if field exists and can be read
        if (IHeader.typeHeaderOk<volScalarField>(true))
        {
            ILambda_.set
            (
                lambdaI,
                new volScalarField(IHeader, mesh_)
            );
        }
        else
        {
            // Demand driven load the IBand field
            if (!IDefaultPtr.valid())
            {
                IDefaultPtr.reset
                (
                    new volScalarField
                    (
                        IOobject
                        (
                            "IBand_" + name(lambdaI),
                            mesh_.time().timeName(),
                            mesh_,
                            IOobject::MUST_READ,
                            IOobject::NO_WRITE
                        ),
                        mesh_
                    )
                );
            }
			
            // Reset the MUST_READ flag
            IOobject noReadHeader(IHeader);
            noReadHeader.readOpt() = IOobject::NO_READ;

            ILambda_.set
            (
                lambdaI,
                new volScalarField(noReadHeader, IDefaultPtr())
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::discreteOrdinate::~discreteOrdinate()
{}



// ************************************************************************* //
