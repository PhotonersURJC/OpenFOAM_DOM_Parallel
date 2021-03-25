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

Class
    Foam::radiation::quadrature

Description

\*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*\
	Scheme to distribute ordinates, with local sphere coordinates system,
	adapted to light sources (detailed in each kind of source).
	
	Ordinates are distributed keeping a constant increment in theta and phi
	angles, resulting in solid angles of different size (omega).
	
	A quadrature is composed of discrete ordinates and a specular guide.

SourceFiles
    quadrature.C

\*---------------------------------------------------------------------------*/

#ifndef radiationModelquadrature_H
#define radiationModelquadrature_H

#include "discreteOrdinate.H"
#include "specularGuide.H"
#include "DORT.H"
#include "radSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                           Class quadrature Declaration
\*---------------------------------------------------------------------------*/

class quadrature
{
    // Private data

		//- Boss Model
		const DORT& DORT_;

		//- Model Number
		label modelNumber_;

		//- Reference to the mesh database
        const fvMesh& mesh_;
		
		//-Is quad being solved serial?
		bool serialQuad_;
		
		//-If so, band being solved
		label actualBand_;

        //- Incident radiation  [W/m2]
        PtrList<volScalarField> Glambda_;

        //- Total radiative heat flux [W/m2]
        PtrList<volScalarField> QrLambda_;

		//- Live radiative heat flux [W/m2]
		PtrList<volScalarField> futQrLambda_;
		
		//- Radiative heat flux of converged ordinates[W/m2]
		PtrList<volScalarField> convQrLambda_;
		
		//- Radiative heat flux in last ordinate[W/m2]
		PtrList<volScalarField> preCalcQrLambda_;

        //- Incident radiative heat flux [W/m2]
        PtrList<volScalarField> QinLambda_;

		//- Live incident radiative heat flux [W/m2]
		PtrList<volScalarField> futQinLambda_;
		
		//- Incident radiative heat flux of converged ordinates[W/m2]
		PtrList<volScalarField> convQinLambda_;
		
		//- Incident radiative heat flux in last ordinate[W/m2]
		PtrList<volScalarField> preCalcQinLambda_;

        //- Number of solid angles in theta
        label nTheta_;

        //- Number of solid angles in phi
        label nPhi_ ;

        //- Total number of rays (1 per direction)
        label nRay_;

        //- Number of wavelength bands
        label nLambda_;

        //- List of pointers to radiative intensity rays
        PtrList<discreteOrdinate> IRay_;

        //- Convergence criterion
        scalar convergence_;

        //- Maximum number of iterations
        scalar maxIter_;

        //- Maximum omega weight
        scalar omegaMax_;

		//- Scattering Matrix
		List<scalar> sctrMatrix;

		//- Determine if sigma is bigger than 0.0 in any point
		bool scatter_;

		//- Cone or Ray direction
		vector CORdirection_;

		//- Cone opening. =0 if Ray
		scalar coneOpening_;    

		//- Way to calculate source term
		radSource& radSource_;
		
		//- List of discrete ordinates directions
		List<vector> directions_;
		
		//- List of specularIndexer positions
		List<label> specPositions_;
		
		//- List of specular faces order(true == secondOrder)
		List<bool> specOrders_;
		
		PtrList<specularGuide> specGuiders_;
		
		label lastSpecked_;

    // Private Member Functions

        

        //- Disallow default bitwise copy construct
        quadrature(const quadrature&);
		
		//- Create band fields
		void createPtrLists(PtrList<volScalarField>& currentField, word Name);
		
		//-Cleaning of temporal fluxes
		void cleanFutures();
		void cleanPreCalculateds();
		void cleanConvergeds();

        //- Disallow default bitwise assignment
        void operator=(const quadrature&);


public:

    // Constructors

        //- Construct from components
        quadrature(const DORT& DORT, const fvMesh& mesh, const label modelNumber,  radSource& radSource, const bool serialQuad);

    //- Destructor
    virtual ~quadrature();


    // Member functions

        // Edit

			//- Initialise
			void initialise();
			
			
			//- Generate specular guiders
			void specularGuiding(const List<vector> mainAxis, const label cut, const scalar opening);
		
			//- Update specular reflection at the begining of each iteration
			void updateSpecular(const label lambdaI);
		
            //- Solve radiation equation(s)
            bool calculate();
			
			//- Solve radiation equation(s) in a single band
			bool calculate(const label band);

            //- Update G and calculate total heat flux on boundary
            void updateG();
			
			//-Same for single band
			void updateBandG(const label lambdaI);

            //- Update inScatter for each Ray in every band
            void updateScatter();
			
			//- Update inScatter for each Ray in single band
            void updateBandScatter(const label lambdaI);

            //- Set the rayId and lambdaId from by decomposing an intensity
            //  field name
            void setRayIdLambdaId
            (
                const word& name,
                label& rayId,
                label& lambdaId
            ) const;
			
			void changeBand(const label lambdaI);
			
			void getBoundaryReport
			(
				const bool specUse,
				const bool specOrder,
				const label patchId,
				const label lambdaId
			);


        // Access
		
			//- Specular reflection (including intensity and omega) in given direction
			scalarField getSpec(label patchId, label rayId, label lambdaId);

            //- Number of wavelengths
            inline label modelNumber() const;

            //- Ray intensity for rayI
            inline const discreteOrdinate& IRay(const label rayI) const;

            //- Ray intensity for rayI and lambda bandwidth
            inline const volScalarField& IRayLambda
            (
                const label rayI,
                const label lambdaI
            ) const;

            //- Number of angles in theta
            inline label nTheta() const;

            //- Number of angles in phi
            inline label nPhi() const;

            //- Number of rays
            inline label nRay() const;

            //- Number of wavelengths
            inline label nLambda() const;

            //- Const access to incident radiation field
            inline const volScalarField& Glambda(const label lambdaI) const;

            //- Const access to total radiative heat flux field
            inline const volScalarField& QrLambda(const label lambdaI) const;

            //- Const access to incident radiative heat flux field
            inline const volScalarField& QinLambda(const label lambdaI) const;

            //- Return non-const access to the boundary heat flux
            inline volScalarField& futQrLambda(const label lambdaI);

            //- Return non-const access to the boundary incident heat flux
            inline volScalarField& futQinLambda(const label lambdaI);

            //- Return omegaMax
            inline scalar omegaMax() const;

			//- Return source
			inline scalar sourceFlux
			(const word patchID, const label face, const label rayId, const label band) const;
			
			//- Const access to list of discrete Ordinates directions
			inline const List<vector> directions() const;
			
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "quadratureI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
