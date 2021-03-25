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
    Foam::radiation::DORT

Description

    Discrete Ordinates Radiative Transfer. Solves the RTE equation for n
    directions in a participating media

    Available extintion models:
        cromConstExt - homogeneous coefficients for each wavelength band
		greyConstExt - homogeneous coefficient
		
	Available phase function models:
		isotropicPhaseFunc			- isotropic phase function
		HenyeyGreensteinPhaseFunc	- Henyey-Greenstein phase function
		coeffPhaseFunc				- multi-band model: same phase function, different coeffs
		bandPhaseFunc				- multi-band model: a specific phase function in each wavelength band

    i.e. dictionary (radiationProperties in region/constant)
    \verbatim
        DORTCoeffs
        {
            nPhi        4;          // azimuthal angles in PI/2 on X-Y. (from Y to X)
            nTheta      0;          // polar angles in PI/2 (from Z to X-Y plane)
			
            convergence 1e-3;       // convergence criteria for radiation iteration
            maxIter     4;          // iterations between inScatter update
            cacheDiv    true;       // cache the div of the RTE equation.
				//NOTE: Caching div is "only" accurate if the upwind 
				//scheme is used in div(Ji,Ii_h)
			
			//ONLY IN SPECIFIC SOLVER: 
			nBands		1;			// number of bands to solve in serial mode
			nQuad		1;			// number of quadratures to solve in serial mode
        }

        solverFreq   1; // Number of flow iterations per radiation iteration
    \endverbatim
	
    The total number of solid angles is  8*nPhi*nTheta.

    In 1D the direction of the rays is X and -X(nPhi and nTheta are ignored)
    In 2D the direction of the rays is on X-Y plane (only nPhi is considered)
    In 3D (nPhi and nTheta are considered)
	
	A DORT is composed of one or more quadratures, 
	and the same number of radSources

SourceFiles
    DORT.C

\*---------------------------------------------------------------------------*/

#ifndef radiationModelDORT_H
#define radiationModelDORT_H

#include "quadrature.H"
#include "radSource.H"
#include "newRadiationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                           Class DORT Declaration
\*---------------------------------------------------------------------------*/

class DORT
:
    public newRadiationModel
{
    // Private data

	    //- Number of wavelength bands
        label nLambda_;
		
		//- bands limits
		List<Vector2D<scalar>> bandsLimits_;
	
	    //- Total absorption coefficient [1/m]
        volScalarField kappa_;
		
        //- Total scattering coefficient [1/m]
        volScalarField sigma_;		
		
		//- Total Absorbed energy [W/m3]
		volScalarField eA_;
		
        //- Incident radiation  [W/m2]
        PtrList<volScalarField> Glambda_;

        //- Total radiative heat flux [W/m2]
        PtrList<volScalarField> QrLambda_;

        //- Incidet radiative heat flux [W/m2]
        PtrList<volScalarField> QinLambda_;
		
		//- Band Absorbed energy [W/m3]
		PtrList<volScalarField> eAlambda_;
		
		//- List of Flux Fields in specific Solver
		PtrList<volScalarField> modelsG_;
		PtrList<volScalarField> modelsQr_;
		PtrList<volScalarField> modelsQin_;

        //- Number of solid angles in theta
        label nTheta_;

        //- Number of solid angles in phi
        label nPhi_;
		
		//- Number of quadratures in the model
		label nQuad_;
		
		//- Flag of actual band
		label actualBand_;

        //- Wavelength total absorption coefficient [1/m]
        PtrList<volScalarField> kappaLambda_;
		
        //- Wavelength total scattering coefficient [1/m]
        PtrList<volScalarField> sigmaLambda_;
		
		//- Temperature emission [W/m3]
		PtrList<volScalarField> Elambda_;

	//- Sub models list
	PtrList<quadrature> submodels_;

	//- List of sources
	PtrList<radSource> sources_;

        //- Convergence criterion
        scalar convergence_;

        //- Maximum number of iterations
        scalar maxIter_;

        //- Determine if extintion coeffs are constant and homogeneus
        bool constantExt_;
		
		//- Flag of thermal calculation
		bool thermo_;

    // Private Member Functions

        //- Initialise
        void initialise(const dictionary& dict);
		
		//- Read own dictionary (DORTproperties)
		void readOwnDict(const dictionary& dict);

        //- Disallow default bitwise copy construct
        DORT(const DORT&);

        //- Disallow default bitwise assignment
        void operator=(const DORT&);
		
		//- Update gray temperature emission
		void setNoBandEmission ();
		
		//- Update band temperature emission
		void setBandEmission (label lambdaI);


public:

    //- Runtime type information
    TypeName("DORT");


    // Constructors
		
		//- Construct from components (General Solver)
        DORT(const volScalarField& T);
		
		//- Construct from components (DORT Solver)
        DORT(const dictionary& dict, const volScalarField& T);


    //- Destructor
    virtual ~DORT();


    // Member functions

        // Edit

            //- Solve radiation equation(s)
            bool calculate();

			//- Solve specific quadrature
            void changeBand(const label model, const label band);
			
			//- Solve specific quadrature
            bool calculateQuad(const label model);
			
			//- Generate specific quadrature
			void createQuad(const label model);
			
			//- Clear submodels list
			void deleteQuad();
			
			//- Generate grey dictionary for specific band called newDortDict
			void generateDicts (const label band, IOdictionary& dict) const;

            //- Read radiation properties dictionary
            bool read();

	    //- Set the ray owner by name
            void setRayModelIndex
            (
                const word& name,
                label& modelIndex
            ) const;

        // Access

            //- Is extintion constant flag
            bool constantExt() const;
			
			//- Is temperature considered flag
			bool thermo() const;

            //- Number of angles in theta
            inline label nTheta() const;

            //- Number of angles in phi
            inline label nPhi() const;
			
            //- Number of directions in quadratures
            inline label nRay() const;

            //- Number of wavelengths
            inline label nLambda() const;
			
			//- Number of quadratures
            inline label nQuad() const;

            //- Const access to absorption coefficient
            inline const volScalarField& kappa() const;
			
	    //- Const access to scattering coefficient
            inline const volScalarField& sigma() const;

            //- Const access to wavelength total absorption coefficient
            inline const volScalarField& kappaLambda(const label lambdaI) const;

	    //- Const access to wavelength total scattering coefficient
            inline const volScalarField& sigmaLambda(const label lambdaI) const;
			
            //- Const access to incident radiation field
            inline const volScalarField& Glambda(const label lambdaI) const;

            //- Const access to total radiative heat flux field
            inline const volScalarField& QrLambda(const label lambdaI) const;

            //- Const access to incident radiative heat flux field
            inline const volScalarField& QinLambda(const label lambdaI) const;
			
			//- Const access to volumetric emitted intensity
			inline const volScalarField& Elambda(const label lambdaI) const;

            //- Const access to specific quadrature
            const quadrature& quadratureById(const label modelIndex) const;

            //- Const access to maxIter
            inline label maxIter() const;

            //- Const acces to convergence
            inline scalar convergence() const;
			
			tmp<volScalarField> Rp() const;
			
			tmp<volScalarField> Ru() const;
			
			scalarField getThermal(label patchI, label lambdaId) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "DORTI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
