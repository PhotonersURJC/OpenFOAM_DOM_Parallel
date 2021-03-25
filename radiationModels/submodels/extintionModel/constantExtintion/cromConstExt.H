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
    Foam::radiation::cromConstExtintion

Description
    Homogeneous radiation extintion coefficient for multi-band
	
	\verbatim
		extintionModel cromConstExt;
		cromConstExtCoeffs
		{
		   absorptivity
		   {
			   band1	absorptivity    [ m^-1 ]         0.2;
			   band2	absorptivity    [ m^-1 ]         0.5;
		   }
		   scatter
		   {
			   band1	scatter    [ m^-1 ]         0.3;
			   band2	scatter    [ m^-1 ]         0.4;
		   }
		   bandLimits
		   {
			   band1 (300 400);
			   band2 (400 500);
		   }
		}
	\endverbatim

SourceFiles
    cromConstExt.C

\*---------------------------------------------------------------------------*/

#ifndef radiationCromConstExt_H
#define radiationCromConstExt_H

#include "radExtintionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                  Class cromConstExt Declaration
\*---------------------------------------------------------------------------*/

class cromConstExt
:
    public radExtintionModel
{
    // Private data

        //- Absorption model dictionary
        dictionary coeffsDict_;

        //- Band absorption coefficient / [1/m]
        PtrList<dimensionedScalar> kappaLambda_;
		
		//- Band scattering coefficient / [1/m]
        PtrList<dimensionedScalar> sigmaLambda_;
		
		//- Number of wavelength bands
		label nBands_;
		
		//- Wavelength limits of bands
		List< Vector2D < scalar > > bandLimits_;

public:

    //- Runtime type information
    TypeName("cromConstExt");


    // Constructors

        //- Construct from components
        cromConstExt(const dictionary& dict, const fvMesh& mesh);


    //- Destructor
    virtual ~cromConstExt();


    // Member Functions

        // Access

            // Coefficients

                //- Absorption coefficient for continuous phase
                tmp<volScalarField> kappaCont(const label bandI) const;

                //- Scattering coefficient for continuous phase
                tmp<volScalarField> sigmaCont(const label bandI) const;
				
				virtual List< Vector2D < scalar > > getBandsLimits() const;
		

    // Member Functions

	        //- Correct absorption coefficients
            virtual void kappaCorrect
            (
                volScalarField& kappa,
                PtrList<volScalarField>& kappaj
            ) const;

			//- Correct scattering coefficients
            virtual void sigmaCorrect
            (
                volScalarField& sigma,
                PtrList<volScalarField>& sigmaj
            ) const;
	
        inline bool isGrey() const
        {
            return false;
        }

    inline bool isHomogeneous() const
    { return true;}
    inline bool isConstant() const
    { return true;}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
