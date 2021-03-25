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
    Foam::radiation::absorptionEmissionModel

Description
    Abstract class to supply absorption and scattering coefficients 
	for radiation modelling.
	Current models:
	 - greyConstExt: homogeneous coefficients in gray simulations
	 - cromConstExt: homogeneous coefficients in multi-band simulations
	 
	Defined in radiationProperties
	i.e. dictionaries in specific models

SourceFiles
    radExtintionModel.C
	radExtintionModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef noTradiationRadExtintionModel_H
#define noTradiationRadExtintionModel_H

#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// Forward declaration of classes
class fvMesh;

namespace radiation
{

/*---------------------------------------------------------------------------*\
                   Class extintionModel Declaration
\*---------------------------------------------------------------------------*/

class radExtintionModel
{

protected:

    // Protected data

        //- Radiation model dictionary
        const dictionary dict_;

        //- Reference to the fvMesh
        const fvMesh& mesh_;

	


public:

    //- Runtime type information
    TypeName("radExtintionModel");

    //- Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            radExtintionModel,
            dictionary,
            (
                const dictionary& dict,
                const fvMesh& mesh
            ),
            (dict, mesh)
        );

    // Constructors

        //- Construct from components
        radExtintionModel
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Selector
    static autoPtr<radExtintionModel> New
    (
        const dictionary& dict,
        const fvMesh& mesh
    );

    //- Destructor
    virtual ~radExtintionModel();


    // Member Functions

        // Access

            //- Reference to the mesh
            inline const fvMesh& mesh() const
            {
                return mesh_;
            }

            //- Reference to the dictionary
            inline const dictionary& dict() const
            {
                return dict_;
            }


            // Absorption coefficient

                //- Absorption coefficient (net)
                virtual tmp<volScalarField> kappa(const label bandI = 0) const;

                //- Absorption coefficient for continuous phase
                virtual tmp<volScalarField> kappaCont(const label bandI = 0) const;

                //- Absorption coefficient for dispersed phase
                virtual tmp<volScalarField> kappaDisp(const label bandI = 0) const;
				
				volScalarField fullKappa() const;

	    // Scattering coefficient

                //- Scattering coefficient (net)
                virtual tmp<volScalarField> sigma(const label bandI = 0) const;

                //- Scattering coefficient for continuous phase
                virtual tmp<volScalarField> sigmaCont(const label bandI = 0) const;

                //- Scattering coefficient for dispersed phase
                virtual tmp<volScalarField> sigmaDisp(const label bandI = 0) const;
				
		//Bands limits
		
				//- Const access to bands limits list
				inline virtual List< Vector2D < scalar > > getBandsLimits() const 
				{return List< Vector2D < scalar > >();}
				
            //- Const access to the number of bands - defaults to 1 for grey
            virtual label nBands() const;

            //- Flag for whether the extintion is for a grey gas
            virtual bool isGrey() const;

            virtual bool isHomogeneous() const;

            virtual bool isConstant() const;

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
			
};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace noTradiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

