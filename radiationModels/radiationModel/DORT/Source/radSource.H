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

Namespace
    Foam::radiation

Description
    Namespace for radiation modelling

Class
    Foam::radiation::radSource

Description
    Abstract class of emmision sources. It can be:
	 - isotropicRadSource
	 - constConeRadSource
	 - parallelRadSource
	 
	Sources are declared in region/constant, in radSource1, radSource2...
	
	i.e. dictionary
	\verbatim
        sourceType constConeRadSource;	//type of radiation source

		Direction	(-1.0 0.0 0.0);		//Central direction (parallel or cone)
		Opening		60;					//Cone opening in degrees (only cone)
		patchesList
		{
			lamp 100.0;					//boundary name and flux (W/m2)
		}
    \endverbatim

SourceFiles
    radSource.C
    radSourceNew.C

\*---------------------------------------------------------------------------*/

#ifndef radSource_H
#define radSource_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "volFieldsFwd.H"
#include "DimensionedField.H"
#include "fvMatricesFwd.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fluidThermo;
class fvMesh;

namespace radiation
{

/*---------------------------------------------------------------------------*\
                       Class radSource Declaration
\*---------------------------------------------------------------------------*/

class radSource
:
    public IOdictionary
{
protected:

    // Protected data

    //- Reference to the mesh database
    const fvMesh& mesh_;

    //- Reference to the time database
    const Time& time_;

	//- Source Index (incremented by DORT)
	label index_;
	
	//- Number of waveLength bands
	label nBands_;

	//- List of afected patches
	List<word> patchesID_;

	//- Subdictionaries containing bands fluxes
	List<dictionary> patchesDict_;
	
	//- List of radiation flux emitted per patch and band
	List<List<scalar> > patchesQ_;

	//- Direction of paralel radiation, central direction of cones, or z NULL if non-oriented sources
	vector direction_;

	//- Cone opening or 0 if non-oriented or paralel radiation
	scalar opening_;

private:

    // Private Member Functions

        //- Create IO object if dictionary is present
        IOobject createIOobject(const fvMesh& mesh, const label sourceIndex) const;

        //- Initialise
        void initialise();

        //- Disallow default bitwise copy construct
        radSource(const radSource&);

        //- Disallow default bitwise assignment
        void operator=(const radSource&);


public:

    //- Runtime type information
    TypeName("radSource");


    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            radSource,
            mesh,
            (
                const fvMesh& mesh,
				const label sourceIndex,
				const label nBands
            ),
            (mesh, sourceIndex, nBands)
        );
        declareRunTimeSelectionTable
        (
            autoPtr,
            radSource,
            dictionary,
            (
                const dictionary& dict,
                const fvMesh& mesh,
				const label sourceIndex,
				const label nBands
            ),
            (dict, mesh, sourceIndex, nBands)
        );
		

    // Constructors

        //- Null constructor
        radSource(const fvMesh& mesh, const label sourceIndex, const label nBands);

        //- Construct from components
        radSource(const word& type, const fvMesh& mesh, const label sourceIndex, const label nBands);

        //- Construct from components
        radSource
        (
            const word& type,
            const dictionary& dict,
            const fvMesh& mesh,
			const label sourceIndex,
			const label nBands
        );


    // Selectors

        //- Return a reference to the selected radiation model
        static autoPtr<radSource> New(const fvMesh& mesh, const label sourceIndex, const label nBands);

        //- Return a reference to the selected radiation model
        static autoPtr<radSource> New
        (
            const dictionary& dict,
            const fvMesh& mesh,
			const label sourceIndex,
			const label nBands
        );


    //- Destructor
    virtual ~radSource();


    // Member Functions

        // Edit


			//- Return main axis to draw the new quadrature
			virtual const List<vector> getAxis(label nTheta, label nPhi) = 0;

            //- Store values of emission
            virtual void generate(const List<vector> directions) = 0;

            //- Read radiationProperties dictionary
            virtual bool read() = 0;


        // Access

			//- Acces to specific stored value
            virtual scalar calculate(const word patchID, const label face, const label direction, const label band) const = 0;
			
            //- Access to source direction
            const vector& direction() const;

	    //- Access to source cone opening
	    const scalar& opening() const;
};


            //- Access to source direction
inline const Foam::vector& Foam::radiation::radSource::direction() const
{ return direction_; }

	    //- Access to source cone opening
inline const Foam::scalar& Foam::radiation::radSource::opening() const
{ return opening_; }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define addToSourceRunTimeSelectionTables(model)                        \
                                                                              \
    addToRunTimeSelectionTable                                                \
    (                                                                         \
        radSource,                                                            \
        model,                                                                \
        dictionary                                                            \
    );                                                                        \
                                                                              \
    addToRunTimeSelectionTable                                                \
    (                                                                         \
        radSource,                                                            \
        model,                                                                \
        mesh                                                              \
    );

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
