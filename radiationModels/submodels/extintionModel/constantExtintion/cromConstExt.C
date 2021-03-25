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

#include "cromConstExt.H"
#include "addToRunTimeSelectionTable.H"
#include "fvm.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(cromConstExt, 0);

        addToRunTimeSelectionTable
        (
            radExtintionModel,
            cromConstExt,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::cromConstExt::cromConstExt
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    radExtintionModel(dict, mesh),
	coeffsDict_(dict.subDict(typeName + "Coeffs")),
	kappaLambda_(1),
	sigmaLambda_(1),
	nBands_(1),
	bandLimits_(0)
{
	
	if(coeffsDict_.isDict("absorptivity")) //Band absorption coefficients
	{
		dictionary kappaDict(coeffsDict_.lookup("absorptivity"));
		wordList bandNames = kappaDict.toc();
		nBands_=bandNames.size();
		kappaLambda_.setSize(nBands_);
		for(label lambdaI = 0; lambdaI<nBands_;lambdaI++)
		{
			kappaLambda_.set(lambdaI,new dimensionedScalar(kappaDict.lookup(bandNames[lambdaI])));
		}
	}
	else //Single absorption coefficient for all bands
	{
		kappaLambda_.set(0,new dimensionedScalar(coeffsDict_.lookup("absorptivity")));
	}
	if(coeffsDict_.isDict("scatter"))
	{
		dictionary sigmaDict(coeffsDict_.lookup("scatter"));
		wordList bandNames = sigmaDict.toc();
		nBands_=bandNames.size();
		sigmaLambda_.setSize(nBands_);
		for(label lambdaI = 0; lambdaI<nBands_;lambdaI++)
		{
			sigmaLambda_.set(lambdaI,new dimensionedScalar(sigmaDict.lookup(bandNames[lambdaI])));
		}
	}
	else
	{
		sigmaLambda_.set(0,new dimensionedScalar(coeffsDict_.lookup("scatter")));
	}
	if(coeffsDict_.isDict("bandLimits"))
	{
		dictionary limitsDict(coeffsDict_.lookup("bandsLimits"));
		wordList bandNames = limitsDict.toc();
		nBands_=bandNames.size();
		bandLimits_.setSize(nBands_);
		for(label lambdaI = 0; lambdaI<nBands_;lambdaI++)
		{
			bandLimits_[lambdaI]=vector2D(limitsDict.lookup(bandNames[lambdaI]));
		}
		if(nBands_ != kappaLambda_.size())
		{
			Info << "Number of bands in extintion do not match" << endl;
		}
	}
	if(kappaLambda_.size() != sigmaLambda_.size())
	{
		Info << "Number of bands in extintion do not match" << endl;
	}
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::cromConstExt::~cromConstExt()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiation::cromConstExt::kappaCont(const label bandI) const
{
    tmp<volScalarField> ta
    (
        new volScalarField
        (
            IOobject
            (
                "kappa_" + Foam::name(bandI),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            kappaLambda_[bandI]
        )
    );

    return ta;
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::cromConstExt::sigmaCont(const label bandI) const
{
    tmp<volScalarField> ts
    (
        new volScalarField
        (
            IOobject
            (
                "sigma_" + Foam::name(bandI),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            sigmaLambda_[bandI]
        )
    );

    return ts;
}

void Foam::radiation::cromConstExt::kappaCorrect
(
    volScalarField& kappa,
    PtrList<volScalarField>& kappaj
) const
{
	kappaj[0] = kappaCont(0);
	kappa = kappaj[0];
	for(label lambdaI = 1; lambdaI<nBands_;lambdaI++)
	{
		kappaj[lambdaI]=kappaCont(lambdaI);
		kappa+=kappaj[lambdaI];
	}
}

void Foam::radiation::cromConstExt::sigmaCorrect
(
    volScalarField& sigma,
    PtrList<volScalarField>& sigmaj
) const
{
	sigmaj[0] = sigmaCont(0);
	sigma = sigmaj[0];
	for(label lambdaI = 1; lambdaI<nBands_;lambdaI++)
	{
		sigmaj[lambdaI]=sigmaCont(lambdaI);
		sigma+=sigmaj[lambdaI];
	}
}

Foam::List< Foam::Vector2D < Foam::scalar > > Foam::radiation::cromConstExt::getBandsLimits() const
{
	return bandLimits_;
}


// ************************************************************************* //
