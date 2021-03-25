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

#include "DORT.H"
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
        defineTypeNameAndDebug(DORT, 0);
        addToRadiationRunTimeSelectionTables(DORT);
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::radiation::DORT::initialise(const dictionary& dict)
{
    // Construct absorption field for each wavelength
    forAll(kappaLambda_, lambdaI)
    {
        kappaLambda_.set
        (
            lambdaI,
            new volScalarField
            (
                IOobject
                (
                    "kappaLambda_" + Foam::name(lambdaI) ,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                kappa_
            )
        );
    }
    // Construct scattering field for each wavelength
    forAll(sigmaLambda_, lambdaI)
    {
        sigmaLambda_.set
        (
            lambdaI,
            new volScalarField
            (
                IOobject
                (
                    "sigmaLambda_" + Foam::name(lambdaI) ,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                sigma_
            )
        );
    }
	forAll(Glambda_,lambdaI)
	{
		Glambda_.set
        (
            lambdaI,
            new volScalarField
            (
				IOobject
				(
					"G_global_" + Foam::name(lambdaI),
					mesh_.time().timeName(),
					mesh_,
					IOobject::NO_READ,
					IOobject::AUTO_WRITE
				),
				mesh_,
				dimensionedScalar("G_global_" + Foam::name(lambdaI), dimMass/pow3(dimTime), 0.0)
			)
        );
		QrLambda_.set
        (
            lambdaI,
            new volScalarField
            (
				IOobject
				(
					"Qr_global_" + Foam::name(lambdaI),
					mesh_.time().timeName(),
					mesh_,
					IOobject::NO_READ,
					IOobject::AUTO_WRITE
				),
				mesh_,
				dimensionedScalar("Qr_global_" + Foam::name(lambdaI), dimMass/pow3(dimTime), 0.0)
			)
        );
		QinLambda_.set
        (
            lambdaI,
            new volScalarField
            (
				IOobject
				(
					"Qin_global_" + Foam::name(lambdaI),
					mesh_.time().timeName(),
					mesh_,
					IOobject::NO_READ,
					IOobject::AUTO_WRITE
				),
				mesh_,
				dimensionedScalar("Qin_global_" + Foam::name(lambdaI), dimMass/pow3(dimTime), 0.0)
			)
        );
		eAlambda_.set
        (
            lambdaI,
            new volScalarField
            (
				IOobject
				(
					"eAlambda_" + Foam::name(lambdaI),
					mesh_.time().timeName(),
					mesh_,
					IOobject::NO_READ,
					IOobject::AUTO_WRITE
				),
				mesh_,
				dimensionedScalar("eAlambda_" + Foam::name(lambdaI), dimMass/(pow3(dimTime)*dimLength), 0.0)
			)
        );
	}
	if(nQuad_>0) //DORTsolver, only radiation
	{
		extintion_.reset
        (
            radExtintionModel::New(dict, mesh_).ptr()
        );
		phaseFunc_.reset
        (
            phaseFuncModel::New(dict, mesh_).ptr()
        );
	}
	if(max(T_).value()>0.0) //temperature taken into account
	{
		bandsLimits_ = extintion_->getBandsLimits();
		thermo_ = true;
		forAll(Elambda_, lambdaI)
		{
			Elambda_.set
			(
				lambdaI,
				new volScalarField
				(
					IOobject
					(
						"E_lambda_" + Foam::name(lambdaI),
						mesh_.time().timeName(),
						mesh_,
						IOobject::NO_READ,
						IOobject::NO_WRITE
					),
					mesh_,
					dimensionedScalar("E_lambda_" + Foam::name(lambdaI), dimMass/(pow3(dimTime)*dimLength), 0.0)
				)
			);
		}
	}
	extintion_->kappaCorrect(kappa_, kappaLambda_);
	extintion_->sigmaCorrect(sigma_, sigmaLambda_);
	bool scatter = false;
	for(label lambdaI = 0; lambdaI < nLambda_; lambdaI++)
		if(max(sigmaLambda_[lambdaI]).value()>0.0)
			scatter=true;
	bool reading = true; label indexer = 0; bool serialQuads = (nQuad_==0);
    while (reading)
    {
        IOobject sourceHeader
        (
            "radSource" + Foam::name(indexer+1),
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ,
			IOobject::NO_WRITE
        );
	    if (sourceHeader.typeHeaderOk<dictionary>(true))
	    {
			indexer++;
	        submodels_.resize(submodels_.size()+1);
	        sources_.resize(sources_.size()+1);
	        sources_.set(sources_.size()-1,Foam::radiation::radSource::New(mesh_,indexer,nLambda_));
			if (nQuad_ == 0) //General Solver
	        {
				submodels_.set(
	            submodels_.size()-1,
	            new quadrature(*this, mesh_, indexer, sources_[sources_.size()-1],serialQuads));
				submodels_[submodels_.size()-1].initialise();
				if(scatter)
				{
					sctrMatrixGen
					(
						submodels_[sources_.size()-1].directions(),
						sources_.size()
					);
				}
		    }
	    }
	    else {reading = false;}
    }
	if(nQuad_ > 0) //Specific Solver
	{
		if (nQuad_ < indexer)
		{
		    FatalErrorIn
            (
                "DORT::initialise()"
            )   << "There are more generated quadratures ("
				<< indexer << ") than specified ones ("
                << nQuad_ << ")" << exit(FatalError);
		}
		if (nQuad_ > indexer)
		{
		    FatalErrorIn
            (
                "DORT::initialise()"
            )   << "There are less generated quadratures ("
				<< indexer << ") than specified ones ("
                << nQuad_ << ")" << exit(FatalError);
		}
	}
}


void Foam::radiation::DORT::readOwnDict(const dictionary& dict)
{
	nTheta_ = readLabel(dict.lookup("nTheta"))*2;
    nPhi_ = readLabel(dict.lookup("nPhi"));
    nLambda_ = dict.lookupOrDefault<label>("nBands",1);	//If no reference to nBands, it's gray
	nQuad_ = readLabel(dict.lookup("nQuad"));
	modelsG_.resize(nQuad_);
	modelsQr_.resize(nQuad_);
	modelsQin_.resize(nQuad_);
	for(label i = 0; i<nQuad_; i++)
	{
		modelsG_.set(
		    i,
			new volScalarField(
			    IOobject
                (
                    "G_Quad_" + Foam::name(i),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
				mesh_)
		);
		modelsQr_.set(
		    i,
			new volScalarField(
			    IOobject
                (
                    "Qr_Quad_" + Foam::name(i),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
				mesh_)
		);
		modelsQin_.set(
		    i,
			new volScalarField(
			    IOobject
                (
                    "Qin_Quad_" + Foam::name(i),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
				mesh_)
		);
	}
    kappaLambda_.resize(nLambda_);
    sigmaLambda_.resize(nLambda_);
	Glambda_.resize(nLambda_);
	QrLambda_.resize(nLambda_);
	QinLambda_.resize(nLambda_);
	eAlambda_.resize(nLambda_);
	Elambda_.resize(nLambda_);
	nLambda_= 1;
    convergence_ = dict.lookupOrDefault<scalar>("convergence", 0.0);
    maxIter_=(dict.lookupOrDefault<label>("inScatterIters", 50));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::DORT::DORT(const volScalarField& T)
:
    newRadiationModel(typeName, T),
	nLambda_(extintion_->nBands()),
	bandsLimits_(0),
    kappa_
    (
        IOobject
        (
            "kappa",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("kappa", dimless/dimLength, 0.0)
    ),
    sigma_
    (
        IOobject
        (
            "sigma",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("sigma", dimless/dimLength, 0.0)
    ),
	eA_
	(
	    IOobject
        (
            "eA",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("eA", dimMass/(pow3(dimTime)*dimLength), 0.0)
    ),
	Glambda_(nLambda_),
	QrLambda_(nLambda_),
	QinLambda_(nLambda_),
	eAlambda_(nLambda_),
	modelsG_(0),
	modelsQr_(0),
	modelsQin_(0),
    nTheta_(readLabel(coeffs_.lookup("nTheta"))*2),
    nPhi_(readLabel(coeffs_.lookup("nPhi"))),
	nQuad_(0),
	actualBand_(0),
    kappaLambda_(nLambda_),
    sigmaLambda_(nLambda_),
	Elambda_(nLambda_),
    submodels_(0),
	sources_(0),
    convergence_(coeffs_.lookupOrDefault<scalar>("convergence", 0.0)),
    maxIter_(coeffs_.lookupOrDefault<label>("maxIter", 50)),
	thermo_(false)
{
    initialise(*this);
}

Foam::radiation::DORT::DORT
(
    const dictionary& dict,
    const volScalarField& T
)
: //Null constructor to take information from DORTproperties instead
    newRadiationModel(T),
    nLambda_(0),
	bandsLimits_(0),
    kappa_
    (
        IOobject
        (
            "kappa",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("kappa", dimless/dimLength, 0.0)
    ),
    sigma_
    (
        IOobject
        (
            "sigma",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("sigma", dimless/dimLength, 0.0)
    ),
	eA_
	(
	    IOobject
        (
            "eA",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("eA", dimMass/(pow3(dimTime)*dimLength), 0.0)
    ),
	Glambda_(0),
	QrLambda_(0),
	QinLambda_(0),
	eAlambda_(0),
    modelsG_(0),
	modelsQr_(0),
	modelsQin_(0),
    nTheta_(0),
    nPhi_(0),
	nQuad_(0),
	actualBand_(0),
    kappaLambda_(0),
    sigmaLambda_(0),
	Elambda_(0),
	submodels_(0),
	sources_(0),
    convergence_(0.0),
    maxIter_(0),
	thermo_(false)
{
	readOwnDict(dict);
    initialise(dict);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::DORT::~DORT()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::radiation::DORT::read()
{
    if (newRadiationModel::read())
    {
        // Only reading solution parameters - not changing ray geometry
        coeffs_.readIfPresent("convergence", convergence_);
        coeffs_.readIfPresent("maxIter", maxIter_);
        return true;
    }
    else
        return false;
}

void Foam::radiation::DORT::createQuad(const label model)
{
	bool serialQuads = (nQuad_==0);
	submodels_.set(model,new quadrature(*this, mesh_, model, sources_[model], serialQuads));
	actualBand_=0;
	bool scatter = false;
	forAll(sigmaLambda_,lambdaI)
	{
		if(max(sigmaLambda_[lambdaI]).value()!=0.0)
			scatter=true;
	}
	if(scatter)
	{
		sctrMatrixGen
		(
			submodels_[sources_.size()-1].directions(),
			sources_.size()
		);
	}
}

void Foam::radiation::DORT::deleteQuad()
{
	submodels_.clear();
	submodels_.resize(nQuad_);
}

void Foam::radiation::DORT::changeBand(const label model, const label band)
{
	Glambda_[actualBand_] += modelsG_[model];
	QrLambda_[actualBand_] += modelsQr_[model];
	QinLambda_[actualBand_] += modelsQin_[model];
	eAlambda_[actualBand_] = Glambda_[actualBand_]*kappaLambda_[actualBand_];
	eA_ += modelsG_[model]*kappaLambda_[actualBand_];
	actualBand_=band;
	submodels_[model].changeBand(band);
}

bool Foam::radiation::DORT::calculateQuad(const label model)
{
	if(thermo_)
	{
		if(bandsLimits_.size()==0)
				setNoBandEmission();
			else
				setBandEmission(actualBand_);
	}
	bool Converged = submodels_[model].calculate();	//Quadrature has already set resolution to single band when created
	modelsG_.set(model,submodels_[model].Glambda(actualBand_));
	modelsQin_.set(model,submodels_[model].QinLambda(actualBand_));
	modelsQr_.set(model,submodels_[model].QrLambda(actualBand_));
	return Converged;
}

bool Foam::radiation::DORT::calculate()
{
	if(thermo_)
	{
		forAll(Elambda_,lambdaI)
		{
			if(bandsLimits_.size()==0)
				setNoBandEmission();
			else
				setBandEmission(lambdaI);
		}
	}
	forAll(Glambda_,lambdaI)
	{
		Glambda_[lambdaI] = dimensionedScalar("zero",dimMass/pow3(dimTime), 0.0);
		QrLambda_[lambdaI] = dimensionedScalar("zero",dimMass/pow3(dimTime), 0.0);
		QinLambda_[lambdaI] = dimensionedScalar("zero", dimMass/pow3(dimTime), 0.0);
	}
	eA_ = dimensionedScalar("zero", dimMass/(pow3(dimTime)*dimLength), 0.0);
    extintion_->kappaCorrect(kappa_, kappaLambda_);
    extintion_->sigmaCorrect(sigma_, sigmaLambda_);
    bool Converged = true;
    forAll(submodels_, quadI)
    {
	    if (!submodels_[quadI].calculate())
            Converged = false;
		forAll(Glambda_,lambdaI)
		{
			Glambda_[lambdaI] += submodels_[quadI].Glambda(lambdaI);
			QrLambda_[lambdaI] += submodels_[quadI].QrLambda(lambdaI);
			QinLambda_[lambdaI] += submodels_[quadI].QinLambda(lambdaI);
		}
    }
	forAll(eAlambda_,lambdaI)
	{
		eAlambda_[lambdaI]=Glambda_[lambdaI]*kappaLambda_[lambdaI];
		eA_+=eAlambda_[lambdaI];
	}
    return Converged;
}

void Foam::radiation::DORT::setRayModelIndex
(
    const word& name,
    label& modelIndex
) const
{
    // assuming name is in the form: CHARS_rayId_lambdaId
    size_t i1 = name.find_first_of("_");
    size_t i2 = name.find_first_of("_",i1+1);
    modelIndex = Foam::readLabel(IStringStream(name.substr(i1+1, i2-1))());
}

const Foam::radiation::quadrature&
Foam::radiation::DORT::quadratureById(const label modelIndex) const
{
    return submodels_[modelIndex-1];
}

bool Foam::radiation::DORT::constantExt() const
{
    return extintion_->isConstant() && extintion_->isHomogeneous();
}

bool Foam::radiation::DORT::thermo() const
{
	return thermo_;
}

void Foam::radiation::DORT::generateDicts
( const label band, IOdictionary& newDict ) const
{
	newDict.add(keyType("nTheta"),nTheta_);
	newDict.add(keyType("nPhi"),nPhi_);
	newDict.add(keyType("nBands"),1);
	newDict.add(keyType("nQuad"),1);
	newDict.add(keyType("inScatterIters"),maxIter_);
	dictionary Parallelization = new dictionary();
	Parallelization.add(keyType("spatial"),"false");
	Parallelization.add(keyType("quadrature"),"false");
	Parallelization.add(keyType("band"),"false");
	newDict.add(keyType("parallelization"),Parallelization);
	newDict.add(keyType("extintionModel"),lookup("extintionModel"));	
}

void Foam::radiation::DORT::setBandEmission
(Foam::label lambdaI)
{
	scalar lambda0 = bandsLimits_[lambdaI][0];
	scalar lambda1 = bandsLimits_[lambdaI][1];
	scalar C10 = 7.49e-16*(lambda1-lambda0)/pi;
	scalar C20 = 1.4385e-2/lambda0;
	scalar C21 = 1.4385e-2/lambda1;
	scalar C11 = C10*pow(lambda1,-5.0);
	C10 *= pow(lambda0,-5.0);
	forAll(T_, i)
		Elambda_[i]=(C10/(exp(C20/T_[i])-1.0) + C11/(exp(C21/T_[i])-1.0))*kappaLambda_[lambdaI][i];
}

void Foam::radiation::DORT::setNoBandEmission ()
{
	Elambda_[0] = dimensionedScalar("zero",dimMass/(pow3(dimTime)*dimLength), 0.0);
	Elambda_[0] += (T_*T_*T_*T_)*Foam::constant::physicoChemical::sigma*kappa_/pi;
}

Foam::scalarField Foam::radiation::DORT::getThermal(Foam::label patchI, Foam::label lambdaId) const
{
	if(bandsLimits_.size()==0){
		scalarField T = T_.boundaryField()[patchI];
		return (T * T * T * T)*Foam::constant::physicoChemical::sigma.value()/pi;
	}
	scalarField T(T_.boundaryField()[patchI].size(),0.0);
	T += T_.boundaryField()[patchI];
	scalar lambda0 = bandsLimits_[lambdaId][0];
	scalar lambda1 = bandsLimits_[lambdaId][1];
	scalar C10 = 7.49e-16*(lambda1-lambda0)/pi;
	scalar C20 = 1.4385e-2/lambda0;
	scalar C21 = 1.4385e-2/lambda1;
	scalar C11 = C10*pow(lambda1,-5.0);
	C10 *= pow(lambda0,-5.0);
	forAll(T, i)
		T[i]=C10/(exp(C20/T[i])-1.0) + C11/(exp(C21/T[i])-1.0);
	return T;
}

Foam::tmp<Foam::volScalarField> Foam::radiation::DORT::Rp() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "Rp",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            4*physicoChemical::sigma*extintion_->fullKappa()
        )
    );
}

Foam::tmp<Foam::volScalarField> Foam::radiation::DORT::Ru() const
{
    volScalarField G = Glambda_[0];

    volScalarField E = Elambda_[0];
	for(int i = 1; i < nLambda_; i++)
	{
		G+=Glambda_[i];
		E+=Elambda_[i];
	}
    const volScalarField a = extintion_->fullKappa();

    return a*G - E;
}

// ************************************************************************* //
