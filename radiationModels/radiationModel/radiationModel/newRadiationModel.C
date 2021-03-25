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

#include "phaseFuncModel.H"
#include "radExtintionModel.H"
#include "newRadiationModel.H"
#include "basicThermo.H"
#include "fvmSup.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(newRadiationModel, 0);
        defineRunTimeSelectionTable(newRadiationModel, T);
        defineRunTimeSelectionTable(newRadiationModel, dictionary);
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //



Foam::IOobject Foam::radiation::newRadiationModel::createIOobject
(
    const fvMesh& mesh
) const
{
    IOobject io
    (
        "radiationProperties",
        mesh.time().constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (io.typeHeaderOk<dictionary>(true))
    {
        io.readOpt() = IOobject::MUST_READ_IF_MODIFIED;
        return io;
    }
    else
    {
        io.readOpt() = IOobject::NO_READ;
        return io;
    }
}


void Foam::radiation::newRadiationModel::initialise()
{
    if (radiation_)
    {
        solverFreq_ = max(1, lookupOrDefault<label>("solverFreq", 1));
	phaseFunc_.reset
        (
            phaseFuncModel::New(*this, mesh_).ptr()
        );
        extintion_.reset
        (
            radExtintionModel::New(*this, mesh_).ptr()
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::newRadiationModel::newRadiationModel(const volScalarField& T)
:
    IOdictionary
    (
        IOobject
        (
            "radiationProperties",
            T.time().constant(),
            T.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(T.mesh()),
    time_(T.time()),
    radiation_(false),
	T_(T),
    coeffs_(dictionary::null),
    solverFreq_(0),
    firstIter_(true),
    phaseFunc_(NULL),
    extintion_(NULL)
{}

Foam::radiation::newRadiationModel::newRadiationModel
(
    const word& type,
    const volScalarField& T
)
:
    IOdictionary(createIOobject(T.mesh())),
    mesh_(T.mesh()),
    time_(T.time()),
    radiation_(lookupOrDefault("radiation", true)),
	T_(T),
    coeffs_(subOrEmptyDict(type + "Coeffs")),
    solverFreq_(1),
    firstIter_(true),
    phaseFunc_(NULL),
    extintion_(NULL)
{
    if (readOpt() == IOobject::NO_READ)
    {
        radiation_ = false;
    }

    initialise();
}

Foam::radiation::newRadiationModel::newRadiationModel
(
    const word& type,
    const dictionary& dict,
    const volScalarField& T
)
:
    IOdictionary
    (
        IOobject
        (
            "radiationProperties",
            T.time().constant(),
            T.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dict
    ),
    mesh_(T.mesh()),
    time_(T.time()),
    radiation_(lookupOrDefault("radiation", true)),
	T_(T),
    coeffs_(subOrEmptyDict(type + "Coeffs")),
    solverFreq_(1),
    firstIter_(true),
    phaseFunc_(NULL),
    extintion_(NULL)
{
    initialise();
}


// * * * * * * * * * * * * * * * * Destructor    * * * * * * * * * * * * * * //

Foam::radiation::newRadiationModel::~newRadiationModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::radiation::newRadiationModel::read()
{
    if (regIOobject::read())
    {
        lookup("radiation") >> radiation_;
        coeffs_ = subOrEmptyDict(type() + "Coeffs");

        solverFreq_ = lookupOrDefault<label>("solverFreq", 1);
        solverFreq_ = max(1, solverFreq_);

        return true;
    }
    else
    {
        return false;
    }
}


bool Foam::radiation::newRadiationModel::correct()
{
    if (!radiation_)
    {
        return false;
    }

    if (firstIter_ || (time_.timeIndex() % solverFreq_ == 0))
    {
        firstIter_ = false;
        return calculate();
    }
    return false;
}

Foam::tmp<Foam::fvScalarMatrix> Foam::radiation::newRadiationModel::Sh
(
    const basicThermo& thermo,
    const volScalarField& he
) const
{
    const volScalarField Cpv(thermo.Cpv());
    const volScalarField T3(pow3(T_));

    return
    (
        Ru()
      - fvm::Sp(4.0*Rp()*T3/Cpv, he)
      - Rp()*T3*(T_ - 4.0*he/Cpv)
    );
}


Foam::tmp<Foam::fvScalarMatrix> Foam::radiation::newRadiationModel::ST
(
    const dimensionedScalar& rhoCp,
    volScalarField& T
) const
{
    return
    (
        Ru()/rhoCp
      - fvm::Sp(Rp()*pow3(T)/rhoCp, T)
    );
}


const Foam::radiation::phaseFuncModel&
Foam::radiation::newRadiationModel::phaseFunc() const
{
    if (!phaseFunc_.valid())
    {
        FatalErrorIn
        (
            "const Foam::noTradiation::phaseFuncModel&"
            "Foam::noTradiation::noTradiationModel::phaseFunc() const"
        )
            << "Requested radiation phase function model, but model is "
            << "not activate" << abort(FatalError);
    }

    return phaseFunc_();
}

const Foam::radiation::radExtintionModel&
Foam::radiation::newRadiationModel::extintion() const
{
    if (!extintion_.valid())
    {
        FatalErrorIn
        (
            "const Foam::noTradiation::radExtintionModel&"
            "Foam::noTradiation::noTradiationModel::extintion() const"
        )
            << "Requested radiation extintion model, but model is "
            << "not activate" << abort(FatalError);
    }

    return extintion_();
}

void Foam::radiation::newRadiationModel::sctrMatrixGen
(
	const List<vector> directions,
	const label model
)
{
	phaseFunc_().MatrixGen(directions,model);
}


// ************************************************************************* //
