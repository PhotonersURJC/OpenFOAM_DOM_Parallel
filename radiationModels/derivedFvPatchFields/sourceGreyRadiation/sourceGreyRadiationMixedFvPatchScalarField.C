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

#include "sourceGreyRadiationMixedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "mappedPatchBase.H"

#include "DORT.H"
#include "quadrature.H"
#include "constants.H"

using namespace Foam::constant;
using namespace Foam::constant::mathematical;

void Foam::radiation::sourceGreyRadiationMixedFvPatchScalarField::initialise()
{
	word name = internalField().name();
	if(name.find_first_of("Lambda")==1)
	{
		const newRadiationModel& radiation =
        db().lookupObject<newRadiationModel>("radiationProperties");
		dort = refCast<const DORT>(radiation);
		
		label quadIndex = -1;
		dort.setRayModelIndex(internalField().name(), quadIndex);
		quad = const_cast<quadrature&>(dort.quadratureById(quadIndex));
		quad.setRayIdLambdaId(internalField().name(), rayId, lambdaId);
		patchI = patch_.index();
		direction = const_cast<vector>(quad.IRay(rayId).dAve());
		
		if(rayId==0)
		{
			//Information is sent to quadrature: booolean to inform if it has specular reflection, order of reflection, and patchId
			bool specUse = true;
			if(refCoeff_==0.0)
			{specUse = false;}
			if(refDifF_ == 1.0)
			{specUse = false;}
			quad.getBoundaryReport(specUse,secondOrder_,patchI,lambdaId);
		}
	}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::sourceGreyRadiationMixedFvPatchScalarField::
sourceGreyRadiationMixedFvPatchScalarField		
(									
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    patch_(p),
    refCoeff_(0.0),
    refDifF_(0.0),
    transCoeff_(0.0),
    transDifF_(0.0),
    sampleMesh_(""),
    samplePatchI_(""),
	secondOrder_(false)
{
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 1.0;
}


Foam::radiation::sourceGreyRadiationMixedFvPatchScalarField::	
sourceGreyRadiationMixedFvPatchScalarField		
(							
    const sourceGreyRadiationMixedFvPatchScalarField& ptf,		
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    patch_(p),
    refCoeff_(ptf.refCoeff_),
    refDifF_(ptf.refDifF_),
    transCoeff_(ptf.transCoeff_),
    transDifF_(ptf.transDifF_),
    sampleMesh_(ptf.sampleMesh_),
    samplePatchI_(ptf.samplePatchI_),
	secondOrder_(ptf.secondOrder_)
{
	initialise();
}


Foam::radiation::sourceGreyRadiationMixedFvPatchScalarField::	
sourceGreyRadiationMixedFvPatchScalarField				
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    patch_(p),
    refCoeff_(dict.lookupOrDefault<scalar>("refCoeff",0.0)),
    refDifF_(dict.lookupOrDefault<scalar>("refDifF",0.0)),
    transCoeff_(dict.lookupOrDefault<scalar>("transCoeff",0.0)),
    transDifF_(dict.lookupOrDefault<scalar>("transDifF",0.0)),
    sampleMesh_(dict.lookupOrDefault<word>("sampleRegion","")),
    samplePatchI_(dict.lookupOrDefault<word>("samplePatchI","")),
	secondOrder_(false)
{
	word orderText(dict.lookupOrDefault<word>("order","first"));
	if((orderText == "2") || (orderText == "second"))
	{ 
		secondOrder_ = true;
	}
	if (dict.found("refValue"))
    {
        fvPatchScalarField::operator=
        (
            scalarField("value", dict, p.size())
        );
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        refValue() = 0.0;
        refGrad() = 0.0;
        valueFraction() = 1.0;

        fvPatchScalarField::operator=(refValue());
    }
	initialise();
}


Foam::radiation::sourceGreyRadiationMixedFvPatchScalarField::	
sourceGreyRadiationMixedFvPatchScalarField		
(
    const sourceGreyRadiationMixedFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    patch_(ptf.patch_),
    refCoeff_(ptf.refCoeff_),
    refDifF_(ptf.refDifF_),
    transCoeff_(ptf.transCoeff_),
    transDifF_(ptf.transDifF_),
    sampleMesh_(ptf.sampleMesh_),
    samplePatchI_(ptf.samplePatchI_),
	secondOrder_(ptf.secondOrder_)
{
	initialise();
}


Foam::radiation::sourceGreyRadiationMixedFvPatchScalarField::	
sourceGreyRadiationMixedFvPatchScalarField			
(
    const sourceGreyRadiationMixedFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    patch_(ptf.patch_),
    refCoeff_(ptf.refCoeff_),
    refDifF_(ptf.refDifF_),
    transCoeff_(ptf.transCoeff_),
    transDifF_(ptf.transDifF_),
    sampleMesh_(ptf.sampleMesh_),
    samplePatchI_(ptf.samplePatchI_),
	secondOrder_(ptf.secondOrder_)
{
	initialise();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::radiation::sourceGreyRadiationMixedFvPatchScalarField::	
updateCoeffs()					
{									
	if (this->updated())
    {
        return;
    }

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.

    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;				
									
    scalarField& Iw = *this;

	//#OMP: Paralelizo??
    const scalarField nAve(patch_.nf() & direction);

	scalarField& Qr = const_cast<fvPatchField<double>&>(quad.futQrLambda(lambdaId).boundaryField()[patchI]);
	Qr += Iw*nAve;
	
	
	scalarField& Qin = const_cast<fvPatchField<double>&>(quad.futQinLambda(lambdaId).boundaryField()[patchI]);
	
    scalarField Ir = quad.QinLambda(lambdaId).boundaryField()[patchI];
    scalarField nbrI(Ir.size(),0.0);
    scalarField nbrFullIn(Ir.size(),0.0);
    if ((isA<mappedPatchBase>(patch_.patch())) & (transCoeff_>0.0))
    {
	const mappedPatchBase& mpp = refCast<const mappedPatchBase>(patch_.patch());

        const polyMesh& nbrMesh = mpp.sampleMesh();

        const radiation::newRadiationModel& nbrRadiation = nbrMesh.lookupObject<radiation::newRadiationModel>
        (
            "radiationProperties"
        );

        const DORT& nbrDort(refCast<const DORT>(nbrRadiation));

        const quadrature& nbrQuad = nbrDort.quadratureById(quadIndex);

        discreteOrdinate& nbrRay =
        const_cast<discreteOrdinate&>(nbrQuad.IRay(rayId));

	scalarField nbrIn(nbrRay.ILambda(lambdaId).boundaryField()[mpp.samplePolyPatch().index()]);

        mpp.distribute(nbrIn);
	
	nbrI += nbrIn;
        
    }
	scalarField specularReflected = quad.getSpec(patchI, rayId, lambdaId)*(refCoeff_*(1.0-refDifF_));
	scalarField thermalEmission = specularReflected;
	if(dort.thermo())
		thermalEmission = (1.0-transCoeff_-refCoeff_)*dort.getThermal(patchI,lambdaId);
    forAll(Iw, faceI)
    {
        if (nAve < 0.0)
        {
            refGrad()[faceI] = 0.0;
            valueFraction()[faceI] = 1.0;
            refValue()[faceI] = quad.sourceFlux(patch_.name(),faceI,rayId,lambdaId);	
			if(transCoeff_>0.0)
            {
                refValue()[faceI] += nbrI[faceI]*transCoeff_*(1.0-transDifF_);	
				refValue()[faceI] += nbrFullIn[faceI]*transCoeff_*transDifF_/pi;
			}
            if(refCoeff_>0.0)
            {
				refValue()[faceI] += Ir[faceI]*refCoeff_*refDifF_/pi;
				refValue()[faceI] += specularReflected[faceI];
			}
			if(dort.thermo())
				refValue()[faceI] += thermalEmission[faceI];
        }
        else
        {
            // direction into the wall
            valueFraction()[faceI] = 0.0;
            refGrad()[faceI] = 0.0;
            refValue()[faceI] = 0.0; //not used

            Qin[faceI] += Iw[faceI]*nAve[faceI];			
        }
    }
    // Restore tag
    UPstream::msgType() = oldTag;

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::radiation::sourceGreyRadiationMixedFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{
    makePatchTypeField
    (
        fvPatchScalarField,
        sourceGreyRadiationMixedFvPatchScalarField
    );
}
}


// ************************************************************************* //
