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

#include "quadrature.H"
#include "fullDO.H"
#include "heteroExtDO.H"
#include "homoExtDO.H"
#include "thermoAbsorptionDO.H"
#include "absorptionDO.H"
#include "cleanDO.H"
#include "phaseFuncModel.H"
#include "radExtintionModel.H"
#include "constants.H"
#include "fvm.H"

using namespace Foam::constant;
using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::radiation::quadrature::createPtrLists
(PtrList<volScalarField>& currentField, word fieldName)
{
	for(label lambdaI = 0; lambdaI < nLambda_; lambdaI++)
	{
		currentField.set
		(
			lambdaI,
			new volScalarField
			(
				IOobject
				(
					fieldName+Foam::name(modelNumber_)+"_lambda_" + Foam::name(lambdaI),
					mesh_.time().timeName(),
					mesh_,
					IOobject::NO_READ,
					IOobject::AUTO_WRITE
				),
				mesh_,
				dimensionedScalar
				(
					fieldName+Foam::name(modelNumber_)+"_lambda_" + Foam::name(lambdaI),
					dimMass/pow3(dimTime),
					0.0
				)
			)
		);
	}
}

void Foam::radiation::quadrature::cleanFutures()
{
	forAll(futQrLambda_,lambdaI)
	{
		futQrLambda_[lambdaI] = dimensionedScalar("zero",dimMass/pow3(dimTime), 0.0);
		futQinLambda_[lambdaI] = dimensionedScalar("zero",dimMass/pow3(dimTime), 0.0);
	}
}
void Foam::radiation::quadrature::cleanPreCalculateds()
{
	forAll(preCalcQrLambda_,lambdaI)
	{
		preCalcQrLambda_[lambdaI] = dimensionedScalar("zero",dimMass/pow3(dimTime), 0.0);
		preCalcQinLambda_[lambdaI] = dimensionedScalar("zero",dimMass/pow3(dimTime), 0.0);
	}
}
void Foam::radiation::quadrature::cleanConvergeds()
{
	forAll(convQrLambda_,lambdaI)
	{
		convQrLambda_[lambdaI] = dimensionedScalar("zero",dimMass/pow3(dimTime), 0.0);
		convQinLambda_[lambdaI] = dimensionedScalar("zero",dimMass/pow3(dimTime), 0.0);
	}
}

void Foam::radiation::quadrature::initialise()
{
	createPtrLists(Glambda_,"G_");
	createPtrLists(QrLambda_,"Qr_");
	createPtrLists(futQrLambda_,"futQr_");
	createPtrLists(preCalcQrLambda_,"preCalcQr_");
	createPtrLists(convQrLambda_,"convQr_");
	createPtrLists(QinLambda_,"Qin_");
	createPtrLists(futQinLambda_,"futQin_");
	createPtrLists(preCalcQinLambda_,"preCalcQin_");
	createPtrLists(convQinLambda_,"convQin_");
	List<bool> conditions(4,false); //heterogeneous, thermal, scattering, absorption
	conditions[0]=!DORT_.constantExt();
	conditions[1]=DORT_.thermo();
	for(label lambdaI = 0; lambdaI < nLambda_; lambdaI++)
	{
		if(max(DORT_.kappaLambda(lambdaI)).value()!=0.0)
			conditions[3]=true;;
		if(max(DORT_.sigmaLambda(lambdaI)).value()!=0.0)
		{
			conditions[2]=true;
			scatter_=true;
		}
	}
    
    scalar deltaPhi = pi/(2.0*nPhi_);
    scalar deltaTheta = pi/nTheta_;
    label i = 0;
    const List<vector> mainAxis = radSource_.getAxis(nTheta_, nPhi_);
	label discOrdBands = nLambda_;
	if(serialQuad_)		//if serial Mode, discrete Ordinates have only one wavelenght band
		discOrdBands = 1;
	scalar opening = radSource_.opening();
    label cut = opening/deltaTheta;	//Number of divisions before the cut
	scalar planeOpening = 0.0;
	if(mesh_.nSolutionD()==2)
	{
		//establezco ejes
		planeOpening = opening;
		opening = 0.0;
		deltaTheta = pi;
		nTheta_ = 1;
	}
	if(mesh_.nSolutionD()==1)
	{
		opening = 0.0;
		deltaTheta = pi;
		nTheta_ = 1;
		nPhi_ = 2;
		deltaPhi = pi;
	}
	nRay_ = 4*nPhi_*nTheta_;
    IRay_.setSize(nRay_);
	if(cut<1){cut=1;}
    for (label n = 1; n <= nTheta_; n++)
    {
        if (opening > 0.0)
        {
			if (n <= cut) { deltaTheta = opening/cut; }
			else { deltaTheta = (pi-opening)/(nTheta_-cut); }
        }
        scalar thetai = (2.0*n - 1.0)*deltaTheta/2.0;
		if (n > cut && opening > 0.0) 
		{ 
			thetai = opening + (2.0*(n-cut) - 1.0)*deltaTheta/2.0;
		}
        for (label m = 1; m <= 4*nPhi_; m++)
        { 
            scalar phii = (2.0*m - 1.0)*deltaPhi/2.0;
			if(conditions[1] && (conditions[0] || conditions[2])) //heterogeneous & thermo or scattering & thermo
			{
				IRay_.set
                (
                    i,
                    new fullDO
                    (
						DORT_,
                        *this,
                        mesh_,
                        phii,
                        thetai,
                        deltaPhi,
                        deltaTheta,
                        discOrdBands,
                        i,
						mainAxis
                    )
                );
			}
			else if (conditions[0]) //heterogeneous, not thermo
			{
				IRay_.set
                (
                    i,
                    new heteroExtDO
                    (
						DORT_,
                        *this,
                        mesh_,
                        phii,
                        thetai,
                        deltaPhi,
                        deltaTheta,
                        discOrdBands,
                        i,
						mainAxis
                    )
                );
			}
			else if(conditions[2]) //homogeneous, with scattering, not thermo
			{
				IRay_.set
                (
                    i,
                    new homoExtDO
                    (
						DORT_,
                        *this,
                        mesh_,
                        phii,
                        thetai,
                        deltaPhi,
                        deltaTheta,
                        discOrdBands,
                        i,
						mainAxis,
                        max(DORT_.kappaLambda(0)).value() + max(DORT_.sigmaLambda(0)).value()
                    )
                );
			}
			else if(conditions[1]) //thermo, with or without homogeneous absorption
			{
				IRay_.set
                (
                    i,
                    new thermoAbsorptionDO
                    (
						DORT_,
                        *this,
                        mesh_,
                        phii,
                        thetai,
                        deltaPhi,
                        deltaTheta,
                        discOrdBands,
                        i,
						mainAxis,
                        max(DORT_.kappaLambda(0)).value()
                    )
                );
			}
			else if(conditions[3]) //only homogeneous absorption
			{
				IRay_.set
                (
                    i,
                    new absorptionDO
                    (
						DORT_,
                        *this,
                        mesh_,
                        phii,
                        thetai,
                        deltaPhi,
                        deltaTheta,
                        discOrdBands,
                        i,
						mainAxis,
                        max(DORT_.kappaLambda(0)).value()
                    )
                );
			}
			else //all conditions false, transparent media
			{
				IRay_.set
               (
                    i,
                    new cleanDO
                    (
						DORT_,
                        *this,
                        mesh_,
                        phii,
                        thetai,
                        deltaPhi,
                        deltaTheta,
                        discOrdBands,
                        i,
						mainAxis
                    )
                );
			}
            i++;
        }
    }
    // Construct scattering Matrix
    Info << "Solving optics..." << endl;
    directions_.setSize(nRay_);
    forAll(IRay_,rayI)
    {
		directions_[rayI]=IRay_[rayI].d();
    }

	if(lastSpecked_>0)
	{
		specularGuiding(mainAxis, cut, opening);
	}
	
    Info<< "quadrature : Allocated " << IRay_.size()
        << " rays with average orientation:" << nl;

    forAll(IRay_, rayId)
    {
        if (omegaMax_ <  IRay_[rayId].omega())
        {
            omegaMax_ = IRay_[rayId].omega();
        }
        Info<< '\t' << "ray" << rayId << " : " << "omega : "
            << '\t' << IRay_[rayId].omega() << '\t' << "dir : " 
	    << '\t' << IRay_[rayId].d() << nl;
    }
    radSource_.generate(directions_);

    Info<< endl;
}

void Foam::radiation::quadrature::specularGuiding
(
	const Foam::List<Foam::vector> mainAxis,
	const Foam::label cut,
	const Foam::scalar opening
)
{
	List<List< label > > facesNeighbours(directions_.size());
	List<scalar> faceMagnitudes(directions_.size(),0.0);
	
	//first row of DO have only five neighbours, as also do the last row
	// first Element in each row has special treatment
	forAll(IRay_, rayI)
	{
		faceMagnitudes[rayI]=IRay_[rayI].omega();
	}
	
	List<label> faceA(5,0); //first face
	faceA[0]=4*nPhi_ - 1;
	faceA[1]=8*nPhi_ - 1;
	faceA[2]=4*nPhi_;
	faceA[3]=4*nPhi_ + 1;
	faceA[4]=1;
	facesNeighbours[0]=faceA;
	
	List<label> faceB(5,0); //last face of first row
	faceB[0]=4*nPhi_ - 2;
	faceB[1]=8*nPhi_ - 2;
	faceB[2]=8*nPhi_ - 1;
	faceB[3]=4*nPhi_;
	faceB[4]=0;
	facesNeighbours[4*nPhi_ - 1]=faceB;
	
	List<label> faceC(5,0); //first face of last row
	faceC[0]= nRay_ - 4*nPhi_ - 1;
	faceC[1]= nRay_ - 1;
	faceC[2]= nRay_ - 4*nPhi_ + 1;
	faceC[3]= nRay_ - 8*nPhi_ + 1;
	faceC[4]= nRay_ - 8*nPhi_;
	facesNeighbours[nRay_ - 4*nPhi_]=faceC;
	
	List<label> faceD(5,0); //last face of last row
	faceD[0]= nRay_ - 4*nPhi_ - 2;
	faceD[1]= nRay_ - 2;
	faceD[2]= nRay_ - 4*nPhi_;
	faceD[3]= nRay_ - 8*nPhi_;
	faceD[4]= nRay_ - 4*nPhi_ - 1;
	facesNeighbours[nRay_ - 1]=faceD;
	
	for(label phiI = 1; phiI < 4*nPhi_-1; phiI++) //first and last row
	{
		List<label> faceE(5,0); //in first row
		faceE[0]= phiI - 1;
		faceE[1]= phiI + 4*nPhi_ - 1;
		faceE[2]= phiI + 4*nPhi_;
		faceE[3]= phiI + 4*nPhi_ + 1;
		faceE[4]= phiI + 1;
		facesNeighbours[phiI]=faceE;
		
		label pP = nRay_ - 4*nPhi_ + phiI;
		List<label> faceF(5,0); //in first row
		faceF[0]= pP - 4*nPhi_ - 1;
		faceF[1]= pP - 1;
		faceF[2]= pP + 1;
		faceF[3]= pP - 4*nPhi_ - + 1;
		faceF[4]= pP - 4*nPhi_;
		facesNeighbours[pP]=faceF;
	}
	//main Rows
	for(label thetaI = 1; thetaI < nTheta_ - 1; thetaI++) 
	{
		//first cell:
		label globalPos = thetaI*4*nPhi_;
		List<label> faceG(8,0);
		faceG[0] = globalPos - 1;
		faceG[1] = faceG[0]+4*nPhi_;
		faceG[2] = faceG[1]+4*nPhi_;
		faceG[3] = globalPos + 4*nPhi_;
		faceG[4] = faceG[3] + 1;
		faceG[5] = faceG[4]-4*nPhi_;
		faceG[6] = faceG[5]-4*nPhi_;
		faceG[7] = faceG[3]-8*nPhi_;
		facesNeighbours[globalPos]=faceG;
		
		//last cell:
		globalPos = (thetaI+1)*4*nPhi_-1;
		List<label> faceH(8,0);
		faceH[0] = globalPos - 4*nPhi_ - 1;
		faceH[1] = faceH[0]+4*nPhi_;
		faceH[2] = faceH[1]+4*nPhi_;
		faceH[3] = globalPos + 4*nPhi_;
		faceH[4] = (thetaI + 1)*4*nPhi_;
		faceH[5] = faceH[4]-4*nPhi_;
		faceH[6] = faceH[5]-4*nPhi_;
		faceH[7] = faceH[3]-8*nPhi_;
		facesNeighbours[globalPos]=faceH;
		
		//remaining cells in row
		for(label phiI = 1; phiI < 4*nPhi_-1; phiI++)
		{
			globalPos = thetaI*4*nPhi_+phiI;
			List<label> faceI(8,0);
			faceI[0] = globalPos - 4*nPhi_ - 1;
			faceI[1] = faceI[0]+4*nPhi_;
			faceI[2] = faceI[1]+4*nPhi_;
			faceI[3] = globalPos + 4*nPhi_;
			faceI[4] = faceI[3] + 1;
			faceI[5] = faceI[4]-4*nPhi_;
			faceI[6] = faceI[5]-4*nPhi_;
			faceI[7] = faceI[3]-8*nPhi_;
			facesNeighbours[globalPos]=faceI;
		}
	}
	
	//Generate a list of quadrature vertex and faces
	label vertexNumber = 2+(nTheta_-1)*4*nPhi_;
	List<vector> vertex(vertexNumber);
	List<face> faces(0);
	vertex[0] = mainAxis[2]; 						//first Vertex pointing to +z
	
	for(label posPhi = 0; posPhi < 4*nPhi_; posPhi++) //first Row of triangled faces
	{
		labelList currentFace(3);
		currentFace[0] = 0;
		currentFace[1] = posPhi + 2;
		currentFace[2] = posPhi + 1;
		if(posPhi==4*nPhi_-1)
		{ currentFace[1]=1; }
		faces.append((face)(currentFace));
	}
	for(label posTheta = 1; posTheta < nTheta_; posTheta++)
	{
		scalar deltaTheta = pi/nTheta_;
		scalar deltaPhi = pi*0.5/nPhi_;
		if(opening > 0.0)
		{
			if(posTheta > cut)
			{
				deltaTheta = (pi - opening)/(nTheta_ - cut);
			}
			else
			{
				deltaTheta = opening / cut;
			}
		}
		for(label posPhi = 0; posPhi < 4*nPhi_; posPhi++)
		{
			label globalPos = 1 + posPhi + 4*nPhi_*(posTheta-1);
			scalar thetai = posTheta*deltaTheta;
			if(opening > 0.0)
			{
				if(posTheta > cut)
				{
					thetai = opening + (posTheta-cut)*deltaTheta;
				} //if posTheta <= cut, it will already be posTheta*deltaTheta
			}
			scalar sinTheta = Foam::sin(thetai);
			scalar cosTheta = Foam::cos(thetai);
			scalar sinPhi = Foam::sin(posPhi*deltaPhi);
			scalar cosPhi = Foam::cos(posPhi*deltaPhi);
			vector currentV = (sinTheta*sinPhi*mainAxis[0])+ (sinTheta*cosPhi*mainAxis[1]) + (cosTheta*mainAxis[2]);
			vertex[globalPos] = currentV;
			if(posTheta < nTheta_ -1) //intermediate rows of squared faces
			{
				labelList currentFace(4);
				currentFace[0] = 4*nPhi_*posTheta+posPhi+1;
				currentFace[1] = 4*nPhi_*(posTheta-1)+posPhi+1;
				currentFace[2] = currentFace[1]+1;
				if(posPhi == 4*nPhi_ - 1)
				{
					currentFace[2] = 4*nPhi_*(posTheta-1)+1;
				}
				currentFace[3] = currentFace[2]+4*nPhi_;
				faces.append((face)(currentFace));
			}
			else //last row of triangled faces
			{
				labelList currentFace(3);
				currentFace[0] = 1+posTheta*4*nPhi_;
				currentFace[1] = 4*nPhi_*(posTheta-1)+posPhi+1;
				currentFace[2] = currentFace[1]+1;
				if(posPhi == 4*nPhi_ - 1)
				{
					currentFace[2] = 4*nPhi_*(posTheta-1)+1;
				}
				faces.append((face)(currentFace));
			}
		}
	}
	vertex[1+(nTheta_-1)*4*nPhi_] = -mainAxis[2]; 	//last Vertex pointing to -z
	
	specGuiders_.setSize(lastSpecked_);
	forAll(specPositions_, patchId)
	{
		if(specPositions_[patchId]>-1)
		{
			specGuiders_.set
			(
				specPositions_[patchId],
				new specularGuide
				(
					mesh_.boundary()[patchId],
					directions_,
					facesNeighbours,
					faceMagnitudes,
					vertex,
					faces,
					specOrders_[patchId],
					nLambda_
				)
			);
		}
	}
}

void Foam::radiation::quadrature::updateSpecular(const Foam::label lambdaI)
{
	forAll(specPositions_, patchI)
	{
		if(specPositions_[patchI]>-1)
		{
			PtrList<scalarField> patchIntensity(nRay_);
			forAll(IRay_, rayI)
			{
				patchIntensity.set
				(
					rayI,
					IRay_[rayI].ILambda(lambdaI).boundaryField()[patchI]
				);
			}
			specGuiders_[specPositions_[patchI]].updateSpecular(patchIntensity, lambdaI);
		}
	}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::quadrature::quadrature(const DORT& dort, const fvMesh& mesh, const label modelNumber, radSource& radSource, const bool serialQuad):

    DORT_(dort),
    modelNumber_(modelNumber),
    mesh_(mesh),
	serialQuad_(serialQuad),
	actualBand_(0),
    Glambda_(dort.nLambda()),
    QrLambda_(dort.nLambda()),
    futQrLambda_(dort.nLambda()),
    convQrLambda_(dort.nLambda()),
    preCalcQrLambda_(dort.nLambda()),
    QinLambda_(dort.nLambda()),
    futQinLambda_(dort.nLambda()),
    convQinLambda_(dort.nLambda()),
    preCalcQinLambda_(dort.nLambda()),
    nTheta_(dort.nTheta()),
    nPhi_(dort.nPhi()),
    nRay_(0),
    nLambda_(dort.nLambda()),
    IRay_(0),
    convergence_(dort.convergence()),
    maxIter_(dort.maxIter()),
    omegaMax_(0),
    sctrMatrix(0),
    scatter_(false),
    radSource_(radSource),
	directions_(0),
	specPositions_(mesh.boundary().size(),-1),
	specOrders_(mesh.boundary().size(),false),
	specGuiders_(0),
	lastSpecked_(0)
{
    //initialise();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::quadrature::~quadrature()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::radiation::quadrature::calculateBand()
{
	List<bool> rayIdConv(nRay_, false);
	if(scatter_)
		updateBandScatter(actualBand_);
	updateSpecular(actualBand_);
	scalar maxResidual = 0.0;
    label radIter = 1;
    Info<< "Radiation solver iter: " << radIter << endl;
	futQrLambda_[actualBand_] = dimensionedScalar("zero",dimMass/pow3(dimTime), 0.0);
	futQinLambda_[actualBand_] = dimensionedScalar("zero",dimMass/pow3(dimTime), 0.0);
	convQrLambda_[actualBand_] = dimensionedScalar("zero",dimMass/pow3(dimTime), 0.0);
	convQinLambda_[actualBand_] = dimensionedScalar("zero",dimMass/pow3(dimTime), 0.0);
	forAll(IRay_, rayI)
    {
		preCalcQrLambda_[actualBand_]=futQrLambda_[actualBand_];
		preCalcQinLambda_[actualBand_]=futQinLambda_[actualBand_];
		maxResidual = max(maxResidual, IRay_[rayI].correct());	//Solves ray
		if(maxRayResidual<convergence_)					//if ray is converged
        {
            rayIdConv[rayI]=true;						//Avoids future calculations over converged ray										
			//Stores fluxes of converged rays
			convQrLambda_[actualBand_]+=futQrLambda_[actualBand_]-preCalcQrLambda_[actualBand_];
			convQinLambda_[actualBand_]+=futQinLambda_[actualBand_]-preCalcQinLambda_[actualBand_];
			//Avoids duplication of converged
			futQrLambda_[actualBand_] = preCalcQrLambda_[actualBand_];
			futQinLambda_[actualBand_] = preCalcQinLambda_[actualBand_];
        }
    }
	QinLambda_[actualBand_]=futQinLambda_[actualBand_]+convQinLambda_[actualBand_];
	if (maxResidual < convergence_)						//All rays converged
    {
		updateBandQrAndG(actualBand_);
        return true;									//Only one iteration till converged->quadrature converges
    }
	do
    {
        Info<< "Radiation solver iter: " << radIter << endl;
        radIter++;
        maxResidual = 0.0;
		cleanFutures();
        forAll(IRay_, rayI)
        {
            if (!rayIdConv[rayI])
            {
				preCalcQrLambda_[actualBand_]=futQrLambda_[actualBand_];
				preCalcQinLambda_[actualBand_]=futQinLambda_[actualBand_];
				maxResidual = max(maxResidual, IRay_[rayI].correct());	//Solves ray
				if(maxRayResidual<convergence_)					//if ray is converged
				{
					rayIdConv[rayI]=true;						//Avoids future calculations over converged ray										
					//Stores fluxes of converged rays
					convQrLambda_[actualBand_]+=futQrLambda_[actualBand_]-preCalcQrLambda_[actualBand_];
					convQinLambda_[actualBand_]+=futQinLambda_[actualBand_]-preCalcQinLambda_[actualBand_];
					//Avoids duplication of converged
					futQrLambda_[actualBand_] = preCalcQrLambda_[actualBand_];
					futQinLambda_[actualBand_] = preCalcQinLambda_[actualBand_];
				}
            }
        }
		QinLambda_[actualBand_]=futQinLambda_[actualBand_]+convQinLambda_[actualBand_];
    } while (maxResidual > convergence_ && radIter < maxIter_);
    updateBandQrAndG(actualBand_);
    return false;
}

bool Foam::radiation::quadrature::calculate()
{
    // Update fluxes
    // Set rays convergence false
    List<bool> rayIdConv(nRay_, false);
    if(scatter_)
		updateScatter();
	forAll(QrLambda_,lambdaI)
	{
		updateSpecular(lambdaI);
	}
    scalar maxResidual = 0.0;
    label radIter = 1;
    Info<< "Radiation solver iter: " << radIter << endl;
    cleanFutures();
	cleanConvergeds();
    forAll(IRay_, rayI)
    {
		forAll(QrLambda_,lambdaI)
		{
			preCalcQrLambda_[lambdaI]=futQrLambda_[lambdaI];
			preCalcQinLambda_[lambdaI]=futQinLambda_[lambdaI];
		}
        maxResidual = max(maxResidual, IRay_[rayI].correct());	//Solves ray
		if(maxRayResidual<convergence_)					//if ray is converged
        {
            rayIdConv[rayI]=true;						//Avoids future calculations over converged ray
			forAll(QrLambda_,lambdaI)
			{											
				//Stores fluxes of converged rays
				convQrLambda_[lambdaI]+=futQrLambda_[lambdaI]-preCalcQrLambda_[lambdaI];
				convQinLambda_[lambdaI]+=futQinLambda_[lambdaI]-preCalcQinLambda_[lambdaI];
				//Avoids duplication of converged
				futQrLambda_[lambdaI] = preCalcQrLambda_[lambdaI];
				futQinLambda_[lambdaI] = preCalcQinLambda_[lambdaI];
			}
        }
    }
	forAll(QrLambda_,lambdaI)
	{
		QinLambda_[lambdaI]=futQinLambda_[lambdaI]+convQinLambda_[lambdaI];
	}
    if (maxResidual < convergence_)						//All rays converged
    {
		updateQrAndG();
        return true;									//Only one iteration till converged->quadrature converges
    }
    do
    {
        Info<< "Radiation solver iter: " << radIter << endl;
        radIter++;
        maxResidual = 0.0;
		cleanFutures();
        forAll(IRay_, rayI)
        {
            if (!rayIdConv[rayI])
            {
				forAll(QrLambda_,lambdaI)
				{
					preCalcQrLambda_[lambdaI]=futQrLambda_[lambdaI];
					preCalcQinLambda_[lambdaI]=futQinLambda_[lambdaI];
				}
                scalar maxRayResidual = IRay_[rayI].correct();
                maxResidual = max(maxRayResidual, maxResidual);
                if(maxRayResidual<convergence_)
                {
                    rayIdConv[rayI]=true;
					forAll(QrLambda_,lambdaI)
					{											
						//Stores fluxes of converged rays
						convQrLambda_[lambdaI]+=futQrLambda_[lambdaI]-preCalcQrLambda_[lambdaI];
						convQinLambda_[lambdaI]+=futQinLambda_[lambdaI]-preCalcQinLambda_[lambdaI];
						//Avoids duplication of converged
						futQrLambda_[lambdaI] = preCalcQrLambda_[lambdaI];
						futQinLambda_[lambdaI] = preCalcQinLambda_[lambdaI];
					}
                }
            }
        }
		forAll(QrLambda_,lambdaI)
		{
			QinLambda_[lambdaI]=futQinLambda_[lambdaI]+convQinLambda_[lambdaI];
		}
    } while (maxResidual > convergence_ && radIter < maxIter_);
    updateQrAndG();
    return false;
}


void Foam::radiation::quadrature::updateQrAndG()
{
	forAll(Glambda_,lambdaI)
	{
		updateBandQrAndG(lambdaI);
	}
}

void Foam::radiation::quadrature::updateBandQrAndG(const label lambdaI)
{
	Glambda_[lambdaI] = dimensionedScalar("zero",dimMass/pow3(dimTime), 0.0);
	label discOrdBand = lambdaI;
	if(serialQuad_)		//if serial Mode, discrete Ordinates have only one wavelenght band
	{discOrdBand=0;}
	forAll(IRay_, rayI)
	{
        Glambda_[lambdaI] += IRay_[rayI].ILambda(discOrdBand)*IRay_[rayI].omega();
	}
	QrLambda_[lambdaI]=futQrLambda_[lambdaI]+convQrLambda_[lambdaI];
}

void Foam::radiation::quadrature::updateScatter()
{
	for(label lambdaI = 0; lambdaI<nLambda_; lambdaI++)
	{
		updateBandScatter(lambdaI);
    }
}

void Foam::radiation::quadrature::updateBandScatter(const label lambdaI)
{
    forAll(IRay_, rayI)
    {
	    volScalarField Cumulant
	    (
		IOobject
		(
        	    "Cumulant" + name(modelNumber_),
        	    mesh_.time().timeName(),
        	    mesh_,
        	    IOobject::NO_READ,
        	    IOobject::NO_WRITE
        	),
        	mesh_,
	        dimensionedScalar("Cumulant", dimMass/pow3(dimTime), 0.0)
	    );
	    Cumulant = (dimensionedScalar("Cumulant", dimMass/pow3(dimTime), 0.0));
		forAll(IRay_, rayCross)
	    {
		Cumulant += DORT_.phaseFunc().sctrFactor(modelNumber_,lambdaI,rayCross,rayI)*IRay_[rayCross].ILambda(lambdaI)*IRay_[rayCross].omega();
	    }
		label discOrdBand = lambdaI;
		if(serialQuad_)			//if serial Mode, discrete Ordinates have only one wavelenght band
		{
			discOrdBand=0;
		}
	    if(DORT_.constantExt())	//if extintion is homogeneus, sigma is included in inScatter term
		{
			IRay_[rayI].setBandScatter(discOrdBand, Cumulant*max(DORT_.sigmaLambda(lambdaI)));
		}
		else
		{
			IRay_[rayI].setBandScatter(discOrdBand, Cumulant);
		}
    }
}


void Foam::radiation::quadrature::setRayIdLambdaId
(
    const word& name,
    label& rayId,
    label& lambdaId
) const
{
    // assuming name is in the form: CHARS_quadId_rayId_lambdaId
    size_t i1 = name.find_first_of("_");
    size_t i2 = name.find_first_of("_",i1+1);
    size_t i3 = name.find_last_of("_");

    rayId = Foam::readLabel(IStringStream(name.substr(i2+1, i3-1))());
    lambdaId = Foam::readLabel(IStringStream(name.substr(i3+1, name.size()-1))());
}

void Foam::radiation::quadrature::changeBand
(
	const label lambdaI
)
{
	actualBand_=lambdaI;
	forAll(IRay_,rayI)
	{
		IRay_[rayI].changeBand(lambdaI);
	}
	
}

void Foam::radiation::quadrature::getBoundaryReport
(
	const bool specUse, const bool specOrder, const label patchId, const label lambdaId
)
{
	if(specUse)
	{
		if(specPositions_[patchId]==-1)
		{
			specPositions_[patchId]=lastSpecked_;
			specOrders_[patchId]=specOrder;
			lastSpecked_++;
		}
	}
}

Foam::scalarField Foam::radiation::quadrature::getSpec
(
	label patchId, label rayId, label lambdaId
)
{
	if(specPositions_[patchId]<0)
	{
		scalarField zeroReflected(mesh_.boundary()[patchId].size(),0.0);
		return zeroReflected;
	} //No specular, shouldn't be here
	return specGuiders_[specPositions_[patchId]].getSpec(rayId,lambdaId);
}


// ************************************************************************* //
