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

inline Foam::label 
Foam::radiation::quadrature::modelNumber() const
{
    return modelNumber_;
}

inline const Foam::radiation::discreteOrdinate&
Foam::radiation::quadrature::IRay(const label rayI) const
{
    return  IRay_[rayI];
}


inline const Foam::volScalarField&
Foam::radiation::quadrature::IRayLambda
(
    const label rayI,
    const label lambdaI
) const
{
    return IRay_[rayI].ILambda(lambdaI);
}


inline Foam::label Foam::radiation::quadrature::nTheta() const
{
    return nTheta_;
}


inline Foam::label Foam::radiation::quadrature::nPhi() const
{
    return nPhi_;
}


inline Foam::label Foam::radiation::quadrature::nRay() const
{
    return nRay_;
}


inline Foam::label Foam::radiation::quadrature::nLambda() const
{
    return nLambda_;
}


inline const Foam::volScalarField& Foam::radiation::quadrature::Glambda
(
	const label lambdaI
) const
{
    return Glambda_[lambdaI];
}


inline const Foam::volScalarField& Foam::radiation::quadrature::QrLambda
(
	const label lambdaI
) const
{
    return QrLambda_[lambdaI];
}

inline Foam::volScalarField& Foam::radiation::quadrature::futQrLambda
(
	const label lambdaI
)
{
    return futQrLambda_[lambdaI];
}

inline const Foam::volScalarField& Foam::radiation::quadrature::QinLambda
(
	const label lambdaI
) const
{
    return QinLambda_[lambdaI];
}

inline Foam::volScalarField& Foam::radiation::quadrature::futQinLambda
(
	const label lambdaI
)
{
    return futQinLambda_[lambdaI];
}

inline Foam::scalar Foam::radiation::quadrature::omegaMax() const
{
    return omegaMax_;
}

inline Foam::scalar Foam::radiation::quadrature::sourceFlux
(const word patchID, const label face, const label rayId, const label band) const
{
    return radSource_.calculate(patchID,face,rayId,band);
}

inline const Foam::List<Foam::vector> Foam::radiation::quadrature::directions() const
{
	return directions_;
}

// ************************************************************************* //
