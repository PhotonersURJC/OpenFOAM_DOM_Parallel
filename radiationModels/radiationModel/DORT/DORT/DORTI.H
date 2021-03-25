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

inline Foam::label Foam::radiation::DORT::nTheta() const
{
    return nTheta_;
}


inline Foam::label Foam::radiation::DORT::nPhi() const
{
    return nPhi_;
}


inline Foam::label Foam::radiation::DORT::nRay() const
{
    return 8*nPhi_*nTheta_;
}


inline Foam::label Foam::radiation::DORT::nLambda() const
{
    return nLambda_;
}

inline Foam::label Foam::radiation::DORT::nQuad() const
{
    return nQuad_;
}

inline const Foam::volScalarField& Foam::radiation::DORT::kappa() const
{
    return kappa_;
}


inline const Foam::volScalarField& Foam::radiation::DORT::kappaLambda
(
    const label lambdaI
) const
{
    return kappaLambda_[lambdaI];
}

inline const Foam::volScalarField& Foam::radiation::DORT::sigma() const
{
    return sigma_;
}


inline const Foam::volScalarField& Foam::radiation::DORT::sigmaLambda
(
    const label lambdaI
) const
{
    return sigmaLambda_[lambdaI];
}

inline const Foam::volScalarField& Foam::radiation::DORT::Glambda
(
	const label lambdaI
) const
{
    return Glambda_[lambdaI];
}


inline const Foam::volScalarField& Foam::radiation::DORT::QrLambda
(
	const label lambdaI
) const
{
    return QrLambda_[lambdaI];
}

inline const Foam::volScalarField& Foam::radiation::DORT::QinLambda
(
	const label lambdaI
) const
{
    return QinLambda_[lambdaI];
}

inline const Foam::volScalarField& Foam::radiation::DORT::Elambda
(
	const label lambdaI
) const
{
    return Elambda_[lambdaI];
}

inline Foam::label Foam::radiation::DORT::maxIter() const
{
    return maxIter_;
}

inline Foam::scalar Foam::radiation::DORT::convergence() const
{
    return convergence_;
}


// ************************************************************************* //
