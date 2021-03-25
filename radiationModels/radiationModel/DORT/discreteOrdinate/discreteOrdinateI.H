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

inline const Foam::vector& Foam::radiation::discreteOrdinate::d() const
{
    return d_;
}


inline const Foam::vector& Foam::radiation::discreteOrdinate::dAve() const
{
    return dAve_;
}


inline Foam::scalar Foam::radiation::discreteOrdinate::nLambda() const
{
    return nLambda_;
}


inline Foam::scalar Foam::radiation::discreteOrdinate::phi() const
{
    return phi_;
}


inline Foam::scalar Foam::radiation::discreteOrdinate::theta() const
{
    return theta_;
}


inline Foam::scalar Foam::radiation::discreteOrdinate::omega() const
{
    return omega_;
}

inline const Foam::volScalarField&
Foam::radiation::discreteOrdinate::ILambda
(
    const label lambdaI
) const
{
    return ILambda_[lambdaI];
}


// ************************************************************************* //
