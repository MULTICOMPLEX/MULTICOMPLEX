#ifndef __ASSOCIATED_LEGENDRE_H__
#define __ASSOCIATED_LEGENDRE_H__

#include <vector>
#include <math.h>

//#include "LegendreExceptions.hpp"

// Uses the recurrence relation:
// (l - m) P_l^m (x) = x (2l - 1) P_l-1^m (x) - (l + m - 1) P_l-2^m (x)
// to calculate the value at x, for starting values it uses the following
// relations:
// P_m^m (x) = (-1)^m (2m - 1)!! (1 - x^2)^m/2
// P_m+1^m (x) = x (2m + 1) P_m^m (x)
// constraints here are that -1 <= x <= 1
// and 0 <= m <= l
// Calculating Associated Legendre Polynomial Values in the interval -1 <= x <= 1
// follows the algorithm:
// (1) Calculate P_m^m (x) with explicit formula as the starting point (since m <= l 
//     this will always work as a starting point)
// (2) If l == m then job is done and we can return the value calculated in (1)
// (3) Then calculate the next value P_m+1^m (x) (again since m <= l 
//     this will always work as our second point for the final recurrence formula)
// (4) If l == m+1 then we already have our value from (3) and can return
// (5) Final step we can calculate the rest of the polynomial values using the recurrence 
//     formula by iterating values of l from m+2 -> l

// For normalised values of Associated Legendre Polynomials the normalising factor is as 
// follows:
// Theta_l^m (x) = sqrt{ (2l + 1)/2 (l-m)!/(l+m)! } P_l^m (x)
// NOTE that the normalising factor for spherical harmonics (with azimuthal part) is:
// Theta_l^m (x) = sqrt{ (2l + 1)/4*pi (l-m)!/(l+m)! } P_l^m (x)

// Legendre polynomials P_l (x) can be obtained by using the Associated Legendre polynomial
// class and setting m = 0

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Associated Legendre Polynomials
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class AssociatedLegendre
{
public:
	AssociatedLegendre(unsigned int l, unsigned int m);
	~AssociatedLegendre() = default;

	unsigned int l() const;
	unsigned int m() const;

	double operator()(double x) const;
	double rf;
	MX0 i;

	unsigned int m_l;
	unsigned int m_m;

	double calculatePolynomialValue(double x) const;
	MX0 SphericalHarmonic(double theta /*cos theta*/, double phi) const;

};

inline AssociatedLegendre::AssociatedLegendre(unsigned int l, unsigned int m) : m_l(l), m_m(m)
{
	//static_assert (m_m < 0, "Upper Index (m) cannot be negative!");
	//static_assert (m_m > m_l, "Upper Index (m) must be less than Lower Index (l)");
	
	rf = std::sqrt((2 * m_l + 1) / (2 * two_pi)) * std::sqrt(Fac<double>(m_l - m_m) / Fac<double>(m_l + m_m));
	i.imag = 1;
}

inline unsigned int AssociatedLegendre::l() const
{
	return m_l;
}

inline unsigned int AssociatedLegendre::m() const
{
	return m_m;
}

inline double AssociatedLegendre::operator()(double x) const
{
	return calculatePolynomialValue(x);
}


inline double AssociatedLegendre::calculatePolynomialValue(double x) const
{
	double pmm = 1.0;
	double pmp1m;

	// Check Inputs
	//static_assert (fabs(x) > 1.0, "Input Out Of Range");

	if (m_m > 0)
	{
		// P_m^m(x) = (-1) ^ m(2m - 1)!!(1 - x ^ 2) ^ m / 2
		double sqrtomx2 = sqrt(1.0 - x*x);
		double oddInt = 1.0;
		for (int i = 1; i <= m_m; ++i)
		{
			pmm *= -1.0*oddInt*sqrtomx2;
			oddInt += 2.0;
		}
	}
	if (m_l == m_m)
	{
		return pmm;
	}
	else
	{
		// P_m + 1 ^ m(x) = x(2m + 1) P_m^m(x)
		pmp1m = x*(2.0*m_m + 1.0) * pmm;
		if (m_l == m_m + 1)
		{
			return pmp1m;
		}
		else
		{
			double pmp2m = 0.0;
			for (int i = m_m + 2; i <= m_l; ++i)
			{
				// (l - m) P_l^m (x) = x (2l - 1) P_l-1^m (x) - (l + m - 1) P_l-2^m (x)
				pmp2m = x*(2.0*m_l - 1.0)*pmp1m - (m_l + m_m - 1.0)*pmm;
				pmp2m /= (m_l - m_m);
				// Rotate variables for next iteration
				pmm = pmp1m;
				pmp1m = pmp2m;
			}
			return pmp2m;
		}
	}
}

inline MX0 AssociatedLegendre::SphericalHarmonic(double theta /*cos theta*/, double phi) const
{	
	return rf * calculatePolynomialValue(cos(theta)) * exp(i * m_m * phi);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Normalised Associated Legendre Polynomials
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class NormalisedAssociatedLegendre
{
public:
	NormalisedAssociatedLegendre(unsigned int l, unsigned int m);
	~NormalisedAssociatedLegendre() = default;

	double operator()(double x) const;
	std::vector<double> operator()(std::vector<double>& x) const;

private:
	AssociatedLegendre m_associatedLegendre;

	double calculateNormalisationValue() const;
};

inline NormalisedAssociatedLegendre::NormalisedAssociatedLegendre(unsigned int l, unsigned int m) : m_associatedLegendre(l, m)
{
}

inline double NormalisedAssociatedLegendre::operator()(double x) const
{
	return calculateNormalisationValue()*m_associatedLegendre(x);
}


inline double NormalisedAssociatedLegendre::calculateNormalisationValue() const
{
	unsigned int lmm = m_associatedLegendre.l() - m_associatedLegendre.m();
	unsigned int lpm = m_associatedLegendre.l() + m_associatedLegendre.m();
	return sqrt(((2.0*m_associatedLegendre.l() + 1.0) / 2.0)*(std::tgamma(lmm + 1) / std::tgamma(lpm + 1)));
}

#endif // __ASSOCIATED_LEGENDRE_H__
