
#include "MULTICOMPLEX.hpp"

#include "rkf78.hpp"
#include "vector_calculus.hpp"
#include <matplotlib.hpp>
#include <Embedded_Verner_6_7.hpp>
#include <Embedded_Verner_7_8.hpp>
#include <Embedded_Verner_8_9.hpp>
#include <Embedded_Fehlberg_7_8.hpp>
#include <Embedded_Fehlberg_3_4.hpp>
#include <fehlberg_4_5.hpp>
#include <embedded_fehlberg_5_6.hpp>
#include <Euler_method.hpp>
#include <Midpoint_method.hpp>

#include <Leapfrog_integration.hpp>
#include <locale>
#include <codecvt>

void Leapfrog_integration();

size_t ODE_test_nl(bool e_plot);
size_t ODE_Lorenz_System();
size_t ODE_harmonic_oscillator();
size_t quantum_harmonic_oscillator();
size_t ODE_Van_der_Pol_oscillator();
size_t ODE_quantum_harmonic_oscillator();
size_t ODE_Predator_Prey();
size_t ODE_Finite_potential_well();
size_t ODE_quantum_harmonic_oscillator_complex();

void trapezoidal();
template <class F, class T>
inline std::pair<T, T> brent_find_minima(F f, T min, T max, int digits);

plot_matplotlib plot;

int main(int argc, char* argv[]) {

	//trapezoidal();
	
	//ODE_test_nl(0);
	//ODE_harmonic_oscillator();

	//ODE_quantum_harmonic_oscillator();
	//ODE_quantum_harmonic_oscillator_complex();
	ODE_Finite_potential_well();
	
	//ODE_Predator_Prey();
	//quantum_harmonic_oscillator();

	//ODE_Lorenz_System();

	//ODE_Van_der_Pol_oscillator();

	//Leapfrog_integration();

	return 0;
}

size_t ODE_test_nl(bool e_plot)
{
	double x = -4.0;
	double tmax = 5.0;
	double h = .05;

	std::vector<double> y(1);

	y[0] = 4.0;

	std::vector<double> dydx(y.size());

	auto func = [&](const auto& x, const auto& y) {

		dydx[0] = x * sqrt(abs(y[0])) + pow(sin(x * pi / 2), 3) - 5. * (x > 2.);
		return dydx; };

	std::vector<double> X = { x }, Y0 = { y[0] };

	size_t steps = 0;
	while (x <= tmax)
	{
		Embedded_Fehlberg_7_8(func, x, y, h);
		x += h;
		X.push_back(x);
		Y0.push_back(y[0]);
		steps++;
	}

	plot.plot_somedata(X, Y0, "k", "test", "blue");

	std::u32string title = U"test";
	std::wstring_convert<std::codecvt_utf8<char32_t>, char32_t> cv;
	plot.set_title(cv.to_bytes(title));

	if (e_plot)plot.show();

	std::cout.setf(std::ios::fixed, std::ios::floatfield);
	std::cout.precision(15);
	auto p = std::minmax_element(begin(Y0), end(Y0));
	std::cout << "minY1 = " << *p.first << ", maxY1 = " << *p.second << '\n';

	std::cout << "test steps = " << steps << '\n';
	return steps;
}


size_t ODE_harmonic_oscillator()
{
	double x = 0;
	double tmax = 5;
	double h = .005;

	std::vector<double> y(2);

	y[0] = 0.5*sqrt(2);
	y[1] = 0.5*sqrt(2);

	std::vector<double> dydx(y.size());

	//y" = -y or y" + y = 0
	// (this is clearly equivalent to y" = −y, after introducing an auxiliary variable) :
	// y' = z
	// z' = -y

	//mxws mxws;
	auto func = [&](const auto& x, const auto& y) {
		const double w0 = 2 * pi * 0.25;
		const double zeta = 0.3;//0.3
		int n = 5;

		dydx[0] =   n * y[1];
		//dydx[1] = -2. * zeta * w0 * y[1] - pow(w0, 2) * y[0];
		dydx[1] = - n * y[0];

		return dydx; };

	std::vector<double> X = { x }, Y0 = { y[0] }, Y1 = { y[1] };

	size_t steps = 0;
	while (x <= tmax)
	{
		Embedded_Fehlberg_7_8(func, x, y, h);
		x += h;
		X.push_back(x);
		Y0.push_back(y[0]);
		Y1.push_back(y[1]);
		steps++;
	}

	plot.plot_somedata(X, Y0, "k", "cosine", "blue");
	plot.plot_somedata(X, Y1, "k", "sine", "red");

	std::u32string title = U"ÿ = -y";
	std::wstring_convert<std::codecvt_utf8<char32_t>, char32_t> cv;
	plot.set_title(cv.to_bytes(title));
	plot.show();

	std::cout.setf(std::ios::fixed, std::ios::floatfield);
	std::cout.precision(15);
	auto p = std::minmax_element(begin(Y0), end(Y0));
	std::cout << "minY0 = " << *p.first << ", maxY0 = " << *p.second << '\n';
	p = std::minmax_element(begin(Y1), end(Y1));
	std::cout << "minY1 = " << *p.first << ", maxY1 = " << *p.second << '\n';

	return steps;
}


size_t ODE_Finite_potential_well()
{

	double tmin = -2;
	double tmax = 2;
	double h = 0.004;

	auto x = tmin;

	std::vector<double> y(2);

	y[0] = 0; 
	y[1] = 1;

	std::vector<double> state(y.size());
	double Vo = 20, E = 19.082127;

	auto V = [&](const auto& x)
	{
		double L = 1;
		if (abs(x) > L)
			return 0.;
		else
			return Vo;
	};

	auto func = [&](const auto& x, const auto& psi) {

		state[0] = psi[1];

		state[1] = -2.0 * (V(x) - E) * psi[0];

		return state; };

	std::vector<double> X = { x }, Y0 = { y[0] }, Y1 = { y[1] };

	size_t steps = 0;


	while (x <= tmax)
	{
		//Embedded_Verner_8_9(func, x, y, h);
		Embedded_Fehlberg_7_8(func, x, y, h);
		//fehlberg_4_5(func, x, y, h);

		//Embedded_Fehlberg_3_4(func, x, y, h);
		//Embedded_Fehlberg_5_6(func, x, y, h);

		//Midpoint_method_explicit(func, x, y, h);

		x += h;
		X.push_back(x);

		Y0.push_back(y[0]);
		Y1.push_back(y[1]);
		steps++;
	}

	plot.plot_somedata(X, Y0, "k", "Y[0]", "red");
	//plot.plot_somedata(X, Y1, "k", "Y[1]", "blue");

	std::u32string title = U"Hermite functions Ψn̈(x) + (2n + 1 - x²) Ψn(x) = 0";
	std::wstring_convert<std::codecvt_utf8<char32_t>, char32_t> cv;
	plot.set_title(cv.to_bytes(title));
	plot.grid_on();
	plot.show();

	plot.plot_somedata(Y0, Y1, "k", "Y[0] vs Y[1]", "green");
	plot.show();

	std::cout.setf(std::ios::fixed, std::ios::floatfield);
	std::cout.precision(15);
	auto p = std::minmax_element(begin(Y0), end(Y0));
	std::cout << "minY0 = " << *p.first << ", maxY0 = " << *p.second << '\n';
	p = std::minmax_element(begin(Y1), end(Y1));
	std::cout << "minY1 = " << *p.first << ", maxY1 = " << *p.second << '\n';

	return steps;
}

size_t ODE_quantum_harmonic_oscillator()
{

	double tmin = -5.5 * pi;
	double tmax = 5.5 * pi;
	double h = 0.001;

	auto x = tmin;

	size_t n = 115;

	std::vector<double> y(2);

	y[0] = pow(-1, n) / (2 * n * pi) / 378.5; //???
	y[1] = 0;


	/*
double tmin = -11;
double tmax = 11;
double h = 0.01;

auto x = tmin;

size_t n = 38;

std::vector<double> y(2);

y[0] = pow(-1, n) * 1e-2;

y[1] = 0;
	*/
	std::vector<double> dydx(y.size());

	auto func = [&](const auto& x, const auto& y) {

		dydx[0] = y[1];

		dydx[1] = - (2 * n + 1 - x * x) * y[0];

		return dydx; };

	std::vector<double> X = { x }, Y0 = { y[0] }, Y1 = { y[1] };

	size_t steps = 0;

	while (x <= tmax)
	{
		//Embedded_Verner_8_9(func, x, y, h);
		Embedded_Fehlberg_7_8(func, x, y, h);
		//fehlberg_4_5(func, x, y, h);

		//Embedded_Fehlberg_3_4(func, x, y, h);
		//Embedded_Fehlberg_5_6(func, x, y, h);

		//Midpoint_method_explicit(func, x, y, h);

		x += h;
		X.push_back(x);

		Y0.push_back(y[0]);
		Y1.push_back(-y[1]);
		steps++;
	}

	plot.plot_somedata(X, Y0, "k", "Y[0], n = " + std::to_string(n) + "", "red");
	plot.plot_somedata(X, Y1, "k", "Y[1]", "blue");

	std::u32string title = U"Hermite functions Ψn̈(x) + (2n + 1 - x²) Ψn(x) = 0";
	std::wstring_convert<std::codecvt_utf8<char32_t>, char32_t> cv;
	plot.set_title(cv.to_bytes(title));
	plot.grid_on();
	plot.show();

	plot.plot_somedata(Y0, Y1, "k", "Y[0] vs Y[1]", "green");
	plot.show();

	std::cout.setf(std::ios::fixed, std::ios::floatfield);
	std::cout.precision(15);
	auto p = std::minmax_element(begin(Y0), end(Y0));
	std::cout << "minY0 = " << *p.first << ", maxY0 = " << *p.second << '\n';
	p = std::minmax_element(begin(Y1), end(Y1));
	std::cout << "minY1 = " << *p.first << ", maxY1 = " << *p.second << '\n';

	return steps;
}

size_t ODE_quantum_harmonic_oscillator_complex()
{

	double tmin = -5.5 * pi;
	double tmax = 5.5 * pi;
	double h = 0.001;

	auto x = tmin;

	size_t n = 115;

	std::vector<MX0> y(2);

	y[0] = pow(-1, n) / (2 * n * pi) / 378.5; //???
	y[1] = 0;

	std::vector<MX0> dydx(y.size());

	MX0 i{ 0,-1 };
	auto func = [&](const auto& x, const auto& y) {

		dydx[0] = i * y[1];

		dydx[1] = i * (2 * n + 1 - x * x) * y[0];

		return dydx; };

	std::vector<double> X = { x }, Y0 = { y[0].real }, Y1 = { y[1].real }, Y2 = { y[0].real }, Y3 = { y[1].real };

	size_t steps = 0;

	while (x <= tmax)
	{
		//Embedded_Verner_8_9(func, x, y, h);
		Embedded_Fehlberg_7_8(func, x, y, h);
		//fehlberg_4_5(func, x, y, h);

		//Embedded_Fehlberg_3_4(func, x, y, h);
		//Embedded_Fehlberg_5_6(func, x, y, h);

		//Midpoint_method_explicit(func, x, y, h);

		x += h;
		X.push_back(x);

		Y0.push_back(y[0].real);
		Y1.push_back(y[0].imag);
		Y2.push_back(y[1].real);
		Y3.push_back(y[1].imag);
		steps++;
	}

	plot.plot_somedata(X, Y0, "k", "y[0].real, n = " + std::to_string(n) + "", "red");
	//plot.plot_somedata(X, Y1, "k", "y[0].imag", "blue");
	//plot.plot_somedata(X, Y2, "k", "y[1].real, n = " + std::to_string(n) + "", "red");
	plot.plot_somedata(X, Y3, "k", "y[1].imag", "blue");

	std::u32string title = U"Hermite functions Ψn̈(x) + (2n + 1 - x²) Ψn(x) = 0";
	std::wstring_convert<std::codecvt_utf8<char32_t>, char32_t> cv;
	plot.set_title(cv.to_bytes(title));
	plot.grid_on();
	plot.show();

	plot.plot_somedata(Y0, Y3, "k", "y[0].real vs y[1].imag", "green");
	plot.show();

	std::cout.setf(std::ios::fixed, std::ios::floatfield);
	std::cout.precision(15);
	auto p = std::minmax_element(begin(Y0), end(Y0));
	std::cout << "minY0 = " << *p.first << ", maxY0 = " << *p.second << '\n';
	p = std::minmax_element(begin(Y3), end(Y3));
	std::cout << "minY3 = " << *p.first << ", maxY3 = " << *p.second << '\n';

	return steps;
}

//https://sam-dolan.staff.shef.ac.uk/mas212/notebooks/ODE_Example.html
size_t ODE_Predator_Prey()
{
	double x = 0;
	double tmax = 12;
	double h = .01;

	double a = 1, b = 1, c = 1, d = 1;

	std::vector<double> y(2);

	y[0] = 1.5;
	y[1] = 1.0;

	std::vector<double> dydx(y.size());

	auto func = [&](const auto& x, const auto& y) {

		dydx[0] = y[0] * (a - b * y[1]);
		dydx[1] = -y[1] * (c - d * y[0]);

		return dydx; };

	std::vector<double> X = { x }, Y0 = { y[0] }, Y1 = { y[1] };

	size_t steps = 0;
	while (x <= tmax)
	{
		Embedded_Fehlberg_7_8(func, x, y, h);
		x += h;
		X.push_back(x);
		Y0.push_back(y[0]);
		Y1.push_back(y[1]);

		steps++;
	}

	plot.plot_somedata(X, Y0, "k", "Y[0]", "blue");
	plot.plot_somedata(X, Y1, "k", "Y[1]", "red");

	std::u32string title = U"Predator-Prey Equations dx/dt = x(a - by) dy/dt = -y(c - dx)";
	std::wstring_convert<std::codecvt_utf8<char32_t>, char32_t> cv;
	plot.set_title(cv.to_bytes(title));
	plot.show();
	plot.plot_somedata(Y0, Y1, "k", "Rabbits vs Foxes", "green");
	plot.show();

	std::cout.setf(std::ios::fixed, std::ios::floatfield);
	std::cout.precision(15);
	auto p = std::minmax_element(begin(Y0), end(Y0));
	std::cout << "minY0 = " << *p.first << ", maxY0 = " << *p.second << '\n';
	p = std::minmax_element(begin(Y1), end(Y1));
	std::cout << "minY1 = " << *p.first << ", maxY1 = " << *p.second << '\n';

	return steps;
}


size_t ODE_Van_der_Pol_oscillator()
{
	double x = -2;
	double tmax = 20;
	double h = .05;

	std::vector<double> y(2);

	y[0] = 2.0;
	y[1] = 2.0;

	std::vector<double> dydx(y.size());

	auto func = [&](const auto& x, const auto& y) {

		const double mu = 5.0;

		dydx[0] = y[1];
		dydx[1] = -y[0] - mu * y[1] * (y[0] * y[0] - 1.0);

		return dydx; };

	std::vector<double> X = { x }, Y0 = { y[0] }, Y1 = { y[1] };

	size_t steps = 0;
	while (x <= tmax)
	{
		Embedded_Fehlberg_7_8(func, x, y, h);
		x += h;
		X.push_back(x);
		Y0.push_back(y[0]);
		Y1.push_back(y[1]);

		steps++;
	}

	plot.plot_somedata(X, Y0, "k", "Y[0]", "blue");
	plot.plot_somedata(X, Y1, "k", "Y[1]", "red");

	std::u32string title = U"x′′+ β(x^2−1)x′ + x = 0";// x = y...
	std::wstring_convert<std::codecvt_utf8<char32_t>, char32_t> cv;
	plot.set_title(cv.to_bytes(title));
	plot.show();

	std::cout.setf(std::ios::fixed, std::ios::floatfield);
	std::cout.precision(15);
	auto p = std::minmax_element(begin(Y0), end(Y0));
	std::cout << "minY0 = " << *p.first << ", maxY0 = " << *p.second << '\n';
	p = std::minmax_element(begin(Y1), end(Y1));
	std::cout << "minY1 = " << *p.first << ", maxY1 = " << *p.second << '\n';

	return steps;
}

//https://www.numbercrunch.de/blog/2014/08/calculating-the-hermite-functions/
size_t quantum_harmonic_oscillator()
{
	double tmin = -5.5 * pi;
	double tmax = 5.5 * pi;
	double h = 0.0001;

	auto x = tmin;

	std::vector<double> X, Y0, Y1;

	size_t steps = 0;

	int n = 115;

	//auto p1 = (1.0 / sqrt(sqrt(pi) * pow(2, n) * factorial<double>(n)));

	while (x <= tmax)
	{

		X.push_back(x);

		auto p1 = ps::Hermite_function(n, x);
		Y0.push_back(p1 * p1);

		//	auto p2 = 
		//	p1 * ps::Hermite(n, x) * exp(-(x * x / 2.0));

		//Y1.push_back(p2 * p2);

		x += h;

		steps++;
	}

	plot.plot_somedata(X, Y0, "k", "y[0], m = " + std::to_string(n) + " ", "blue");
	//plot.plot_somedata(X, Y1, "k", "y[1]", "red");


	std::u32string title = U"ÿ - 2xẏ + 2my = 0";
	std::wstring_convert<std::codecvt_utf8<char32_t>, char32_t> cv;
	plot.set_title(cv.to_bytes(title));
	plot.show();

	std::cout.setf(std::ios::fixed, std::ios::floatfield);
	std::cout.precision(15);
	auto p = std::minmax_element(begin(Y0), end(Y0));
	std::cout << "minY0 = " << *p.first << ", maxY0 = " << *p.second << '\n';


	return steps;
}

size_t ODE_Lorenz_System()
{
	double x = 0;
	double tmax = 25.0;
	double h = .001;

	std::vector<double> y(3);

	y[0] = 10;
	y[1] = 1;
	y[2] = 1;

	std::vector<double> dydx(y.size());

	auto func = [&](const auto& x, const auto& y) {
		const double sigma = 10.0;
		const double R = 28.0;
		const double b = 8.0 / 3.0;

		dydx[0] = sigma * (y[1] - y[0]);
		dydx[1] = R * y[0] - y[1] - y[0] * y[2];
		dydx[2] = -b * y[2] + y[0] * y[1];

		return dydx;
	};

	std::vector<double> X = { x }, Y0 = { y[0] }, Y1 = { y[1] }, Y2 = { y[2] };

	size_t steps = 0;
	while (x <= tmax)
	{
		Embedded_Fehlberg_7_8(func, x, y, h);
		x += h;
		X.push_back(x);
		Y0.push_back(y[0]);
		Y1.push_back(y[1]);
		Y2.push_back(y[2]);
		steps++;
	}

	plot.plot_somedata_3D(Y0, Y1, Y2, "k", "Lorenz System", "blue");
	plot.show();

	return steps;
}


template<typename F, typename T>
T Trapezoidal_Quadrature
(
	F y, 
	const T& a, 
	const T& b, 
	const T& n
)
{
	// Grid spacing
	T h = (b - a) / n;

	// Computing sum of first and last terms
	// in above formula
	T s = y(a) + y(b);

	// Adding middle terms in above formula
	for (int i = 1; i < n; i++)
		s += 2 * y(a + i * h);

	// h/2 indicates (b-a)/2n. Multiplying h/2
	// with s.
	return (h / 2) * s;
}

//https://www.geeksforgeeks.org/trapezoidal-rule-for-approximate-value-of-definite-integral/
void trapezoidal()
{
	// Range of definite integral
	double x0 = 0;
	double xn = 1;

	// Number of grids. Higher value means
	// more accuracy
	double  n = 6;

	auto func = [](const auto& x) {
		// Declaring the function f(x) = 1/(1+x*x)
		return 1 / (1 + x * x);
	};
	//Value of integral is 0.7842
	std::cout << "Value of integral is" << std::endl 
		<< Trapezoidal_Quadrature(func, x0, xn, n) << std::endl;
}

template <class F, class T>
std::pair<T, T> brent_find_minima(F f, T min, T max, int bits, size_t& max_iter)
{
	T tolerance = static_cast<T>(ldexp(1.0, 1 - bits));
	T x;  // minima so far
	T w;  // second best point
	T v;  // previous value of w
	T u;  // most recent evaluation point
	T delta;  // The distance moved in the last step
	T delta2; // The distance moved in the step before last
	T fu, fv, fw, fx;  // function evaluations at u, v, w, x
	T mid; // midpoint of min and max
	T fract1, fract2;  // minimal relative movement in x

	static const T golden = 0.3819660f;  // golden ratio, don't need too much precision here!

	x = w = v = max;
	fw = fv = fx = f(x);
	delta2 = delta = 0;

	uintmax_t count = max_iter;

	do {
		// get midpoint
		mid = (min + max) / 2;
		// work out if we're done already:
		fract1 = tolerance * fabs(x) + tolerance / 4;
		fract2 = 2 * fract1;
		if (fabs(x - mid) <= (fract2 - (max - min) / 2))
			break;

		if (fabs(delta2) > fract1)
		{
			// try and construct a parabolic fit:
			T r = (x - w) * (fx - fv);
			T q = (x - v) * (fx - fw);
			T p = (x - v) * q - (x - w) * r;
			q = 2 * (q - r);
			if (q > 0)
				p = -p;
			q = fabs(q);
			T td = delta2;
			delta2 = delta;
			// determine whether a parabolic step is acceptable or not:
			if ((fabs(p) >= fabs(q * td / 2)) || (p <= q * (min - x)) || (p >= q * (max - x)))
			{
				// nope, try golden section instead
				delta2 = (x >= mid) ? min - x : max - x;
				delta = golden * delta2;
			}
			else
			{
				// whew, parabolic fit:
				delta = p / q;
				u = x + delta;
				if (((u - min) < fract2) || ((max - u) < fract2))
					delta = (mid - x) < 0 ? (T)-fabs(fract1) : (T)fabs(fract1);
			}
		}
		else
		{
			// golden section:
			delta2 = (x >= mid) ? min - x : max - x;
			delta = golden * delta2;
		}
		// update current position:
		u = (fabs(delta) >= fract1) ? T(x + delta) : (delta > 0 ? T(x + fabs(fract1)) : T(x - fabs(fract1)));
		fu = f(u);
		if (fu <= fx)
		{
			// good new point is an improvement!
			// update brackets:
			if (u >= x)
				min = x;
			else
				max = x;
			// update control points:
			v = w;
			w = x;
			x = u;
			fv = fw;
			fw = fx;
			fx = fu;
		}
		else
		{
			// Oh dear, point u is worse than what we have already,
			// even so it *must* be better than one of our endpoints:
			if (u < x)
				min = u;
			else
				max = u;
			if ((fu <= fw) || (w == x))
			{
				// however it is at least second best:
				v = w;
				w = u;
				fv = fw;
				fw = fu;
			}
			else if ((fu <= fv) || (v == x) || (v == w))
			{
				// third best:
				v = u;
				fv = fu;
			}
		}

	} while (--count);

	max_iter -= count;

	return std::make_pair(x, fx);
}

template <class F, class T>
inline std::pair<T, T> brent_find_minima(F f, T min, T max, int digits)
{
	auto m = (std::numeric_limits<size_t>::max)();
	return brent_find_minima(f, min, max, digits, m);
}

