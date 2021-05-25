

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

plot_matplotlib plot;

int main(int argc, char* argv[]) {
 
	//ODE_test_nl(0);
	//ODE_harmonic_oscillator();
	
	ODE_quantum_harmonic_oscillator();
	//ODE_Predator_Prey();
	//quantum_harmonic_oscillator();

	//ODE_Lorenz_System();
	
	//ODE_Van_der_Pol_oscillator();

	//Leapfrog_integration();
	//test_boost();
	//contour_integral();

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
	
	if(e_plot)plot.show();

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

	y[0] = 1.0;
	y[1] = 1.0;	

	std::vector<double> dydx(y.size());

	//y" = -y or y" + y = 0
	// (this is clearly equivalent to y" = −y, after introducing an auxiliary variable) :
	// y' = z
	// z' = -y

	mxws mxws;
	auto func = [&](const auto& x, const auto& y) {
		const double w0 = 2 * pi * 0.25;
		const double zeta = 0.3;

		dydx[0] = 4 * y[1];
		//dydt[1] = -2. * zeta * w0 * y[1] - pow(w0, 2) * y[0];
		dydx[1] = -4 * y[0];
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

size_t ODE_quantum_harmonic_oscillator()
{
	
	double tmin = -5.5 * pi;
	double tmax = 5.5 * pi;
	double h = 0.001;
	
	auto x = tmin;

	size_t n = 115;

	std::vector<double> y(2);

	y[0] = pow(-1, n) / (2*n*pi) / 378.5; //???
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

		dydx[0] = -y[1];

		dydx[1] = (2 * n + 1 - x * x) * y[0];

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
		Y1.push_back(y[1]);
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

		dydx[0] =  y[0] * (a - b * y[1]);
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
	double tmin = -5.5*pi;
	double tmax = 5.5*pi;
	double h = 0.0001;

	auto x = tmin;

	std::vector<double> X, Y0, Y1;

	size_t steps = 0;

	int n = 115;

	while (x <= tmax)
	{

		X.push_back(x);
		
		auto p1 = ps::Hermite_function(n, x);
		Y0.push_back(p1 * p1);
		
		//	auto p2 = 
		//	(1.0 / sqrt(sqrt(pi) * pow(2, n) * factorial<double>(n))) * ps::Hermite(n, x) * exp(-(x * x / 2.0));

		//Y1.push_back(p2 * p2);

		x += h;

		steps++;
	}

	plot.plot_somedata(X, Y0, "k", "y[0], m = "+std::to_string(n)+" ", "blue");
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

