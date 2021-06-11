
#include "MULTICOMPLEX.hpp"
#include "rkf78.hpp"
#include "vector_calculus.hpp"
#include "matplotlib.hpp"
#include "Embedded_Verner_6_7.hpp"
#include "Embedded_Verner_7_8.hpp"
#include "Embedded_Verner_8_9.hpp"
#include "Embedded_Fehlberg_7_8.hpp"
#include "Embedded_Fehlberg_3_4.hpp"
#include "embedded_fehlberg_5_6.hpp"
#include "Euler_method.hpp"
#include "Midpoint_method.hpp"
#include "Leapfrog_integration.hpp"
#include "brent.hpp"
#include "dekker.hpp"
#include "secant.hpp"
#include "Laplace_transform.hpp"
#include <codecvt>

void Leapfrog_integration();
void ODE_test_nl(bool e_plot);
void ODE_Lorenz_System();
void ODE_harmonic_oscillator();
void quantum_harmonic_oscillator();
void ODE_Van_der_Pol_oscillator();
void ODE_quantum_harmonic_oscillator();
void ODE_Predator_Prey();
void ODE_Quantum_Solver(int mode = 0);
void ODE_quantum_harmonic_oscillator_complex();

void trapezoidal();

template <typename T>
int sign(const T& x);

template<typename T>
std::vector<T> linspace(const T start_in, const T end_in, std::size_t num_in);

template <typename F, typename T>
std::vector<T> Find_all_zeroes
(
	const F& Wave_function,
	const std::vector<T>& en, 
	bool groundstate
);

template <typename T>
std::vector<T> zeroCrossing(const std::vector<T>& s, const std::vector<T>& en);

plot_matplotlib plot;
std::string colours(const int& t);

int main(int argc, char* argv[]) {

	//trapezoidal();

	//ODE_test_nl(1);
	//ODE_harmonic_oscillator();

	//ODE_quantum_harmonic_oscillator();

	//ODE_quantum_harmonic_oscillator_complex();

	ODE_Quantum_Solver(4);

	//ODE_Predator_Prey();
	//quantum_harmonic_oscillator();

	//ODE_Lorenz_System();

	//ODE_Van_der_Pol_oscillator();

	//Leapfrog_integration();

	return 0;
}

void ODE_test_nl(bool e_plot)
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

	while (x <= tmax)
	{
		Embedded_Fehlberg_7_8(func, x, y, h);
		x += h;
		X.push_back(x);
		Y0.push_back(y[0]);
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
}

void ODE_harmonic_oscillator()
{
	double x = 0;
	double tmax = 5;
	double h = .005;

	std::vector<double> y(2);

	y[0] = 0;//0.5*sqrt(2);
	y[1] = 1;//0.5*sqrt(2)0

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

		dydx[0] = n * y[1];
		//dydx[1] = -2. * zeta * w0 * y[1] - pow(w0, 2) * y[0];
		dydx[1] = -n * y[0];

		return dydx; };

	std::vector<double> X = { x }, Y0 = { y[0] }, Y1 = { y[1] };

	while (x <= tmax)
	{
		Embedded_Fehlberg_7_8(func, x, y, h);
		x += h;
		X.push_back(x);
		Y0.push_back(y[0]);
		Y1.push_back(y[1]);
	}

	plot.plot_somedata(X, Y0, "k", "cosine", "blue", 1);
	plot.plot_somedata(X, Y1, "k", "sine", "red", 1);

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

}

template <typename T>
std::vector<T> gaussian_wave_packet(const std::vector<T>& en, const T& sigma = 1.0, const T& mu = 0.0)
{
	std::vector<T> v;
	T a = 1. / (sigma * sqrt(2 * pi));
	a *= 1.198585e4;

	for (auto& x : en)
	{
		v.push_back(a * exp(-0.5 * pow((x - mu) / sigma, 2)));
	}

	return v;
}

template <typename T>
std::tuple<std::vector<std::vector<T>>, std::vector<T>> ODE_Q_sine_cosine
(
	const T& Vo_b,
	const T& Vo_e,
	const T& tmin,
	const T& tmax,
	const T& h
)
{
	std::vector<T> y(4);
	std::vector<T> state(y.size());
	std::vector<std::vector<T>> Y(y.size());

	T x = tmin;
	T E = 0;
	auto en = linspace(Vo_b, Vo_e + h, int(2 * abs(Vo_e - Vo_b)));

	auto SE = [&](const auto& x, const auto& psi) {

		state[0] = psi[1];
		state[1] = -2 * E * psi[0];

		state[2] = psi[3];
		state[3] = -2 * E * psi[2];

		return state;
	};

	auto Wave_function = [&](const auto& energy) {
		E = energy;

		Y.clear();
		Y.resize(y.size());

		x = tmin;

		y[0] = 1;
		y[1] = 0;

		y[2] = 0;
		y[3] = 1;

		for (size_t i = 0; i < y.size(); i++)
			Y[i].emplace_back(y[i]);

		while (x <= tmax)
		{
			Midpoint_method_explicit(SE, x, y, h);
			//Midpoint_method_implicit(SE, x, y, h);
			//Euler_method(SE, x, y, h);
			//Embedded_Fehlberg_3_4(SE, x, y, h);
			//Embedded_Fehlberg_7_8(SE, x, y, h);

			x += h;
			for (size_t i = 0; i < y.size(); i++)
				Y[i].emplace_back(y[i]);
		}

		return *Y[0].rbegin();
	};

	std::vector<T> E_zeroes1, E_zeroes2;

	//E_zeroes1 = Find_all_zeroes(Wave_function, en);

	E_zeroes2 = Find_all_zeroes(Wave_function, en, 0);

	//std::cout << E_zeroes;

	std::vector<std::vector<T>> psi_sol;

	for (auto& E : E_zeroes1) {
		//Wave_function(E);
		//psi_sol.push_back(Y0);
	}

	for (auto& E : E_zeroes2) {
		Wave_function(E);
		psi_sol.emplace_back(Y[3]);
		std::reverse(Y[3].begin(), Y[3].end());
		psi_sol.emplace_back(Y[3]);
	}

	E_zeroes2.insert(E_zeroes2.end(), E_zeroes2.begin(), E_zeroes2.end());
	return { psi_sol, E_zeroes2 };
}

void ODE_Quantum_Solver(int mode)
{
	double tmin = -1;
	double tmax = 1;

	if (mode == 0) {
		tmin = -2; tmax = 2;
	}
	else if (mode == 1) {
		tmin = -1; tmax = 1;
	}
	else if (mode == 2) {
		tmin = -6; tmax = 6;
	}

	double h = 0.005;

	std::vector<double> y(2);

	std::vector<double> state(y.size());
	std::vector<double> X;
	std::vector<std::vector<double>> Y(y.size());

	double sigma = 1, mu = 0;
	double Vo = 20, E = 0;
	int physicist = 0;

	if (mode == 2)Vo = 10;

	bool wave_packet = 0;
	
	double L1 = -1., L2 = 1.;

	if (mode == 2) {
		physicist = 1;
	}

	if (mode == 3 || mode == 4) {
		wave_packet = true;
		sigma = 32, mu = 1;
		Vo = 0; tmin = -9; tmax = 11;
		h = 0.005;
		E = 0;
		physicist = 0;
	}

	double Vb=0;
	if (mode == 4) {
		double ws = .4;
		L1 = 2 * mu, L2 = 2 * mu + ws;
		plot.line(L1, L1, 0, 2060);
		plot.line(L2, L2, 0, 2060);
		Vb = 20;
		plot.text(L2 + 0.5, 1900, "Vb = " + std::to_string(int(Vb)) + " ", "black", 12);
		std::ostringstream out;
		out.precision(2);
		out << std::fixed << "L   = " << ws;
		plot.text(L2 + 0.5, 1800, out.str(), "black", 12);
	}

	auto x = tmin;

	auto V = [&](const auto& x)
	{
		if (mode == 4) {
			if (x > (L1 - mu) && x < (L2 - mu))
				return  Vb;
		}

		if (mode == 0 || mode == 1) {
			if (x > L1 && x < L2) {
				return 0.;
			}
			else return Vo;
		}

		else return -(2 * E + physicist - (x * x) / sigma);
	};

	auto SE = [&](const auto& x, const auto& psi) {

		state[0] = psi[1];

		if (mode == 0 || mode == 1)
			state[1] = 2 * (V(x) - E) * psi[0];

		else state[1] = V(x - mu) * psi[0];

		return state;
	};

	X.emplace_back(x);
	while (x <= tmax)
	{
		x += h;
		X.emplace_back(x);
	}

	auto Wave_function = [&](const auto& energy) {
		E = energy;

		x = tmin;

		y[0] = 0;
		y[1] = 1;

		Y.clear();
		Y.resize(y.size());

		Y[0].emplace_back(y[0]);
		Y[1].emplace_back(y[1]);

		while (x <= tmax)
		{
			Midpoint_method_explicit(SE, x, y, h);
			//Midpoint_method_implicit(SE, x, y, h);
			//Euler_method(SE, x, y, h);
			//Embedded_Fehlberg_3_4(SE, x, y, h);
			//Embedded_Fehlberg_7_8(SE, x, y, h);

			x += h;
			Y[0].emplace_back(y[0]);
			Y[1].emplace_back(y[1]);

		}
		return *Y[0].rbegin();
	};

	std::cout.setf(std::ios::fixed, std::ios::floatfield);
	std::cout.precision(8);

	int t = 0;
	std::ostringstream oss;
	oss.setf(ios::fixed);
	oss.precision(2);

	std::vector<std::vector<double>> psi_sola, psi_solb, psi_sols;
	std::vector<double> E_zeroes, en;

	if (Vo == 0) {
		en = linspace(0., 1 / sqrt(sigma), 2);
	}
	else en = linspace(E, Vo, int(2 * Vo));

	E_zeroes = Find_all_zeroes(Wave_function, en, mode == 4);
	if (E_zeroes.empty()) { std::cout << "No roots found !\n\n"; return; }

	for (auto& E : E_zeroes) {
		Wave_function(E);

		psi_sola.emplace_back(Y[0]);
		psi_solb.emplace_back(Y[1]);
		std::vector<double> Ys;
		for (auto& k : Y[0])
			Ys.emplace_back(k * k);
		psi_sols.emplace_back(Ys);
		t++;
	}

	if (wave_packet)
	{
		auto gwp = gaussian_wave_packet(X, (1 + 2 * 0.0072973525693) / (1 + sqrt(2)) * sqrt(sigma), mu);//σ μ
	//gwp -= gaussian_wave_packet(X, 1 / (1 + sqrt(2)) * sqrt(sigma + 1 + 2 * 0.0072973525693), mu);
		std::u32string text = U"Gaussian wave packet( (1 + 2 * 0.0072973525693) / (1 + sqrt(2)) * sqrt(32), 1 )\\n\
Normal distribution(σ = 2.377343271833, μ = 1)\\n\
0.0072973525693 = Fine-structure constant\\n\
2.377343271833 ≈ 4 sqrt(C_HSM), C_HSM = Hafner-Sarnak-McCurley Constant";//4 sqrt(C_HSM)≈2.3773476711
	//http://www.totemconsulting.ca/FineStructure.html
	//https://mathworld.wolfram.com/Hafner-Sarnak-McCurleyConstant.html

		std::wstring_convert<std::codecvt_utf8<char32_t>, char32_t> cv1;
		text = U"Normal distribution (σ = (1 + 2 * 0.0072973525693) / (1 + sqrt(2)) * sqrt(32), μ = 1 )";
		plot.plot_somedata(X, gwp, "k", cv1.to_bytes(text), "Red", 1.0);

		t = 0;
		auto v = ODE_Q_sine_cosine(0., 8., tmin, tmax, h);
		for (auto& i : std::get<0>(v))
		{
			oss.str(std::string());
			oss << std::get<1>(v)[t];
			//std::cout << oss.str() << std::endl;
			//auto k = zeroCrossing(i, X);
			//std::cout << k << std::endl;
			//auto s = i * 1000;
			//plot.plot_somedata(X, s, "k", "E = " + oss.str() + " ", colours(t), 1.0);
			t++;
		}

		auto rf = psi_sola[0] * *(std::get<0>(v).rbegin());
		plot.plot_somedata(X, rf, "k", "E = " + oss.str() + " ", "Red", 1.0);
		rf = psi_sola[0] * *(std::get<0>(v).rbegin() + 1);
		plot.plot_somedata(X, rf, "k", "E = " + oss.str() + " ", "Green", 1.0);
	}

	t = 0;
	for (auto& E : E_zeroes) {
		Wave_function(E);
		oss.str(std::string());
		oss << E;
		std::cout << "E " << E << std::endl;
		plot.plot_somedata(X, psi_sola[t], "k", "E = " + oss.str() + " ", colours(t), 1.0);
		//plot.plot_somedata(X, psi_solb[t], "k", "E = " + oss.str() + " ", colours(t), 1.0);
		t++;
	}

	std::u32string title;
	if (mode == 0)title = U"Finite potential well";
	else if (mode == 1)title = U"Infinite potential well = Particle in 1-D Box";
	else if (mode == 2)title = U"Quantum harmonic oscillator";
	else if (mode == 3)title = U"Quantum Gaussian wave packet";
	else title = U"Quantum Gaussian wave packet + tunnelling";
	
	std::wstring_convert<std::codecvt_utf8<char32_t>, char32_t> cv;
	plot.set_title(cv.to_bytes(title));
	plot.grid_on();
	plot.show();

	t = 0;
	for (auto& E : E_zeroes) {
		Wave_function(E);
		oss.str(std::string());
		oss << E;
		plot.plot_somedata(Y[1], Y[0], "k", "E = " + oss.str() + " ", colours(t++), 1.0);
	}
	plot.set_title(cv.to_bytes(title));
	plot.show();

	std::cout.setf(std::ios::fixed, std::ios::floatfield);
	std::cout.precision(15);
	auto p = std::minmax_element(begin(Y[0]), end(Y[0]));
	std::cout << "minY0 = " << *p.first << ", maxY0 = " << *p.second << '\n';
	p = std::minmax_element(begin(Y[1]), end(Y[1]));
	std::cout << "minY1 = " << *p.first << ", maxY1 = " << *p.second << '\n';
}

void ODE_quantum_harmonic_oscillator()
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

		dydx[1] = -(2 * n + 1 - x * x) * y[0];

		return dydx; };

	std::vector<double> X = { x }, Y0 = { y[0] }, Y1 = { y[1] };

	while (x <= tmax)
	{
		//Embedded_Verner_8_9(func, x, y, h);
		Embedded_Fehlberg_7_8(func, x, y, h);
		//Embedded_Fehlberg_3_4(func, x, y, h);
		//Embedded_Fehlberg_5_6(func, x, y, h);
		//Midpoint_method_explicit(func, x, y, h);

		x += h;
		X.push_back(x);

		Y0.push_back(y[0] * y[0]);
		Y1.push_back(-y[1]);
	}

	plot.plot_somedata(X, Y0, "k", "Y[0], n = " + std::to_string(n) + "", "red");
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

}

void ODE_quantum_harmonic_oscillator_complex()
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

	while (x <= tmax)
	{
		//Embedded_Verner_8_9(func, x, y, h);
		Embedded_Fehlberg_7_8(func, x, y, h);
		//Embedded_Fehlberg_3_4(func, x, y, h);
		//Embedded_Fehlberg_5_6(func, x, y, h);
		//Midpoint_method_explicit(func, x, y, h);

		x += h;
		X.push_back(x);

		Y0.push_back(y[0].real);
		Y1.push_back(y[0].imag);
		Y2.push_back(y[1].real);
		Y3.push_back(y[1].imag);
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
}

//https://sam-dolan.staff.shef.ac.uk/mas212/notebooks/ODE_Example.html
void ODE_Predator_Prey()
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

	while (x <= tmax)
	{
		Embedded_Fehlberg_7_8(func, x, y, h);
		x += h;
		X.push_back(x);
		Y0.push_back(y[0]);
		Y1.push_back(y[1]);
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
}


void ODE_Van_der_Pol_oscillator()
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

	while (x <= tmax)
	{
		Embedded_Fehlberg_7_8(func, x, y, h);
		x += h;
		X.push_back(x);
		Y0.push_back(y[0]);
		Y1.push_back(y[1]);
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
}

//https://www.numbercrunch.de/blog/2014/08/calculating-the-hermite-functions/
void quantum_harmonic_oscillator()
{
	double tmin = -5.5 * pi;
	double tmax = 5.5 * pi;
	double h = 0.0001;

	auto x = tmin;

	std::vector<double> X, Y0, Y1;


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

}

void ODE_Lorenz_System()
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

	while (x <= tmax)
	{
		Embedded_Fehlberg_7_8(func, x, y, h);
		x += h;
		X.push_back(x);
		Y0.push_back(y[0]);
		Y1.push_back(y[1]);
		Y2.push_back(y[2]);
	}

	plot.plot_somedata_3D(Y0, Y1, Y2, "k", "Lorenz System", "blue");
	plot.show();
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

template <typename T>
int sign(const T& x)
{
	if (x > 0)return 1;
	else if (x < 0) return -1;
	else return 0;
}

template <typename T>
std::vector<T> linspace(const T start_in, const T end_in, std::size_t num_in)
{
	std::vector<T> linspaced(num_in);

	T start = start_in;
	T end = end_in;
	T num = T(num_in);

	if (num == 0) { return linspaced; }
	if (num == 1)
	{
		linspaced[0] = (start);
		return linspaced;
	}

	T delta = (end - start) / (num - 1);

	for (size_t i = 0; i < num - 1; ++i)
	{
		linspaced[i] = (start + delta * i);
	}
	linspaced[num_in - 1] = end; // I want to ensure that start and end
														// are exactly the same as the input
	return linspaced;
}

template <typename F, typename T>
std::vector<T> Find_all_zeroes
(
	const F& Wave_function,
	const std::vector<T>& en,
	bool groundstate
)
{
	//Gives all zeroes in y = Psi(x)
	const T epsilon = 1e-8;
	const auto brent = new Brent(epsilon, Wave_function);
	const auto secant = new Secant(epsilon, Wave_function);
	const auto dekker = new Dekker(epsilon, Wave_function);

	std::vector<T> s;
	std::vector<T> all_zeroes;

	for (auto& e1 : en)
	{
		auto psi_b = Wave_function(e1);

		s.push_back(sign(psi_b));

		all_zeroes.clear();

		if (groundstate)
		{
			T zero = secant->solve(e1, e1 + epsilon);
			all_zeroes.push_back(zero);
			return all_zeroes;
		}

		else {
			//static int t = 0;
			for (size_t i = 0; i < s.size() - 1; i++)
			{
				if ((s[i] + s[i + 1]) == 0)
				{
					//T zero = secant->solve(en[i], en[i + 1]);
					//T zero = dekker->solve(en[i], en[i + 1]);
					T zero = brent->solve(en[i], en[i + 1]);
					//std::cout << zero << " " << t++ << std::endl;
					all_zeroes.push_back(zero);
				}
			}
		}
	}
	return all_zeroes;
}

std::string colours(const int& t)
{
	std::string colours[24] = { "Blue", "Green",
																"Red", "Cyan", "Magenta", "Yellow", "Black", "Silver",
		"Blue", "Green",
																"Red", "Cyan", "Magenta", "Yellow", "Black", "Silver",
		"Blue", "Green",
																"Red", "Cyan", "Magenta", "Yellow", "Black", "Silver" };
	return colours[t];
}

template <typename T>
std::vector<T> zeroCrossing(const std::vector<T>& s, const std::vector<T>& en)
{
	std::vector<size_t> zerCrossi;
	std::vector<T> zerCross;

	for (size_t i = 0; i < s.size() - 1; i++)     /* loop over data  */
	{
		if ((sign(s[i]) + sign(s[i + 1])) == 0) /* set zero crossing location */
			zerCrossi.push_back(i);
	}

	for (size_t i = 0; i < zerCrossi.size(); i++)
		zerCross.push_back(en[zerCrossi[i]]);

	return zerCross;
}