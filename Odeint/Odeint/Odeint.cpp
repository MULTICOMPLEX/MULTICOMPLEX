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
#include "embedded_prince_dormand_v3_4_5.hpp"

#include "NumCpp.hpp"

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
void ODE_test_poly();
template <typename T>
T tildePlm(const int l, const int m, const T x);
template <typename T>
T Plm(const int l, const int m, const T x);
void tal();
template <typename T>
std::complex<T> Ylm(const int l, const int m, const T x, const T y, const T z);

// calculate the spherical harmonics for theta and phi
template <typename T>
std::complex<T> Ylm(const int l, const int m, const T theta, const T phi);

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
	bool tunnel
);

template <typename T>
std::vector<T> zeroCrossing(const std::vector<T>& s, const std::vector<T>& en);

plot_matplotlib plot;
std::string colours(const int& t);

void test_fillhermites();
void Haidingers_brush();
void ODE_Bessels_equation();
void ODE_Bessels_equationQ();

int main(int argc, char* argv[]) {

	//trapezoidal();

	//ODE_test_nl(1);
	//ODE_harmonic_oscillator();

	//ODE_quantum_harmonic_oscillator();

	//ODE_quantum_harmonic_oscillator_complex();

	//testk1();
	//for (int x = 0; x <= 8; x++)
	//ODE_Quantum_Solver(2);

	//ODE_Bessels_equation();
	ODE_Bessels_equationQ();

	//test_fillhermites();
	//Haidingers_brush();

	//ODE_test_poly();
	//tal();

	//ODE_Predator_Prey();
	//quantum_harmonic_oscillator();

	//ODE_Lorenz_System();

	//ODE_Van_der_Pol_oscillator();

	//Leapfrog_integration();

	return 0;
}

void Haidingers_brush()
{
	double E = 3;

	auto mg = nc::meshgrid(nc::linspace(-E, E, 1000), nc::linspace(-E, E, 1000));

	//auto zz = nc::exp(-(nc::power(mg.first, 2) + nc::power(mg.second, 2))) *
		//(nc::power(mg.second, 2) + nc::power(mg.second, 2));
	//zz = np.exp(-(xx * *2 + yy * *2)) * (yy * *2 + yy * *2)
	auto zz = nc::exp(-(nc::power(mg.first, 2) + nc::power(mg.second, 2))) *
		(nc::power(mg.first, 2) - nc::power(mg.second, 2));
	//zz = np.exp(-(xx * *2 + yy * *2)) * (xx * *2 - yy * *2)

	plot.set_title(U"Haidinger's Brush");
	plot.imshow(zz.str(), "viridis", E);

	plot.show();
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
	auto en = linspace(Vo_b, Vo_e, int(8 * abs(Vo_e - Vo_b)));

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
	if (E_zeroes2.empty()) { std::cout << "No roots found !\n\n"; exit(0); }

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

template <typename T>
std::tuple<std::vector<std::vector<T>>, std::vector<T>> ODE_Q_sine_cosine_Rectangular_potential_barrier
(
	const T& Vo_b,
	const T& Vo_e,
	const T& tmin,
	const T& tmax,
	const T& h,
	const T& B1,
	const T& B2,
	const T& Reset
)
{
	std::vector<T> y(4);
	std::vector<T> state(y.size());
	std::vector<std::vector<T>> Y(y.size());

	T x = tmin;
	T E = 0;
	auto en = linspace(Vo_b, Vo_e, int(2 * abs(Vo_e - Vo_b)));

	auto SE = [&](const auto& x, const auto& psi) {

		state[0] = psi[1];
		state[2] = psi[3];

		state[1] = -2 * E * psi[0];
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
			//EI2(x, y);
			//Euler_method(SE, x, y, h);
			//Embedded_Fehlberg_3_4(SE, x, y, h);
			//Embedded_Fehlberg_7_8(SE, x, y, h);

			if (x > B1 && x < B2)
			{
				y *= pow(Reset, 2);
			}

			else if (x > B2) y *= Reset;

			x += h;
			for (size_t i = 0; i < y.size(); i++)
				Y[i].emplace_back(y[i]);
		}

		return *Y[2].rbegin();
	};

	std::vector<T> E_zeroes1, E_zeroes2;

	//E_zeroes1 = Find_all_zeroes(Wave_function, en);

	E_zeroes2 = Find_all_zeroes(Wave_function, en, 0);
	if (E_zeroes2.empty()) { std::cout << "No roots found !\n\n"; exit(0); }

	std::cout << E_zeroes2;

	std::vector<std::vector<T>> psi_sol;

	for (auto& E : E_zeroes1) {
		//Wave_function(E);
		//psi_sol.push_back(Y0);
	}

	for (auto& E : E_zeroes2) {
		Wave_function(E);
		psi_sol.emplace_back(Y[3]);
		//std::reverse(Y[3].begin(), Y[3].end());
		//psi_sol.emplace_back(Y[3]);
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

	double h = 0.005;

	std::vector<double> y(2);

	std::vector<double> state(y.size());
	std::vector<double> X;
	std::vector<std::vector<double>> Y(y.size());

	double sigma = 1, mu = 0;
	double Vo = 20, E = 0;

	if (mode == 2 || mode == 7 || mode == 8) {

		tmin = -13; tmax = 13;

		E = 50; Vo = 60;//50
	}

	bool wave_packet = 0;

	double L1 = -1., L2 = 1.;

	bool tunnel = false;
	if (mode == 3 || mode == 4) {
		tunnel = true;
		wave_packet = true;
		sigma = 32, mu = 1;
		Vo = 0; tmin = -9; tmax = 11;
		h = 0.005;
		E = 0;
	}

	double Vb = 0;
	if (mode == 4) {
		Vo = 0;
		double ws = .1;
		L1 = 2 * mu, L2 = 2 * mu + ws;
		plot.line(L1, L1, 0, 2060);
		plot.line(L2, L2, 0, 2060);
		Vb = 20;
		std::ostringstream out;
		out.precision(2);
		out << std::fixed << "W   = " << ws;
		plot.text(L2 + 1.4, 2000, out.str(), "black", 12);
		plot.arrow(3.4, 2020, 2.1, 2020, "green");
		out.str(std::string());
		out << std::fixed << "Vb  = " << Vb;
		plot.text(L2 + 1.4, 1850, out.str(), "black", 12);
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

		else return -2 * E + (x * x) / sigma;
	};

	double MU = 5;

	double a = -0, b = -1, c = -0, d = -1;//default
	a = -0, b = -1, c = -sqrt(3), d = -sqrt(3);

	auto SE = [&](const auto& x, const auto& psi) {

		if (mode == 8)
			state[0] = a - b * psi[1];
		else
			state[0] = psi[1];

		if (mode == 0 || mode == 1)
			state[1] = 2 * (V(x) - E) * psi[0];

		else if (mode == 7)
			state[1] = V(x - mu) * (psi[0]) - MU * psi[1] * (psi[0] * psi[0] - 1.0);
		// dydx[1] =          -y[0]  - mu *   y[1] *   (y[0] *   y[0] - 1.0);
		else if (mode == 8)
			state[1] = V(x - mu) * (c - d * psi[0]) - MU * psi[1] * (psi[0] * psi[0] - 1.0);


		else state[1] = V(x - mu) * psi[0];

		return state;
	};

	for (x = tmin; x <= tmax; x += h)
		X.emplace_back(x);
	//X = linspace(tmin, tmax, 5201);

	auto Wave_function = [&](const auto& energy) {
		E = energy;

		x = tmin;

		y[0] = 0;
		y[1] = 1;

		Y.clear();
		Y.resize(y.size());

		Y[0].emplace_back(y[0]);
		Y[1].emplace_back(y[1]);

		for (auto& x : X)
		{
			//Midpoint_method_explicit(SE, x, y, h);

			Embedded_Prince_Dormand_v3_4_5(SE, x, y, h);
			//Embedded_Fehlberg_3_4(SE, x, y, h);

			//Midpoint_method_implicit(SE, x, y, h);
			//Euler_method(SE, x, y, h);

			//Embedded_Fehlberg_7_8(SE, x, y, h);

			//Y[0].emplace_back(y[0] * y[0]);
			//Y[1].emplace_back(y[1] * y[1]);

			Y[0].emplace_back(y[0]);
			Y[1].emplace_back(y[1]);

		}
		return y[0];
	};

	std::cout.setf(std::ios::fixed, std::ios::floatfield);
	std::cout.precision(8);//8

	int t = 0;
	std::ostringstream oss;
	oss.setf(std::ios::fixed);
	oss.precision(3);

	std::vector<std::vector<double>> psi_sola, psi_solb, psi_sols;
	std::vector<double> E_zeroes, en;

	if (Vo == 0) {
		en = linspace(0.0, 0.6, 32);
	}
	else en = linspace(E, Vo, int(2 * (Vo - E)));

	E_zeroes = Find_all_zeroes(Wave_function, en, tunnel);
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
		auto v = ODE_Q_sine_cosine(4., 8., tmin, tmax, h);
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

	if (mode != 5 && mode != 6) {
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
	}

	else
	{
		if (mode == 5) {
			tmin = -4;
			tmax = 4;
		}

		else {
			tmin = -8;
			tmax = 8;
			//h = 1e-3;
		}

		X.clear();
		for (x = tmin; x <= tmax; x += h)
			X.emplace_back(x);

		std::tuple<std::vector<std::vector<double>>, std::vector<double>> v;
		if (mode == 5)v = ODE_Q_sine_cosine(0., 4., tmin, tmax, h);

		else {
			double B1 = -.5, B2 = .5;

			v = ODE_Q_sine_cosine_Rectangular_potential_barrier(0., 8., tmin, tmax, h, B1, B2, 0.995);
			plot.line(B1, B1, -1., 1.);
			plot.line(B2, B2, -1., 1.);
		}

		t = 0;
		for (auto& i : std::get<0>(v))
		{
			oss.str(std::string());
			oss << std::get<1>(v)[t];
			plot.plot_somedata(X, i, "k", "E = " + oss.str() + " ", colours(t), 1.0);
			t++;
		}
	}

	std::u32string title;
	if (mode == 0)title = U"Finite potential well";
	else if (mode == 1)title = U"Infinite potential well = Particle in 1-D Box";
	else if (mode == 2)title = U"Quantum Harmonic oscillator";
	else if (mode == 7) title = U"Quantum Harmonic oscillator, van der Pol";
	else if (mode == 8) title = U"Quantum Harmonic oscillator, van der Pol, Predator-Prey";
	//title = U"Quantum Harmonic oscillator, chained 𓂀 = 𓄍 𓄎";
	else if (mode == 3)title = U"Gaussian wave packet";
	else if (mode == 4) title = U"Quantum Gaussian wave packet + tunnelling through a rectangular potential barrier";
	else if (mode == 5) title = U"Quantum Sine and cosine wave";

	else title = U"Quantum Cosine wave + tunnelling through a rectangular potential barrier";

	//plot.set_title(title, "Segoe UI Historic", 20);
	plot.set_title(title);
	plot.grid_on();
	plot.show();

	Wave_function(E_zeroes.front());
	auto v = Y[1];

	if (mode != 5 && mode != 6) {
		t = 0;
		for (auto& E : E_zeroes) {
			Wave_function(E);
			oss.str(std::string());
			oss << E;
			plot.plot_somedata(Y[0], v, "k", "E = " + oss.str() + " ", colours(t++), 1.0);
			v = Y[1];
		}
		//plot.set_title(title, "Segoe UI Historic", 20);
		plot.set_title(title);
		plot.show();

		std::cout.setf(std::ios::fixed, std::ios::floatfield);
		std::cout.precision(15);
		auto p = std::minmax_element(begin(Y[0]), end(Y[0]));
		std::cout << "minY0 = " << *p.first << ", maxY0 = " << *p.second << '\n';
		p = std::minmax_element(begin(Y[1]), end(Y[1]));
		std::cout << "minY1 = " << *p.first << ", maxY1 = " << *p.second << '\n';
	}
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
	plot.set_title(title);

	if (e_plot)plot.show();

	std::cout.setf(std::ios::fixed, std::ios::floatfield);
	std::cout.precision(15);
	auto p = std::minmax_element(begin(Y0), end(Y0));
	std::cout << "minY1 = " << *p.first << ", maxY1 = " << *p.second << '\n';
}


void ODE_Bessels_equationQ()
{
	double tmin = 1e-15;
	double x = tmin;
	double tmax = 10;
	double h = .005;
	double E = 0;

	std::vector<double> y(2);

	y[0] = 1;
	y[1] = 0;

	std::vector<double> dydx(y.size());
	std::vector<double> X, en;
	std::vector<std::vector<double>> Y(y.size());

	for (x = 0; x <= 10; x += 0.1)
		en.emplace_back(x);

	for (x = tmin; x <= tmax; x += h)
		X.emplace_back(x);

	auto SE = [&](const auto& x, const auto& y) {

		const double 	nu = 0.0;

		dydx[0] = y[1];
		dydx[1] = 1.0 / pow(x, 2) * (-x * y[1] - (pow(x, 2) - pow(nu, 2)) * E * y[0]);

		return dydx; };

	auto Wave_function = [&](const auto& energy) {
		E = energy;

		Y.clear();
		Y.resize(y.size());

		x = tmin;

		y[0] = 1;
		y[1] = 0;

		for (size_t i = 0; i < y.size(); i++)
			Y[i].emplace_back(y[i]);

		for (auto& x : X)
		{
			Midpoint_method_explicit(SE, x, y, h);
			//Midpoint_method_implicit(SE, x, y, h);
			//Euler_method(SE, x, y, h);
			//Embedded_Fehlberg_3_4(SE, x, y, h);
			//Embedded_Fehlberg_7_8(SE, x, y, h);

			for (size_t i = 0; i < y.size(); i++)
				Y[i].emplace_back(y[i]);
		}

		return y[0];
	};

	std::vector<double> E_zeroes;

	E_zeroes = Find_all_zeroes(Wave_function, en, 0);
	if (E_zeroes.empty()) { std::cout << "No roots found !" << std::endl << std::endl; exit(0); }

	std::vector<std::vector<double>> psi_sola, psi_solb;

	for (auto& E : E_zeroes) {
		Wave_function(E);
		psi_sola.emplace_back(Y[0]);
		psi_solb.emplace_back(Y[1]);
	}

	int t = 0;
	std::ostringstream oss;
	oss.setf(std::ios::fixed);
	oss.precision(6);

	for (auto& E : E_zeroes) {
		Wave_function(E);
		oss.str(std::string());
		oss << E;
		plot.plot_somedata(X, Y[0], "k", "E = " + oss.str() + " ", colours(t++), 1.0);

		std::cout.setf(std::ios::fixed, std::ios::floatfield);
		std::cout.precision(15);
		std::cout << "E = " << oss.str() << std::endl;
	}

	std::u32string title = U"Quantum Bessel function";
	plot.set_title(title);
	plot.grid_on();

	plot.show();
}

void ODE_Bessels_equation()
{
	double tmin = 1e-15;
	double x = tmin;
	double tmax = 10;
	double h = .005;

	std::vector<double> y(2);

	y[0] = 1;
	y[1] = 0;

	std::vector<double> dydx(y.size());

	auto func = [&](const auto& x, const auto& y) {

		const double 	nu = 0.0;

		dydx[0] = y[1];
		dydx[1] = 1.0 / pow(x, 2) * (-x * y[1] - (pow(x, 2) - pow(nu, 2)) * y[0]);

		return dydx; };

	std::vector<double> X = { tmin }, Y0 = { y[0] }, Y1 = { y[1] };

	while (x <= tmax)
	{
		Embedded_Fehlberg_7_8(func, x, y, h);
		x += h;
		X.push_back(x);
		Y0.push_back(y[0]);
		Y1.push_back(y[1]);
	}

	plot.plot_somedata(X, Y0, "k", "y^1", "blue", 1);
	plot.plot_somedata(X, Y1, "k", "y^2", "red", 1);

	std::u32string title = U"Bessels equation";
	plot.set_title(title);

	plot.show();

	std::cout.setf(std::ios::fixed, std::ios::floatfield);
	std::cout.precision(15);
	auto p = std::minmax_element(begin(Y0), end(Y0));
	std::cout << "minY0 = " << *p.first << ", maxY0 = " << *p.second << '\n';
	p = std::minmax_element(begin(Y1), end(Y1));
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

	plot.plot_somedata(X, Y0, "k", "sine", "blue", 1);
	plot.plot_somedata(X, Y1, "k", "cosine", "red", 1);

	std::u32string title = U"ÿ = -y";
	plot.set_title(title);

	plot.show();

	std::cout.setf(std::ios::fixed, std::ios::floatfield);
	std::cout.precision(15);
	auto p = std::minmax_element(begin(Y0), end(Y0));
	std::cout << "minY0 = " << *p.first << ", maxY0 = " << *p.second << '\n';
	p = std::minmax_element(begin(Y1), end(Y1));
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
		//Midpoint_method_explicit(func, x, y, h);
		//Midpoint_method_implicit(func, x, y, h);
		//Euler_method(func, x, y, h);
		//Embedded_Fehlberg_3_4(func, x, y, h);
		Embedded_Fehlberg_7_8(func, x, y, h);

		x += h;
		X.push_back(x);

		Y0.push_back(y[0] * y[0]);
		Y1.push_back(-y[1]);
	}

	plot.plot_somedata(X, Y0, "k", "Y[0], n = " + std::to_string(n) + "", "red");
	//plot.plot_somedata(X, Y1, "k", "Y[1]", "blue");

	std::u32string title = U"Hermite functions Ψn̈(x) + (2n + 1 - x²) Ψn(x) = 0";

	plot.set_title(title);
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
	plot.set_title(title);
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

	plot.set_title(title);
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
	plot.set_title(title);
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
	plot.set_title(title);
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
	bool tunnel
)
{
	//Gives all zeroes in y = Psi(x)
	const T epsilon = 1e-12;
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
				if (tunnel)break;
			}
		}
	}
	delete brent, dekker, secant;
	return all_zeroes;
}

std::string colours(const int& t)
{
	std::vector<std::string> colours = { "Blue", "Green",
																"Red", "Cyan", "Magenta", "Yellow", "Black", "Silver" };

	auto v = colours;
	for (int i = 0; i < 16; i++)
		colours.insert(colours.end(), v.begin(), v.end());

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

// calculate \tilde{P}_m^m
template <typename T>
T tildePmm(const int m, const T x) {
	T factor = ((m % 2) == 0) ? 1.0 : -1.0;
	const T tmp = 1.0 - x * x;
	T tmp2 = 1.0;
	for (int i = 1; i <= m * 2 - 1; i += 2) {
		tmp2 *= tmp * static_cast<T>(i) / static_cast<T>(T(i) + 1);
	}
	factor *= std::sqrt((2 * T(m) + 1) * tmp2 / (4.0 * pi));
	return factor;
}

template <typename T>
T tildePmm_derivative(const int m, const T x) {
	const T tilde_pmm = tildePmm(m, x);
	const T result = 0.5 * m * (-2.0) * x / (1 - x * x) * tilde_pmm;
	return result;
}

// calculate \tilde{P}_{(m+1)}^m
template <typename T>
double tildePmm1(const int m, const T x) {
	return x * std::sqrt(2.0 * m + 3) * tildePmm(m, x);
}

template <typename T>
T tildePmm1_derivative(const int m, const T x) {
	return x * std::sqrt(2.0 * m + 3) * tildePmm_derivative(m, x) + std::sqrt(2.0 * m + 3) * tildePmm(m, x);
}

// calculate \tilde{P}_m^l
template <typename T>
T tildePlm(const int l, const int m, const T x) {
	if (m < 0 || m > l || std::abs(x) > 1.0) {
		std::cerr << "Bad arguments in tildePlm.\n";
		return 0;
	}
	const T factor1 = std::sqrt((4.0 * l * l - 1) / (T(l) * l - T(m) * m));
	const T factor2 = std::sqrt(((T(l) - 1) * (T(l) - 1) - T(m) * m) / (4.0 * (T(l) - 1) * (T(l) - 1) - 1));
	if (l == m) return tildePmm(m, x);
	if (l == m + 1) return tildePmm1(m, x);
	return factor1 * (x * tildePlm(l - 1, m, x) - factor2 * tildePlm(l - 2, m, x));
}

template <typename T>
T tildePlm_derivative(const int l, const int m, const T x) {
	if (m < 0 || m > l || std::abs(x) > 1.0) {
		std::cerr << "Bad arguments in tildePlm.\n";
		return 0;
	}
	const T factor1 = std::sqrt((4.0 * l * l - 1) / (T(l) * l - T(m) * m));
	const T factor2 = std::sqrt(((T(l) - 1) * (T(l) - 1) - T(m) * m) / (4.0 * (T(l) - 1) * (T(l) - 1) - 1));
	if (l == m) return tildePmm_derivative(m, x);
	if (l == m + 1) return tildePmm1_derivative(m, x);
	return factor1 * (tildePlm(l - 1, m, x) + x * tildePlm_derivative(l - 1, m, x) - factor2 * tildePlm_derivative(l - 2, m, x));
}

// calculate the spherical harmonics for theta and phi
template <typename T>
std::complex<T> Ylm(const int l, const int m, const T theta, const T phi) {
	const std::complex<T> tmp(0.0, static_cast<T>(m) * phi);
	return tildePlm(l, m, std::cos(theta)) * std::exp(tmp);
}

// calculate the spherical harmonics for a vector in Cartesian coordinates
template <typename T>
std::complex<T> Ylm(const int l, const int m, const T x, const T y, const T z) {
	const T norm = x * x + y * y + z * z;
	const T cosine = z / std::sqrt(norm);
	const std::complex<T> eiphi(x / std::sqrt(x * x + y * y), y / std::sqrt(x * x + y * y));
	return tildePlm(l, m, cosine) * std::exp(m) * eiphi;
}

template <typename T>
T Plm(const int l, const int m, const T x) {
	T tmp1 = 1.0;
	T tmp2 = 1.0;
	if (m < 0 || m > l || std::abs(x) > 1.0) {
		std::cerr << "Bad arguments in Plm.\n";
		return 0;
	}
	for (int i = 1; i <= (l - m); ++i) {
		tmp1 *= static_cast<T>(i);
	}
	for (int i = 1; i <= (l + m); ++i) {
		tmp2 *= static_cast<T>(i);
	}
	const T factor = std::sqrt((4.0 * pi * tmp2) / ((2.0 * l + 1) * tmp1));
	return factor * tildePlm(l, m, x);
}

template <typename T>
T Plm_derivative(const int l, const int m, const T x) {
	T tmp1 = 1.0;
	T tmp2 = 1.0;
	if (m < 0 || m > l || std::abs(x) > 1.0) {
		std::cerr << "Bad arguments in Plm.\n";
		return 0;
	}
	for (int i = 1; i <= (l - m); ++i) {
		tmp1 *= static_cast<T>(i);
	}
	for (int i = 1; i <= (l + m); ++i) {
		tmp2 *= static_cast<T>(i);
	}
	const T factor = std::sqrt((4.0 * pi * tmp2) / ((2.0 * l + 1) * tmp1));
	return factor * tildePlm_derivative(l, m, x);
}

template <typename F, typename T>
T numericalDerivative(F& f, T x, T width) {
	return (f(x + width) - f(x - width)) / (width * 2.0);
}

template <typename F1, typename F2, typename T>
void checkDerivative(F1& f, F2& df, T x, const std::string& function_name) {
	static const int magnitude = 5;
	double delta_x = 0.1;
	const T analytical_derivative = df(x);
	std::cout << "Checking derivative of " << function_name << " for x = " << x << '\n';
	std::cout << "Analytical derivative: " << analytical_derivative << std::endl;
	for (int i = 1; i <= magnitude; ++i) {
		const T numerical_derivative = numericalDerivative(f, x, delta_x);
		const T rmse = std::sqrt((numerical_derivative - analytical_derivative) * (numerical_derivative - analytical_derivative));
		std::cout << "Delta_x = " << delta_x << " ; "
			<< "numerical derivative = " << numerical_derivative << " ; "
			<< "RMSE = " << rmse << '\n';
		delta_x = delta_x / 10;
	}
	std::cout << "=========================================\n";
}

//https://github.com/HanatoK/Associated-Legendre-Polynomials/blob/main/associated_legendre_polynomials.cpp
void tal() {
	int l = 6;
	int m = 2;
	double x = 0.5;

	std::cout << std::fixed << std::setprecision(9);
	std::cout << "tildePlm(l, m, x) = " << tildePlm(l, m, x) << std::endl;

	std::cout << "Plm: " << Plm(l, m, x) << std::endl;

	std::cout << "STL (no Condon-Shortley phase term) Plm: " << std::assoc_legendre(l, m, x) << std::endl;

	using namespace std::placeholders;
	auto f1 = std::bind(tildePmm<double>, m, _1);
	auto df1 = std::bind(tildePmm_derivative<double>, m, _1);
	checkDerivative(f1, df1, x, "tildePmm(m, x)");
	auto f2 = std::bind(tildePmm1<double>, m, _1);
	auto df2 = std::bind(tildePmm1_derivative<double>, m, _1);
	checkDerivative(f2, df2, x, "tildePmm1(m, x)");
	auto f3 = std::bind(tildePlm<double>, l, m, _1);
	auto df3 = std::bind(tildePlm_derivative<double>, l, m, _1);
	checkDerivative(f3, df3, x, "tildePlm(l, m, x)");
	auto f4 = std::bind(Plm<double>, l, m, _1);
	auto df4 = std::bind(Plm_derivative<double>, l, m, _1);
	checkDerivative(f4, df4, x, "Plm(l, m, x)");
	std::cout << "Spherical harmonics:\n";
	const double theta = 50.0 / 180.0 * pi;
	const double phi = -30.0 / 180.0 * pi;
	std::cout << Ylm(l, m, theta, phi) << std::endl;

}


#pragma warning( disable : 4129 )
using namespace std;

const int MAXHERMITES = 10;
const int SIZe = 400;
const int X = 4;//4
const int Y = 3;//3
const double W = 80;

int hermites[MAXHERMITES][MAXHERMITES];

double hermite(double x, int n) {
	if (n < 0 || n > MAXHERMITES) {
		cout << "Warning: Hermite polynomial not computed for degree " << n << "." << endl;
		return 0;
	}
	double h = 0;
	for (int i = 0; i <= n; i++) {
		h += pow(x, i) * hermites[n][i];
	}
	return h;
}

double I(double x, int n) {
	return pow(hermite(sqrt(2) * x / W, n) * exp(-x * x / W / W), 2);
}

void svg() {
	char filename[200];
	sprintf(filename, "hermites.svg");
	fstream fout(filename, fstream::out);
	fout << "<svg xmlns='http://www.w3.org/2000/svg' version='1.1' width='" << X * SIZe << "px' height='" << Y * SIZe << "px'>" << endl;
	fout << "<defs>" << endl;
	for (int i = 0; i < X; i++) {
		fout << "<linearGradient id='grad" << i << "' x1='0%' y1='0%' x2='100%' y2='0%'>" << endl;
		double intensities[SIZe], maxintensity = 0;
		for (int j = 0; j < SIZe; j++) {
			intensities[j] = I(j - SIZe / 2, i);
			if (intensities[j] > maxintensity) maxintensity = intensities[j];
		}
		std::vector<double> Y, X;

		for (int j = 0; j < SIZe; j++) {
			fout << fixed << setprecision(4) << "<stop offset='" << j * 100.0L / SIZe << "\%' style='stop-color:#000; stop-opacity:" << 1 - intensities[j] / maxintensity << ";' />" << endl;

			Y.push_back(1 - intensities[j] / maxintensity);
			X.push_back(j * 100.0L / SIZe);
		}

		std::u32string title;
		title = U"Test hermite 2d";
		std::wstring_convert<std::codecvt_utf8<char32_t>, char32_t> cv1;
		plot.plot_somedata(X, Y, "o", cv1.to_bytes(title), "Red", 1.0, 1);
		plot.show();

		fout << "</linearGradient>" << endl;
	}
	fout << "</defs>" << endl;
	//fout << "<rect x='0' y='0' height='1024' width='1024' fill='url(#grad4)' />" << endl;
	for (int x = 0; x < X; x++) {
		for (int y = 0; y < Y; y++) {
			fout << "<rect x='0' y='0' width='" << SIZe << "' height='" << SIZe << "' fill='url(#grad" << x << ")' transform='matrix(1,0,0,1," << x * SIZe << "," << y * SIZe << ")'/>" << endl;
			fout << "<rect x='0' y='0' width='" << SIZe << "' height='" << SIZe << "' fill='url(#grad" << y << ")' transform='matrix(0,1,1,0," << x * SIZe << "," << y * SIZe << ")'/>" << endl;
		}
	}
	fout << "</svg>" << endl;
}

void fillhermites() {
	for (int i = 1; i < MAXHERMITES; i++) {
		hermites[0][i] = 0;
	}
	hermites[0][0] = 1;
	for (int j = 1; j < MAXHERMITES; j++) {
		hermites[j][0] = 0;
		for (int i = 1; i <= MAXHERMITES; i++) { // up to j would suffice, but it's better to set higher coefficients to zero
			hermites[j][i] = 2 * hermites[j - 1][i - 1];
			hermites[j][i - 1] -= i * hermites[j - 1][i];
		}
	}
	return;
}


void test_fillhermites() {
	fillhermites();
	svg();
}


template <typename T>
inline T AssociatedLegendre
(
	const int  m_l = 3,
	const int  m_m = 8,
	const T& x = 0
)
{
	T pmm = 1.0;
	T pmp1m;

	// Check Inputs
	//static_assert (fabs(x) > 1.0, "Input Out Of Range");

	if (m_m > 0)
	{
		// P_m^m(x) = (-1) ^ m(2m - 1)!!(1 - x ^ 2) ^ m / 2
		T sqrtomx2 = sqrt(1.0 - x * x);
		T oddInt = 1.0;
		for (int i = 1; i <= m_m; ++i)
		{
			pmm *= -1.0 * oddInt * sqrtomx2;
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
		pmp1m = x * (2.0 * T(m_m) + 1.0) * pmm;
		if (m_l == m_m + 1)
		{
			return pmp1m;
		}
		else
		{

			T pmp2m = 0.0;
			for (int i = m_m + 2; i <= m_l; ++i)
			{

				// (l - m) P_l^m (x) = x (2l - 1) P_l-1^m (x) - (l + m - 1) P_l-2^m (x)
				pmp2m = x * (2.0 * T(m_l) - 1.0) * pmp1m - (T(m_l) + T(m_m) - 1) * pmm;
				pmp2m /= T(m_l) - T(m_m);
				// Rotate variables for next iteration
				pmm = pmp1m;
				pmp1m = pmp2m;
			}
			//std::cout << x << " " << pmp2m << std::endl;
			return pmp2m;
		}
	}
}

void ODE_test_poly()
{
	double tmin = -2;
	double tmax = 2;
	double h = 0.00008;

	std::vector<double> y(2);
	std::vector<double> state(y.size());

	double k = 100, L = 1, E = 49, m = 100, H = 1;

	auto en = linspace(E, E + 1, 4);

	auto V = [&](const auto& x) {

		if (abs(x) < L)
			return  0.5 * k * x * x;

		else return  0.5 * k * (L * L);
	};

	auto SE = [&](const auto& x, const auto& psi) {

		state[0] = psi[1];

		state[1] = (2.0 * m) * (V(x) - E) * psi[0];

		return state;
	};

	//std::vector<double> X = { x }, Y0 = { tildePlm(m, l, x) }, Y1 = { tildePlm(m, l, x) };
	std::vector<double> Y0, Y1, X;

	for (auto x = tmin; x <= tmax; x += h)
		X.emplace_back(x);

	//X = linspace(tmin, tmax, 50001);
	//auto Xnorm = linspace(-1., 1., 1000);

	auto Wave_function = [&](const auto& energy) {
		E = energy;

		Y0.clear();
		Y1.clear();

		y[0] = .001;
		y[1] = 0;

		//while (x <= tmax)
		for (auto& x : X)
		{
			Midpoint_method_explicit(SE, x, y, h);
			//Midpoint_method_implicit(SE, x, y, h);
			//Euler_method(fSE, x, y, h);
			//Embedded_Fehlberg_3_4(SE, x, y, h);
			//Embedded_Fehlberg_7_8(SE, x, y, h);

			//Y0.push_back(pow(std::hermite(m, x),2));

			//Y0.push_back(pow(ps::Hermite_function(m, std::cosh(-pi*x)), 1)+ pow(ps::Hermite_function(l, std::cosh(pi*x)), 1));
			//Y0.push_back(pow(std::assoc_legendre(m, l, std::cos(x*pi)), 2));
			//Y1.push_back(pow(AssociatedLegendre(m, l, std::cos(x)), 1));
			//Y0.push_back(pow(std::legendre(m-2, x), 1));
			//Y0.push_back(pow(Ylm(m, l,x, 0.).real(),1));


			Y0.push_back(pow((y[0]), 2));
			Y1.push_back(pow((y[1]), 2));

		}

		return y[0];
	};


	bool mp = 1;
	if (mp) {
		std::vector<double> E_zeroes;
		E_zeroes = Find_all_zeroes(Wave_function, en, 0);
		if (E_zeroes.empty()) { std::cout << "No roots found !\n\n"; exit(0); }

		std::stringstream oss;
		int t = 0;
		for (auto& E : E_zeroes) {
			Wave_function(E);
			oss.str(std::string());
			oss << E;
			std::cout << "E " << E << std::endl;
			//Y0 *= pow(t, 100);
			plot.plot_somedata(X, Y0, "k", "E = " + oss.str() + " ", colours(t), 1.0);
			//plot.plot_somedata(X, Y1, "k", "E = " + oss.str() + " ", colours(t+1), 1.0);
			//if (t == 2)break;
			t++;
		}
	}
	else {
		Wave_function(1);
		//plot.plot_somedata(X, Y0, "k", "Y[0], m = " + std::to_string(m) + ", l = " + std::to_string(l) + "", "red", 1);
		plot.plot_somedata(X, Y1, "k", "Y[1]", "blue", 1);
	}
	//https://root.cern.ch/doc/v610/LegendreAssoc_8C.html
	std::u32string title = U"Wave function";
	plot.set_title(title);
	plot.grid_on();
	plot.show();

	//plot.plot_somedata(Y0, Y1, "k", "Y[0] vs Y[1]", "green", 1);
	//plot.show();

	std::cout.setf(std::ios::fixed, std::ios::floatfield);
	std::cout.precision(15);
	auto p = std::minmax_element(begin(Y0), end(Y0));
	std::cout << "minY0 = " << *p.first << ", maxY0 = " << *p.second << '\n';
	p = std::minmax_element(begin(Y1), end(Y1));
	std::cout << "minY1 = " << *p.first << ", maxY1 = " << *p.second << '\n';

}


