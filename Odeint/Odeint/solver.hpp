#pragma once
#include <matplotlib.hpp>

static plot_matplotlib plot;

template <typename F, typename A, typename B>
inline A MC1(F func, const A& a, const B& b, const int& M)
{
	A sum = 0;
	mxws rng;

	for (int i = 1; i <= M; i++)
	{
		sum += func(rng(a, b));
	}
	return ((b - a) / A(M)) * sum;
}

namespace np
{
	template<size_t n_samples, typename T>
	std::array<T, n_samples> random_uniform(const std::array<T, 2>& domain)
	{
		mxws mxws;

		std::array<T, n_samples> arr;

		auto low = domain[0];
		auto high = domain[1];

		for (auto& i : arr)
			i = mxws(low, high);

		return arr;
	}

	template<size_t n_samples, typename F, typename T>
	T mean1(F& func, const std::array<T, n_samples>& samples1)
	{
		T sum = 0;

		for (size_t i = 0; i < n_samples; i++)
			sum += func(samples1[i]);

		return sum / T(n_samples);
	}

	template<size_t n_samples, typename F, typename T>
	T ode(F& func, const std::array<T, 2>& domain)
	{
		mxws mxws;
		auto low = domain[0];
		auto high = domain[1];
		T sum = 0;

		for (size_t i = 0; i < n_samples; i++)
			sum += func(mxws(low, high));

		return sum / T(n_samples);
	}

	template<size_t n_samples, typename F, typename T>
	T ode2b(F& func, const std::array<T, 2>& domain)
	{
		mxws mxws;
		auto low = domain[0];
		auto high = domain[1];
		T sum = 0;

		for (size_t i = 0; i < n_samples; i++) {
			
			sum += func(mxws(low, high), mxws(low, high));
		}
		auto volume = abs(domain[1] - domain[0]);

		return (sum / T(n_samples)) * volume;
	}
}

template <size_t n_samples, typename F, typename T>
T mc_int1(F& func, const std::array<T, 2>& domain) {

	auto samples1 = np::random_uniform<n_samples>(domain);
	auto volume = abs(domain[1] - domain[0]);
	return np::mean1(func, samples1) * volume;
}

template <size_t n_samples, typename F, typename T>
T mc_int2(F& func, const std::array<T, 2>& domain) {

	auto volume = abs(domain[1] - domain[0]);
	return np::ode<n_samples>(func, domain) * volume;
}

template <typename F, typename T, size_t n_samples>
std::array<T, n_samples> mc_solver1(F& func, const T& y0, const std::array<T, n_samples>& x)
{
	//: param func : Lambda function of the gradient field.Signature should be `func(x)`
	//: param y0 : Initial function value corresponding to the x0.The latter is given as the first element of x
	//: param x : Domain points over which the field `func` should be integrated
	//: param n_samples : Number of random samples for the estimate
	std::array<T, 2> domain;

	std::array<T, n_samples> vals = {};
	vals[0] = y0;

	for (size_t i = 0; i < n_samples-1; i++) {
		domain[0] = x[i];
		domain[1] = x[i+1];

		vals[i+1] = vals[i] + mc_int1<n_samples>(func, domain);
	}

	return vals;
}

template <typename F, typename T, size_t n_samples>
std::array<T, n_samples> mc_solver1b(F& func, const T& y0, const std::array<T, n_samples>& x)
{
	//: param func : Lambda function of the gradient field.Signature should be `func(x)`
	//: param y0 : Initial function value corresponding to the x0.The latter is given as the first element of x
	//: param x : Domain points over which the field `func` should be integrated
	//: param n_samples : Number of random samples for the estimate
	std::array<T, 2> domain;

	std::array<T, n_samples> vals = {};
	vals[0] = y0;

	for (size_t i = 0; i < n_samples-1; i++) {
		domain[0] = x[i];
		domain[1] = x[i+1];

		vals[i+1] = vals[i] + mc_int2<n_samples>(func, domain);
	}

	return vals;
}

template <typename F, typename T, size_t n_samples>
std::array<T, n_samples> mc_solver2b
(
	F& func,
	const T& y0,
	const std::array<T, n_samples>& x
)
{
	//: param func: Lambda function of the nonlinear gradient field. Signature should be `func(y0, x)`
	//: param y0 : Initial function value corresponding to the x0.The latter is given as the first element of x
	//: param x : Domain points over which the field `func` should be integrated
	//: param n_samples : Number of random samples for the estimate
	std::array<T, 2> domain;

	std::array<T, n_samples> vals = {};
	vals[0] = y0;

	for (size_t i = 0; i < n_samples - 1; i++) {
		domain[0] = x[i];
		domain[1] = x[i + 1];

		auto part_func = [&](const auto& x, const auto& y) { return func(x, vals[i]); };
		vals[i + 1] = vals[i] + np::ode2b<1000>(part_func, domain);
	}

	return vals;
}

template <typename F, typename elem, int order, size_t n_samples>
std::array<multicomplex<elem, order>, n_samples> mc_solver2c
(
	F& func,
	const elem& y0,
	const std::array<multicomplex<elem, order>, n_samples>& x
)
{
	std::array<multicomplex<elem, order>, 2> domain;

	std::array<multicomplex<elem, order>, n_samples> vals = {};
	vals[0] = y0;

	for (size_t i = 0; i < n_samples - 1; i++) {
		domain[0] = x[i];
		domain[1] = x[i + 1];

		auto part_func = [&](const auto& x, const auto& y) { return func(x, vals[i]); };
		vals[i + 1] = vals[i] + Generalized_midpoint(part_func, domain[0], domain[1], 12, 4);
	}

	return vals;
}

template <typename F, typename elem, size_t n_samples>
std::array<elem, n_samples> mc_solver2d
(
	F& func,
	const elem& y0,
	std::array<elem, n_samples>& t,
	const int& n,
	const std::array<elem, 2>& tspan
)
{
	std::array<elem, n_samples> vals = {};
	
	vals[0] = y0;//initial
	
	t[0] = tspan[0];

	elem dt = (tspan[1] - tspan[0]) / elem(n);

	mxws mxws;

	const size_t totali = 100000;

	for (size_t j = 0; j < n; j++) 
	{
		
		elem sum = 0;
		
		for (int i = 0; i < totali; i++) {

			sum += func(mxws(t[j], t[j] + dt), vals[j]);
		}
		
		sum /= elem(totali);
		
		t[j + 1] = t[j] + dt;

		vals[j + 1] = vals[j] + dt * sum;
	}

	return vals;
}

template <typename F, typename elem, size_t n_samples>
std::array<elem, n_samples> mc_solver2e
(
	F& func,
	const elem& y0,
	std::array<elem, n_samples>& t,
	const int& n,
	const std::array<elem, 2>& tspan
)
{
	std::array<elem, n_samples> vals = {};

	vals[0] = y0;//initial

	t[0] = tspan[0];

	elem dt = (tspan[1] - tspan[0]) / elem(n);

	for (size_t j = 0; j < n; j++)
	{
		elem sum = 0;
		
		sum = Generalized_midpoint(func, t[j], t[j] + dt, vals[j], vals[j], 4, 4); 

		t[j + 1] = t[j] + dt;

		vals[j + 1] = vals[j] + sum;
	}

	return vals;
}



/* The type of container used to hold the state vector */
typedef std::vector< long double > state_type;

struct push_back_state_and_time
{
	std::vector< state_type >& m_states;
	std::vector<double >& m_times;

	push_back_state_and_time(std::vector< state_type >& states, std::vector< double >& times)
		: m_states(states), m_times(times) { }

	void operator()(const state_type& x, double t)
	{
		m_states.push_back(x);
		m_times.push_back(t);
	}
};


void test_boost()
{
	using namespace boost::numeric::odeint;

	//[ state_initialization
	state_type x(2);
	//x[0] = 1.0; // start at x=1.0, p=0.0
	//x[1] = 0.0;

	x[0] = 4.0; // start at x=1.0, p=0.0
	x[1] = 0.0;

	//[ integrate_observ
	std::vector<state_type> x_vec;

	//runge_kutta_fehlberg78< state_type > stepper;
	typedef runge_kutta_fehlberg78< state_type > error_stepper_type;
	typedef controlled_runge_kutta< error_stepper_type > controlled_stepper_type;
	//controlled_stepper_type controlled_stepper;

	std::vector<double> times;
	size_t steps;

	//steps = integrate_const(stepper, [](const state_type& x, state_type& dxdt, double t) {
		//dxdt[0] = x[1]; dxdt[1] = -x[0] - gam * x[1]; }
	//, x, 0.0, 10., 1., push_back_state_and_time(x_vec, times));

	steps = integrate_const(make_controlled< error_stepper_type >(1.0e-12, 1.0e-12), [](const state_type& x, state_type& dxdt, double t) {

		dxdt[0] = t * sqrt(abs(x[0])) + pow(sin(t * pi / 2), 3) - 5. * (t > 2.);
		}
	, x, -4.0, 5.0, .01, push_back_state_and_time(x_vec, times));

	/* output */
	std::vector<double> Y, X;
	for (size_t i = 0; i <= steps; i++)
	{
		//std::cout << times[i] << '\t' << x_vec[i][0] << '\t' << x_vec[i][1] << '\n';
		//std::cout << times[i] << '\t' << x_vec[i][0] << '\t' << '\n';
	}

	for (size_t i = 0; i < x_vec.size(); i++)
	{
		Y.push_back(x_vec[i][0]);
		X.push_back(times[i]);
	}

	// https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta%E2%80%93Fehlberg_method
	std::cout << "runge_kutta_fehlberg78 steps " << steps << std::endl;
	plot.plot_somedata(X, Y, "k", "runge_kutta_fehlberg78 ", "brown");
	plot.show();

	x[0] = 1.0; // start at x=1.0, p=0.0
	x[1] = 0.0;

	x_vec.clear();
	times.clear();
	steps = integrate_const(make_controlled< error_stepper_type >(1.0e-15, 1.0e-15),
		[](const auto& x, auto& dxdt, auto t) {

			//dxdt[0] = t * sqrt(abs(x[0])) + pow(sin(t * pi / 2), 3) - 5. * (t > 2.);

			long double w0 = 2 * pi * 0.25;
			long double zeta = 0;

			auto X = x[0];
			auto p = x[1];
			auto dx = p;
			auto dp = -2. * zeta * w0 * p - pow(w0, 2) * X;
			dxdt[0] = dx;
			dxdt[1] = dp;


		}
	, x, 0., 10., .1, push_back_state_and_time(x_vec, times));

	Y.clear();
	X.clear(); 
	
	for (size_t i = 0; i < x_vec.size(); i++)
	{
		Y.push_back(x_vec[i][1]);
		X.push_back(times[i]);
	}

	std::cout.precision(15);
	std::cout << "runge_kutta_fehlberg78 steps " << steps << std::endl;
	std::cout << Y << std::endl;
	std::cout << "cos(2.5*pi) " << cos(2.5*pi) << std::endl;
	
	plot.plot_somedata(X, Y, "k", "sin ", "red");
	Y.clear();
	X.clear();
	for (size_t i = 0; i < x_vec.size(); i++)
	{
		Y.push_back(x_vec[i][0]);
		X.push_back(times[i]);
	}
	plot.plot_somedata(X, Y, "k", "cos ", "blue");
	plot.show();

}

void ode_python()
{
	plot.set_title("Ordinary differential equation");

	plot.PyRun_Simple("from typing import Tuple, Callable");
	plot.PyRun_Simple("import numpy as np");
	plot.PyRun_Simple("import matplotlib.colors as mcolors");

	plot.PyRun_Simple("from scipy.integrate import quadrature, odeint, solve_ivp");

	plot.PyRun_Simple("def func(y, x):\n\
			return  x * np.sqrt(np.abs(y)) + np.sin(x * np.pi/2)**3 - 5 * (x > 2)");

	plot.PyRun_Simple("def func2(x, y):\n\
			return  x * np.sqrt(np.abs(y)) + np.sin(x * np.pi / 2)**3 - 5 * (x > 2)");

	plot.PyRun_Simple("def mc_int(func: Callable, domain : Tuple, n_samples : int) :\n\
			samples = np.random.uniform(low = domain[0], high = domain[1], size = (n_samples, ))\n\
			#print(samples)\n\
			volume = abs(domain[1] - domain[0])\n\
			return np.mean(func(samples)) * volume");

	plot.PyRun_Simple("def mc_ode_solve(func, y0, t, n_samples = 100000) :\n\
		sols = [y0]\n\
		for lo, hi in zip(t[:-1], t[1:]) :\n\
				part_func = lambda v : func(x = v, y = sols[-1])\n\
				sols.append(sols[-1] + mc_int(part_func, (lo, hi), n_samples = n_samples))\n\
		return np.asarray(sols)");

	plot.PyRun_Simple("base2 = np.linspace(-4, 5, 501)");
	plot.PyRun_Simple("y0 = 4.");
	plot.PyRun_Simple("ys_mc = []");

	plot.PyRun_Simple("ys_mc.append(mc_ode_solve(func, y0, base2))");

	plot.PyRun_Simple("y_ode3 = solve_ivp(func2, (-4, 5), [y0], method = 'RK23', t_eval = base2, dense_output = True)");
	plot.PyRun_Simple("z1 = y_ode3.sol(base2)");
	plot.PyRun_Simple("plt.plot(base2, z1.T, 'k', label= 'solve_ivp method RK23', linewidth=1, color = 'orange')");
	plot.PyRun_Simple("plt.plot(base2, np.asarray(ys_mc).flatten(), 'k', label = 'Python Monte Carlo ode solve', linewidth = 1, color = 'xkcd:sky blue')");
	plot.PyRun_Simple("base = np.linspace(-5, 5, 30)");
	plot.PyRun_Simple("xx, yy = np.meshgrid(base, base)");

	plot.PyRun_Simple("U = 1");
	plot.PyRun_Simple("V = func(yy, xx)");
	plot.PyRun_Simple("N = np.sqrt(U **2 + V **2)");
	plot.PyRun_Simple("U2, V2 = U / N, V / N");
	plot.PyRun_Simple("plt.quiver(xx, yy, U2, V2, color = 'lightgray')");

}

