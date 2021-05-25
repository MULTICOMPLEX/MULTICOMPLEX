
#pragma warning( disable : 6262 6246 6326 6244) 

#include "MULTICOMPLEX.hpp"
#include <rkf78.hpp>

#include "vector_calculus.hpp"
#include <boost/numeric/odeint.hpp>
#include <solver.hpp>

using namespace std;
using namespace boost::numeric::odeint;
void test_boost();

template <typename T>
void rk4(void dyDT(T x, T y[], T yp[]), T tspan[2],
	T y0[], const int n, const size_t m, T t[], T y[]);

template <typename T>
void EULER(void dyDT(T x, T y[], T yp[]), T tspan[2],
	T y0[], const int n, const size_t m, T t[], T y[]);

template <typename T>
void EULER_mc(void dyDT(T x, T y[], T yp[]), T tspan[2],
	T y0[], const int n, const size_t m, T t[], T y[]);

template <typename T>
void midpoint_explicit_mc(void dyDT(T x, T y[], T yp[]),
	T tspan[2], T y0[], const int n, const size_t m, T t[], T y[]);

template <typename T>
void midpoint_explicit2(void dyDT(T x, T y[], T yp[]),
	T tspan[2], T y0[], const int n, const size_t m, T t[], T y[]);

template <typename T>
void midpoint_implicit(void dyDT(T x, T y[], T yp[]),
	T tspan[2], T y0[], const int n, const size_t m, T t[], T y[]);

template <typename T>
void rkf(void dyDT(T x, T y[], T yp[]),
	T tspan[2], T y0[], const int n, const size_t m, T t[], T y[]);

int central_config_test();


template <typename T>
void dydt(T t, T u[], T f[]);
template <typename T>
void dydt2(T t, T u[], T f[]);
template <typename T>
void dydt3(T t, T u[], T f[]);
//template <typename T>
//T* dydt3(T x, T y[]);


template<typename F>
size_t ode2d(F& func)
{
	return func(1, 1);
}


template <typename F, typename T, size_t n_samples>
std::array<T, n_samples> mc_solver2d
(
	F& func,
	const std::array<T, n_samples>& x
)
{
	std::array<T, n_samples> vals = {};

	for (size_t i = 0; i < n_samples - 1; i++) {

		auto part_func = [&](const auto& x, const auto& y) { return func(x, vals[i]); };
		vals[i + 1] = vals[i] + ode2d(func);
	}

	return vals;
}

//0, 1, 3, 7, 15, 31, 63, 127, 255, 511, 1023
//0, 2, 4, 6,  8, 10, 12,  14,  16,  18,   20

int main()
{
	
	return central_config_test();
	

	std::cout.setf(std::ios::fixed, std::ios::floatfield);
	std::cout.precision(15);
	{
		auto mfunc = [](const auto& z, const auto& z2) { return sin(z); };
		REAL t = -12.7;
		MX0 a = -12.7;
		MX0 b = -8.5;
		//(integrate(a),a = 1to10) + (integrate(b),b=1to10)
		//std::cout << "Generalized_midpoint" << Generalized_midpoint(mfunc, a, b, 10, 4) << std::endl;//min 3,4 sin(z)
	//1.593096774499077
	//1.593096774565459
	}
	std::cout.precision(12);
	//return 0;

	ode_python();

	MX0 x,y;

	auto f = [](const auto& x) { return sqrt(x) - pow(x,2); };

	x.random(-10, 10);
	y.random(-10, 10);
	
	//std::cout << Critical_point(f, x, 25);
	Critical_point(f, x, 30);
	
	std::cout << std::endl << std::endl;
	std::cout << "max = " << f(x) << std::endl;
	
	//root(f, x, 50);

	auto f2 = [](const auto& x, const auto& y) { return sqrt(x) - pow(x, 2); };
	//std::vector<MX0> grad = VX::gradient (f2,x,f2,y);//∇f 

	//VX::Hessian2x2(f2, /*x*/x, /*y*/y);
	auto m = VX::inverse_Hessian2x2(f2, x, y);

	x.random(-10, 10);
	y.random(-10, 10);

	std::cout << m;
	std::cout << std::endl << std::endl;
	//auto f3 = [](const auto& x, const auto& y) { return  ((x + y) * (x * y + x * pow(y,2)));};
	auto f4 = [](const auto& x, const auto& y) { return 2 *y + x * (x + y) * (x * cos(y) + x * pow(y, 1.9)); };
	//critical points(2 y + x(x + y) (x y + x y^1.9))

	auto f3 = [](const auto& x, const auto& y) { return pow(x - 1, 2) + 10 * pow(y - pow(x, 2), 2); }; //Rosenbrock function
	//(x - 1)** 2 + b * (y - x * *2) * *2;

	//VX::Critical_point_2(f3, x, y, 30);
	std::cout << std::endl;

	//x.random(-1, 1);
	//y.random(-1, 1);

	x.real =  1;
	x.imag = -1;
	y.real = -0.5;
	y.imag = -0.2;

	std::cout << " x  " << x << "y" << y << std::endl;
	VX::Critical_point_2(f3, x, y, 30);
	
	std::cout << std::endl;
	VX::Critical_point_2(f4, x, y, 30);

	
	auto mc2 = [](const auto& x) { if (x == 0) return 1.; return pow(sin(pi*x) / (pi*x), 2.); };

	int M = 1000000;
	REAL a = -3.5;
	REAL b =  3.5;


	std::cout << std::endl;
	std::cout << "MC1 " << MC1(mc2, a, b, M) << std::endl << std::endl;

	auto xs = linspace<REAL, 500>(-2.5, 2.05);
	REAL y0 = -11.2;
	std::array<REAL, xs.size()> ys = {};
	ys[0] = y0;

	auto func = [](const auto& x) { return pow(x,3) + 2 * pow(x,2) - pow(3,x); };
	//{y'(x) = x^3 + 2  x^2 - 3^x, y(0)=-11.2} from -2.5 to 2.05 by implicit midpoint
	//y(x) = -0.910239 (-0.274653 x^4 - 0.732408 x^3 + 3^x + 11.3045)
	//auto func_exact = [](const auto& x) { return 
		//-0.910239 * (-0.274653 * pow(x,4) - 0.732408 * pow(x,3) + pow(3,x) + 11.3045); };

	
	//std::cout << mc_solver1b(func, y0, xs);
	std::cout << std::endl;

	auto func2 = [](const auto& x, const auto& y) { return x * sqrt(abs(y)) + pow(sin(x * pi / 2), 3) - 5. * (x > 2.); };

	{
		auto func = [](const auto& x, const auto& y) { return x + y ; };
		auto xs = linspace<size_t, 11>(0, 10);
		std::cout << mc_solver2d(func,xs) << std::endl;
	}
	
	//return 0;

	y0 = 4;

	std::vector<double> X;
	std::vector<double> Y;
	
	{
		std::cout << "mc_solver2b(func2, y0, xs2)" << std::endl;

		auto xs2 = linspace<REAL, 501>(-4, 5);
		auto p = mc_solver2b(func2, y0, xs2);

		X.clear();
		Y.clear();
		for (int i = 0; i < p.size(); i++) {
			X.push_back(xs2[i]);
			Y.push_back(p[i]);
		}

		plot.plot_somedata(X, Y, "k", "mc_solver2b", "red");

		//std::cout << Y;
	}
	
	{
		
		std::cout << "mc_solver2c(func2, y0, xs2)" << std::endl;

		auto xs2 = linspace<REAL, 0, 501>(-4, 5);
		std::array<REAL, xs2.size()> xs2r;

		auto p = mc_solver2c(func2, y0, xs2);

		for (int i = 0; i < p.size(); i++)
			xs2r[i] = p[i].real;

		X.clear();
		Y.clear();
		for (int i = 0; i < p.size(); i++) {
			X.push_back(xs2[i].real);
			Y.push_back(p[i].real);
		}

		plot.plot_somedata(X, Y, "k", "Generalized_midpoint", "xkcd:deep lilac");

		//std::cout << Y;
	}

	{
		std::cout << "mc_solver2d(func2, y0, xs2)" << std::endl;

		const int n = 1500; //int n: the number of steps to take. 
		std::array<REAL, n+1> t = {};

		std::array<REAL, 2> tspan;
		tspan[0] = -4;
		tspan[1] = 5;
		
		auto p = mc_solver2d(func2, y0, t, n, tspan);

		X.clear();
		Y.clear();
		for (int i = 0; i < p.size(); i++) {
			X.push_back(t[i]);
			Y.push_back(p[i]);
		}

		plot.plot_somedata(X, Y, "k", "mc_solver2d", "yellow");

		//std::cout << Y;
	}

	{
		std::cout << "mc_solver2e(func2, y0, xs2)" << std::endl;

		const int n = 500; //int n: the number of steps to take. 
		std::array<REAL, n + 1> t = {};

		std::array<REAL, 2> tspan;
		tspan[0] = -4;
		tspan[1] = 5;

		//auto p = mc_solver2e(func2, y0, t, n, tspan);

		X.clear();
		Y.clear();
		//for (int i = 0; i < p.size(); i++) {
			//X.push_back(t[i]);
			//Y.push_back(p[i]);
		//}

		//plot.plot_somedata(X, Y, "k", "mc_solver2e", "black");

		//std::cout << Y;
	}

	{
		REAL tspan[2]; // double tspan[2]: contains the initial and final times.

		tspan[0] = -4;
		tspan[1] = 5;

		const int n = 1500; //int n: the number of steps to take. 
		const size_t m = 1; //int m: the number of variables.

		REAL y0[m]; // double y0[m]: a column vector containing the initial condition.
		y0[0] = 4;
		//y0[1] = 0; 

		// double t[n + 1], y[m * (n + 1)]: the times and solution values.
		REAL t[n + 1];
		REAL y[m * (n + 1)];

		int arrSize = sizeof(y) / sizeof(y[0]);

		std::cout << "rk4" << std::endl;
		rk4(dydt3, tspan, y0, n, m, t, y);

		X.clear();
		Y.clear();
		for (int p = 0, i = 0; i < arrSize; i += 1, p++) {
			X.push_back(t[p]);
			Y.push_back(y[i]);
			//std::cout << " " << t[p] << " " << y[i] << std::endl;
		}
		//plot.plot_somedata(X, Y, "k", "rk4", "red");


		std::cout << "midpoint_implicit" << std::endl;
		midpoint_implicit(dydt3, tspan, y0, n, m, t, y);

		X.clear();
		Y.clear();
		for (int p = 0, i = 0; i < arrSize; i += 1, p++) {
			X.push_back(t[p]);
			Y.push_back(y[i]);
			//std::cout << " " << t[p] << " " << y[i] << std::endl;
		}
		plot.plot_somedata(X, Y, "k", "midpoint_implicit", "blue");


		std::cout << "rkf" << std::endl;
		rkf(dydt3, tspan, y0, n, m, t, y);

		X.clear();
		Y.clear();
		for (int p = 0, i = 0; i < arrSize; i += 1, p++) {
			X.push_back(t[p]);
			Y.push_back(y[i]);
			//std::cout << " " << t[p] << " " << y[i] << std::endl;
		}
		plot.plot_somedata(X, Y, "k", "rkf", "xkcd:kiwi");

	}
	std::cout << std::endl << std::endl;

	test_boost();

}

template <typename T>
void dydt(T t, T u[], T f[])
{
	f[0] = u[1];
	f[1] = -u[0] - 0.15 * u[1]; 
}

double w0 = 2 * pi * 0.25;
double zeta = 0;
//tspan[0] = 0;
//tspan[1] = 10;
//y0[0] = 1; 
//y0[1] = 0; 

template <typename T>
void dydt2(T t, T u[], T f[])
{
	auto x = u[0];
	auto p = u[1];
	auto dx = p;
	auto dp = -2. * zeta * w0 * p - pow(w0, 2) * x;
	f[0] = dx;
	f[1] = dp;
}

//y0[0] = 4; 
//tspan[0] = -4;
//tspan[1] = 5;
//n = 500; //int n: the number of steps to take. 
//m = 1; //int m: the number of variables.
template <typename T>
void dydt3(T t, T u[], T f[])
{
	auto x = t;
	f[0] = x * sqrt(abs(u[0])) + pow(sin(x * pi / 2),3) - 5. * (x > 2.);
}

void dydt3(int t)
{
	
}

//https://people.sc.fsu.edu/~jburkardt/cpp_src/rk4/rk4.cpp
//****************************************************************************80
template <typename T>
void rk4(void dyDT(T x, T y[], T yp[]), T tspan[2],
	T y0[], const int n, const size_t m, T t[], T y[])

	/******************************************************************************/
	//
	//  Purpose:
	//
	//    euler approximates the solution to an ODE using Euler's method.
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    25 April 2020
	//
	//  Author:
	//
	//    John Burkardt
	//
	//  Input:
	//
	//    dydt: a function that evaluates the right hand side of the ODE.
	//
	//    T tspan[2]: contains the initial and final times.
	//
	//    T y0[m]: a column vector containing the initial condition.
	//
	//    int n: the number of steps to take.
	//
	//    int m: the number of variables.
	//
	//  Output:
	//
	//    T t[n+1], y[m*(n+1)]: the times and solution values.
	//
{
	T dt;
	T* f0;
	T* f1;
	T* f2;
	T* f3;
	size_t i;
	size_t j;
	T t0;
	T t1;
	T t2;
	T t3;
	T* u0;
	T* u1;
	T* u2;
	T* u3;

	f0 = new T[m];
	f1 = new T[m];
	f2 = new T[m];
	f3 = new T[m];
	u0 = new T[m];
	u1 = new T[m];
	u2 = new T[m];
	u3 = new T[m];

	dt = (tspan[1] - tspan[0]) / (T)(n);

	j = 0;
	t[0] = tspan[0];
	for (i = 0; i < m; i++)
	{
		y[i + j * m] = y0[i];
	}

	for (j = 0; j < n; j++)
	{
		t0 = t[j];
		for (i = 0; i < m; i++)
		{
			u0[i] = y[i + j * m];
		}
		dyDT(t0, u0, f0);

		t1 = t0 + dt / 2.0;
		for (i = 0; i < m; i++)
		{
			u1[i] = u0[i] + dt * f0[i] / 2.0;
		}
		dyDT(t1, u1, f1);

		t2 = t0 + dt / 2.0;
		for (i = 0; i < m; i++)
		{
			u2[i] = u0[i] + dt * f1[i] / 2.0;
		}
		dyDT(t2, u2, f2);

		t3 = t0 + dt;
		for (i = 0; i < m; i++)
		{
			u3[i] = u0[i] + dt * f2[i];
		}
		dyDT(t3, u3, f3);

		t[j + 1] = t[j] + dt;
		for (i = 0; i < m; i++)
		{
			y[i + (j + 1) * m] = u0[i] + dt * (f0[i] + 2.0 * f1[i] + 2.0 * f2[i] + f3[i]) / 6.0;
		}
	}
	/*
		Free memory.
	*/
	delete[] f0;
	delete[] f1;
	delete[] f2;
	delete[] f3;
	delete[] u0;
	delete[] u1;
	delete[] u2;
	delete[] u3;

	return;
}

//****************************************************************************80
template <typename T>
void EULER(void dyDT(T x, T y[], T yp[]), T tspan[2],
	T y0[], const int n, const size_t m, T t[], T y[])

	/******************************************************************************/
	//
	//  Purpose:
	//
	//    euler approximates the solution to an ODE using Euler's method.
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    25 April 2020
	//
	//  Author:
	//
	//    John Burkardt
	//
	//  Input:
	//
	//    dydt: a function that evaluates the right hand side of the ODE.
	//
	//    T tspan[2]: contains the initial and final times.
	//
	//    T y0[m]: a column vector containing the initial condition.
	//
	//    int n: the number of steps to take.
	//
	//    int m: the number of variables.
	//
	//  Output:
	//
	//    T t[n+1], y[m*(n+1)]: the times and solution values.
	//
{
	T dt;
	T* dy;
	int i;
	size_t j;
	T tfirst;
	T tlast;

	dy = new T[m];

	tfirst = tspan[0];
	tlast = tspan[1];
	dt = (tlast - tfirst) / (T)(n);
	j = 0;
	t[j] = tspan[0];
	for (i = 0; i < m; i++)
	{
		y[i + j * m] = y0[i];
	}

	for (j = 0; j < n; j++)
	{
		dyDT(t[j], y + j * m, dy);
		t[j + 1] = t[j] + dt;
		for (i = 0; i < m; i++)
		{
			y[i + (j + 1) * m] = y[i + j * m] + dt * dy[i];
		}
	}

	delete[] dy;

	return;
}

//****************************************************************************80


template <typename T>
void EULER_mc(void dyDT(T x, T y[], T yp[]), T tspan[2],
	T y0[], const int n, const size_t m, T t[], T y[])
{
	T dt;
	T* dy;
	int i;
	size_t j;
	T tfirst;
	T tlast;
	mxws mxws;

	dy = new T[m];
	T* mc;
	mc = new T[m];

	tfirst = tspan[0];
	tlast = tspan[1];
	
	dt = (tlast - tfirst) / (T)(n);
	j = 0;
	t[j] = tspan[0];
	for (i = 0; i < m; i++)
	{
		y[i + j * m] = y0[i];
	}

	const int totali = 1000;
	for (j = 0; j < n; j++)
	{
		
		for (int i = 0; i < m; i++)
			mc[i] = 0;

		for (int x = 0; x < totali; x++)
		{
			for (i = 0; i < m; i++) 
			{
				dyDT(mxws(t[j], t[j] + dt), y + j * m, dy);	
				mc[i] += dy[i];
			}
		}
		
		t[j + 1] = t[j] + dt;
		
		for (i = 0; i < m; i++)
		{
			mc[i] /= T(totali);
			y[i + (j + 1) * m] = y[i + j * m] + dt * mc[i];
		}
	}

	delete[] dy;
	delete[] mc;

	return;
}


/******************************************************************************/

//https://people.sc.fsu.edu/~jburkardt/cpp_src/midpoint_explicit/midpoint_explicit.cpp

void midpoint_explicit(double* dydt(double x, double y[]),
	double tspan[2], double y0[], int n, int m, double t[], double y[])

	/******************************************************************************/
	/*
		Purpose:

			midpoint_explicit uses an explicit midpoint method to solve an ODE.

		Licensing:

			This code is distributed under the GNU LGPL license.

		Modified:

			06 April 2021

		Author:

			John Burkardt

		Input:

			dydt: a function that evaluates the right hand side of the ODE.

			double tspan[2]: contains the initial and final times.

			double y0[m]: a column vector containing the initial condition.

			int n: the number of steps to take.

			int m: the number of variables.

		Output:

			double t[n+1], y[m*(n+1)]: the times and solution values.
	*/
{
	double dt;
	double* f;
	size_t i;
	size_t j;
	double tm;
	double* ym;

	ym = new double[m];

	dt = (tspan[1] - tspan[0]) / (double)(n);

	t[0] = tspan[0];
	j = 0;
	for (i = 0; i < m; i++)
	{
		y[i + j * m] = y0[i];
	}

	for (j = 0; j < n; j++)
	{
		f = dydt(t[j], y + j * m);
		tm = t[j] + 0.5 * dt;
		for (i = 0; i < m; i++)
		{
			ym[i] = y[i + j * m] + 0.5 * dt * f[i];
		}
		delete[] f;

		f = dydt(tm, ym);
		t[j + 1] = t[j] + dt;
		for (i = 0; i < m; i++)
		{
			y[i + (j + 1) * m] = y[i + j * m] + dt * f[i];
		}
		delete[] f;
	}

	delete[] ym;

	return;
}
//****************************************************************************80
template <typename T>
void midpoint_explicit2(void dyDT(T x, T y[], T yp[]),
	T tspan[2], T y0[], const int n, const size_t m, T t[], T y[])
{
	T dt;

	size_t i;
	size_t j;

	T* ym;
	ym = new T[m];

	T* dy;
	dy = new T[m];

	T tm;

	dt = (tspan[1] - tspan[0]) / (T)(n);

	t[0] = tspan[0];
	j = 0;
	for (i = 0; i < m; i++)
	{
		y[i + j * m] = y0[i];
	}

	for (j = 0; j < n; j++)
	{
		dyDT(t[j], y + j * m, dy);

		for (i = 0; i < m; i++)
		{
			ym[i] = y[i + j * m] + 0.5 * dt * dy[i];
		}
	
		tm = t[j] + 0.5 * dt;

		dyDT(tm, ym, dy);
		
		t[j + 1] = t[j] + dt;
		for (i = 0; i < m; i++)
		{
			y[i + (j + 1) * m] = y[i + j * m] + dt * dy[i];
		}
	
	}

	delete[] ym;
	delete[] dy;

	return;
}


template <typename T>
void midpoint_implicit(void dyDT(T x, T y[], T yp[]),
	T tspan[2], T y0[], const int n, const size_t m, T t[], T y[])
	//implicit midpoint method 
	//https://en.wikipedia.org/wiki/Midpoint_method
{
	T dt;

	size_t i;
	size_t j;

	T* ym;
	ym = new T[m];

	T* dy;
	dy = new T[m];

	T tm;

	dt = (tspan[1] - tspan[0]) / (T)(n);

	t[0] = tspan[0];
	j = 0;
	for (i = 0; i < m; i++)
	{
		y[i + j * m] = y0[i];
	}

	for (j = 0; j < n; j++)
	{
		tm = t[j] + dt;
		dyDT(tm, y + j * m, dy);
		
		for (i = 0; i < m; i++)
		{
			ym[i] = 0.5 * (y[i + j * m] + (y[i + j * m] + dt * dy[i]));
		}
	
		tm = t[j] + 0.5 * dt;
		dyDT(tm, ym, dy);
		
		t[j + 1] = t[j] + dt;
		for (i = 0; i < m; i++)
		{
			y[i + (j + 1) * m] = y[i + j * m] + dt * dy[i];
		}
	
	}

	delete[] ym;
	delete[] dy;

	return;
}


template <typename T>
void rkf(void dyDT(T x, T y[], T yp[]),
	T tspan[2], T y0[], const int n, const size_t m, T t[], T y[])
{
	T dt;

	size_t i;
	size_t j;

	T* dy;
	dy = new T[m];

	T* k1,*k2,*k3,*k4,*k5,*k6;
	k1 = new T[m];
	k2 = new T[m];
	k3 = new T[m];
	k4 = new T[m];
	k5 = new T[m];
	k6 = new T[m];

	dt = (tspan[1] - tspan[0]) / (T)(n);

	t[0] = tspan[0];
	j = 0;
	for (i = 0; i < m; i++)
	{
		y[i + j * m] = y0[i];
	}

	for (j = 0; j < n; j++)
	{
		
		for (i = 0; i < m; i++)
		{
			auto w = y[i + j * m];
			
			auto h = dt;
			
			dyDT(t[j], &w, dy);
			k1[i] = dt * dy[i];
			
			auto W = (w + k1[i] * 0.25);
			dyDT((t[j] + h * 0.25), &W , dy);
			k2[i] = dt * dy[i];
			W = (w + k1[i] * 3. / 32. + k2[i] * 9. / 32.);
			dyDT((t[j] + h * 3. / 8.), &W, dy);
			
			k3[i] = dt * dy[i];
			W = (w + k1[i] * 1932. / 2197. - k2[i] * 7200. / 2197. + k3[i] * 7296. / 2197.);
			dyDT((t[j] + h * 12. / 13), &W, dy);
			k4[i] = dt * dy[i];
			W = (w + k1[i] * 439. / 216. - k2[i] * 8. + k3[i] * 3680. / 513. - k4[i] * 845. / 4104.);
			dyDT((t[j] + h), &W, dy);
			k5[i] = dt * dy[i];
			W = (w - k1[i] * 8. / 27. + k2[i] * 2. - k3[i] * 3544. / 2565. + k4[i] * 1859. / 4104. - k5[i] * 11. / 40.);
			dyDT((t[j] + h * .5), &W , dy);
			k6[i] = dt * dy[i];

		}
		
		t[j + 1] = t[j] + dt;		
		for (i = 0; i < m; i++)
		{
			y[i + (j + 1) * m] = y[i + j * m] + k1[i] * 25. / 216. + k3[i] * 1408. / 2565. + k4[i] * 2197. / 4104. - k5[i] / 5.;
			
		}

	}

	delete[] k1;
	delete[] k2;
	delete[] k3;
	delete[] k4;
	delete[] k5;
	delete[] k6;
	delete[] dy;

	return;
}

template <typename T>
void midpoint_explicit_mc(void dyDT(T x, T y[], T yp[]),
	T tspan[2], T y0[], const int n, const size_t m, T t[], T y[])
{
	T dt;

	int i;
	size_t j;

	T* ym;
	ym = new T[m];

	T* dy;
	dy = new T[m];
	
	T* mc;
	mc = new T[m];
	
	mxws mxws;

	dt = (tspan[1] - tspan[0]) / (T)(n);

	t[0] = tspan[0];
	j = 0;
	for (i = 0; i < m; i++)
	{
		y[i + j * m] = y0[i];
	}

	const size_t totali = 10000;
	for (j = 0; j < n; j++)
	{
		
		for (int i = 0; i < m; i++)
			mc[i] = 0;

		for (int x = 0; x < totali; x++)
		{
			dyDT(mxws(t[j], t[j] + dt), y + j * m, dy);

			for (i = 0; i < m; i++)
				mc[i] += dy[i];
		}

		for (i = 0; i < m; i++)
		{
			mc[i] /= T(totali);
			//ym[i] = y[i + j * m] + 0.5 * dt * dy[i];
			ym[i] = y[i + j * m] + 0.5 * dt * mc[i];
			mc[i] = 0;
		}
	
		//auto tm = t[j] + 0.5 * dt;
		//dyDT(tm, ym, dy);
		
		for (int x = 0; x < totali; x++)
		{
			dyDT(mxws(t[j], t[j] + dt), ym, dy);

			for (i = 0; i < m; i++)
				mc[i] += dy[i];
		}

		t[j + 1] = t[j] + dt;
		for (i = 0; i < m; i++)
		{
			mc[i] /= T(totali);
			y[i + (j + 1) * m] = y[i + j * m] + dt * mc[i];
		}

	}

	delete[] ym;
	delete[] dy;
	delete[] mc;

	return;
}

void Generalized_midpoint3
(
	void dyDT(int t)
)
{
}

typedef void (*functiontype)(int t);

void Generalized_midpoint_ode2
(
	void dyDT(int t)
)
{	
	//functiontype func = dyDT;
	auto func = dyDT;
	Generalized_midpoint3(func);
}


std::vector< std::array<long double, 3> > rkf_calon(
	long double a, 
	long double b, 
	long double alpha, 
	long double TOL, 
	long double hmin, 
	long double hmax, 
	long double (*f)(long double, long double), 
	std::ostream& os) {

	long double h, del, R, t, w;
	long double k1, k2, k3, k4, k5, k6;
	bool flag = true, fail = false;
	std::vector< std::array<long double, 3> > out;
	std::array<long double, 3> temp = { 0,0,0 };

	t = a;
	w = alpha;
	h = hmax;
	temp[0] = t; temp[1] = w; temp[2] = h;
	out.emplace_back(temp);

	while (flag) {
		k1 = h * f(t, w);
		k2 = h * f((t + h * 0.25), (w + k1 * 0.25));
		k3 = h * f((t + h * 3 / 8), (w + k1 * 3 / 32 + k2 * 9 / 32));
		k4 = h * f((t + h * 12 / 13), (w + k1 * 1932 / 2197 - k2 * 7200 / 2197 + k3 * 7296 / 2197));
		k5 = h * f((t + h), (w + k1 * 439 / 216 - k2 * 8 + k3 * 3680 / 513 - k4 * 845 / 4104));
		k6 = h * f((t + h * .5), (w - k1 * 8 / 27 + k2 * 2 - k3 * 3544 / 2565 + k4 * 1859 / 4104 - k5 * 11 / 40));

		R = abs(k1 / 360 - k3 * 128 / 4275 - k4 * 2197 / 75240 + k5 / 50 + k6 * 2 / 55) / h;

		if (R <= TOL) {
			t = t + h;
			w += k1 * 25 / 216 + k3 * 1408 / 2565 + k4 * 2197 / 4104 - k5 / 5;
			temp[0] = t; temp[1] = w; temp[2] = h;
			out.push_back(temp);
		}

		del = 0.84 * pow((TOL / R), .25);
		if (del <= .1) h *= .1;
		else if (del >= 4.0) h *= 4.0;
		else h *= del;

		if (h > hmax) h = hmax;

		if (t >= b) flag = false;
		else if (t + h > b) h = b - t;
		else if (h < hmin) {
			flag = false;
			os << "minimum \'h\' exceeded, failure";
			fail = true;
		}

	}
	os << endl;
	if (fail) out.clear();

	return out;
}

//0.4596976941
//0.4596976941 31860
//0.459697694131860

//exp(-1/t^2)
//wolfram fail : series(exp(-1 / sin(t) ^ 2))

//z = f(x,y) = (x + y) (x y + x y^2)
//∂z / ∂x((x + y) (x y ^ 2 + x y)) = y(y + 1) (2 x + y) == y(2x + y) (y + 1)
//∂z / ∂y((x + y) (x y ^ 2 + x y)) = x(2 x y + x + y(3 y + 2)) == x(3 y ^ 2 + 2y(x + 1) + x)

//Reduce[D[(x + y) (x y + x y^2), {{x, y}}] == 0, {x, y}]
//(x == 0 && y == -1) || (x == 1 && y == -1) || (x == 3/8 && y == -3/4) || (x == 0 && y == 0)


