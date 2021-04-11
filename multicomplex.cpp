
#include "stdafx.hpp"
#include "fast_hadamard_transform_test.hpp"

#include <cmath>

int64_t next_dyck_word
(
	int64_t w
);
 
template<typename T>
std::vector<std::complex<T>> mul_tx
(
	std::vector<std::complex<T>>&,
	std::vector<std::complex<T>>&
);

template<typename A, typename T>
A Binomial_Coefficient(const T& n, const T& k);

const REAL x_max = 2 * pi;
const REAL x_min = 0;

void travellingSalesmanProblem_driver();
void rankOfMatrix_driver();
void Ackermann();

int main(int argc, char** argv)
{
	//diverse tests :
	//rankOfMatrix_driver();
	
	//fast_hadamard_transform_driver<MX0>();
	{
		Matrix<REAL> a(2,2);
		Matrix<REAL> b(2,2);

		a.matrix[0][0] = 6;
		a.matrix[0][1] = -7;
		a.matrix[1][0] = 0;
		a.matrix[1][1] = 3;

		b.matrix[0][0] = 1;
		b.matrix[0][1] = 2;
		b.matrix[1][0] = 3;
		b.matrix[1][1] = 5;

		Matrix<REAL>  c = a / b;
		
		std::cout << c;
	}
	std::cout << std::endl;
	
	std::cout.setf(std::ios::fixed, std::ios::floatfield);
	std::cout.precision(12);

	{
		MX0 x(0.5,1);
		//x = -0.5;
		std::cout << "Riemann Zeta(x) = " << Riemann_Zeta(x) << std::endl << std::endl;
		//std::cout << 1/x << std::endl << std::endl;

	}

	{
		MX0 x(0.4, -0.5);
		MX5 y;
		sh(y, x);


		std::cout << "d^5/dx^5(sin(sqrt(y))) = -16.4786 - 18.9303 i, x = 0.4 - 0.5i =\n ";
		
		auto begin = std::chrono::steady_clock::now();
		
		auto d = dv(sin(sqrt(y)));

		auto end = std::chrono::steady_clock::now();
		std::cout << "\nduration d^5/dx^5(sin(sqrt(y))) " <<
			std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << " uS\n" << std::endl;

		std::cout << d << std::endl << std::endl;

		MX0 a = { -1,0 };
		std::cout << log(a) << std::endl << std::endl; //i pi
	}

	MX2 aa = { {{3,0},{4,0}},{{6,0},{8,0}} };
	std::cout << 1 / aa << std::endl << std::endl;


	Matrix<MX0> qg(2, 2);
	qg = { {1,1},{1,-1} };
	MX0 qg_in{ 0,1 };

	//std::cout << qg_in * (1 / sqrt(2)) * qg << std::endl;

	//hep_driver(x_min,x_max);//vegas integrator

	//auto fx = [](const auto& x) {auto d = di(x); auto k = d-sin(d); auto k2 = 1-cos(d);
	auto fx = [](const auto& x)
	{
		auto d = di(x); auto k = d - sin(d); auto k2 = 1 - cos(d);
		auto k3 = dv(k);
		auto k4 = dv(k2);
		k3 *= k3;
		k4 *= k4;
		return sqrt(k3 + k4);
	};

	std::cout << "Generalized midpoint rule formula\n";
	std::cout << midpoint(fx, MX0{ x_min }, MX0{ x_max }) << "\n\n";

	std::cout << "pi test log(-1)\n";
	std::cout << log(MX0{ -1 }).imag << "\n\n";

	integrator<REAL, 0> qt;

	//std::cout << "integral 1-exp(x) dx from x=0 to pi" << std::endl;/
	//std::cout << "integral 1/sqrt(1-x*x) dx from x = -1 to 1" << std::endl;
	//std::cout << "integral (cos(x)-1)/x" << std::endl;
	//auto f = [](const auto& x) { return sin(log(1-x*x)); };
	//auto f = [](const auto& x) { return (cos(x)-1)/x; };
	//auto f = [](const auto& x) { return sin(x)/x; };

	auto fx2 = [](const auto& x)
	{
		return (cos(x) - 1) / x;//-0.23974116548
	};

	std::cout << "tanh sinh quadrature\n";
	std::cout << qt.ix(/*function*/ fx2, /*a*/{ 0,0 }, /*b*/{ 1,0 }) << std::endl << std::endl;
	//std::cout << qt.ix(/*function*/ f, /*a*/ {-1,0}, /*b*/ {1,0}) << std::endl << std::endl;
	//std::cout << euler + log(xcin) + qt.ix(/*function*/ f, /*a*/ {0,0}, /*b*/ {2,0}) << std::endl << std::endl;

	//std::cout << "Feigenbaum constant\n";
	//std::cout << Feigenbaum<REAL,0>() << std::endl << std::endl;
	MX0 p = { .2, 0.9 };

	std::cout << Riemann_Siegel_Z2(p) << std::endl;

	auto rz = [](const auto& x) { return Riemann_Siegel_Z2(x); };

	root(rz, p, 50);

	//std::cout << std::endl << std::endl;
	//Logistic_Map(p,10);

	std::cout << std::endl << std::endl;

	std::bitset<8> xb(next_dyck_word(0b1010));

	std::cout << "next_dyck_word " << xb << "\n\n";

	{
		{
			MX0 x;
			sh(x, { 2,0 });
			std::cout << "ps::Hermite(0,2) " << dv(ps::Hermite(0, x)) << std::endl;
			std::cout << "ps::Hermite(1,2) " << dv(ps::Hermite(1, x)) << std::endl;
			std::cout << "ps::Hermite(2,2) " << dv(ps::Hermite(2, x)) << std::endl;
			std::cout << "ps::Hermite(3,2) " << dv(ps::Hermite(3, x)) << std::endl << std::endl;
		}

		const int n = 3;
		
		Matrix<REAL> A(n, n), I(n, n), t1(n, n), t2(n, n);
		I.identity();

		A = {
			{2,lambda,lambda},
			{lambda,4,5 },
			{-1,4,-3.353422543}
		};

		auto a = A.tr();
		auto b = a * a - (A * A).tr();

		auto determinant = det(A, n);

		std::cout << "Characteristic polynomial\n" << a << " " << -0.5 * b << " " << determinant << std::endl;

		//2.646577457000 32.120535258000 -66.827380344000 = coefficients
		//=-λ^3+2.646577457λ^2+32.12053\dots λ-66.82738
		
		std::cout << std::endl;

		//eigen vector
		Vec bk{ 1,1,1 };
		auto r = A - (2 * I);
		determinant = det(r, n);
		r.inverse(t1, n, determinant);
		// std::cout << t1;
	
		MX0 sx{ -5.882346754565423L,0 };
		auto fx1 = [](const auto& x) { return pow(-x, 3) + 2.646577457 * pow(x, 2) + 32.120535258000 * x - 66.827380344000; };

		//2.00000 -5.46620 6.11277 roots
		std::cout << bk;
		std::cout << std::endl;

		root(fx1, sx, 15);
		std::cout << std::endl;

		//Inverse iteration
		// is an iterative eigenvalue algorithm. 
		//It allows one to find an approximate eigenvector when an approximation to a corresponding eigenvalue is already known

		bk = bk * t1;
		bk = normalize(bk);
		// 2.00000 -5.46620 6.11277

		std::cout << "Eigen value 2 " << "Eigen vector " << bk;
	
	}

	std::cout << std::endl;

	//Linked_List_introduction<MX0,MX1,MX5,MX5>(); 

	//Construct_Complete_Binary_Tree_from_its_Linked_List_Representation<MX0,MX1,MX5,MX5>();
	//std::cout << std::endl;

	//flatten_binary_tree<MX0,MX1>();
	//std::cout << std::endl;

	MX4 x;//multicomplex derivative order

	//std::vector<MX15> xv(32);      
	std::cout << std::endl;
	//xv[31] *= 6; 
	//- 0.165395312473            

	multicomplex<REAL, 0> y{ 0.7,-0.5 };
	const int dv = 0;//multicomplex order

	mcdv mcdv;
	//mcdv.sh<dv>(x, {{{{-1,0.7},{-0.7,-0.6}},{{-1,0.5},{-0.7,-0.6}}},{{{-1,0.23446},{-0.7,-0.0865}},{{-1.654,0.5},{0.7,-0.6}}}});
	//mcdv.sh<dv>(x, {{{-1,0.7},{-0.7,-0.6}},{{-1,0.5},{-0.7,-0.6}}});
	mcdv.sh<dv>(x, { 2,0 });

	std::cout << "multicomplex derivatives..." << std::endl;

	//auto jb = mcdv.dv<dv>(sqrt(x/y));
	//std::cout <<  jb*jb << std::endl << std::endl;

	auto begin = std::chrono::steady_clock::now();
	auto kb = mcdv.dv<dv>(sin(x / y));

	auto end = std::chrono::steady_clock::now();
	std::cout << "\nDuration sin(x/y) " <<
		std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << " uS\n" << std::endl;
	std::cout << kb << std::endl;

	//util::StartCounter();
	begin = std::chrono::steady_clock::now();
	kb = mcdv.dv<dv>(Sin(x / y));
	//std::this_thread::sleep_for(std::chrono::microseconds(100) );

	end = std::chrono::steady_clock::now();
	std::cout << "\nDuration Sin(x/y) " <<
		std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << " uS\n" << std::endl;
	std::cout << kb << std::endl << std::endl;

	MX0 x2{ 1.5,0 };
	//sh(x2, {1.5,0});
	//auto D1 = gamma(x2);
	std::cout << "gamma" << std::endl;
	std::cout << gamma(x2) << std::endl;
	std::cout << std::tgamma(1.5) << std::endl << std::endl;

	std::cout.precision(13);
	std::cout << "Exp integral " << std::endl;
	MX0 xcin = { 2, 0 };
	//sh(xcin, {1.5,0});
	std::cout << std::expint(2) << std::endl;
	std::cout << ps::Ei(xcin) << std::endl << std::endl;

	std::cout.precision(13);
	std::cout << "Exp integral N1 " << std::endl;
	std::cout << ps::E1(xcin) << std::endl << std::endl;

	std::cout << "Log integral " << std::endl;
	std::cout << ps::li(xcin) << std::endl << std::endl;

	MX0 xc1 = { 0, 2 };
	MX0 xc2 = { 0,-2 };
	MX0 xc3 = { 0, 2 };

	std::cout << "Sine integral " << std::endl;
	//std::cout << (1 / xc3) * (ps::E1(xc1) - ps::E1(xc2)) + half_pi << std::endl << std::endl;

	std::cout << "Cosine integral " << std::endl;
	//std::cout << half * (ps::Ei(xc1) + ps::Ei(xc2)) << std::endl << std::endl;

	std::cout << "Sinh integral " << std::endl;
	std::cout << ps::Shi(xcin) << std::endl << std::endl;

	std::cout << "Cosh integral " << std::endl;
	std::cout << ps::Chi(xcin) << std::endl << std::endl;

	std::cout.precision(10);
	std::cout << "atan" << std::endl;
	std::cout << std::atan(1.5) << std::endl;
	MX0 xa{ 1.5,0 };
	std::cout << atan(xa) << std::endl << std::endl;

	std::cout.precision(12);
	//MX0 sx {7,0};//(+ 1.1673039783 + 0.0000000000*i1)//initial guess
	MX2 sx;//(+ 1.1673039783 + 0.0000000000*i1)//initial guess
	//MX0 sx {0,11};//-ProductLog[-3,1]  {4, 3} [1.53391331979357 4.37518515306190] -ProductLog[-1,1] 
	//std::cout << -ps::LambertW(4, MX0{1,0}) << std::endl;//exp(x) + x
	//std::cout << -ps::LambertW(4, MX0{-1,0}) << std::endl;//exp(x) - x
	//std::cout << ps::LambertW(4, MX0{-1,0}) << std::endl;//exp(-x) + x
	std::cout << "LambertW(4, MX0{1,0}) = " << LambertW(4, MX0{-1, 2}) << std::endl;//exp(-x) - x

	std::cout << std::endl;

	sx.random(-10.,10.);

	std::cout << "multicomplex roots..." << std::endl;
	std::cout << "ini" << sx << std::endl << std::endl;

	//auto fx1 = [](const auto& x) { return Wilkinsons_polynomial(x, 20); };
	auto fx1 = [](const auto& x) { return 2. * pow(x, 2.) - 10. * x + 5.; };

	root(fx1, sx, 15);
	
	std::cout << std::endl;

	//(+ 1.531754163448*j0 + 0.000000000000*j1 - 0.000000000000*j2 - 0.968245836552*j3 - 0.000000000000*j4 + 0.968245836552*j5 + 0.968245836552*j6 - 0.000000000000*j7)  
	//(+ 1.531754163448*j0 - 0.000000000000*j1 + 0.000000000000*j2 + 0.968245836552*j3 + 0.000000000000*j4 + 0.968245836552*j5 - 0.968245836552*j6 + 0.000000000000*j7
	//(+ 1.531754163448*j0 + 0.000000000000*j1 + 0.000000000000*j2 - 0.968245836552*j3 - 0.000000000000*j4 - 0.968245836552*j5 - 0.968245836552*j6 - 0.000000000000*j7)
	//(+ 1.531754163448*j0 + 0.000000000000*j1 - 0.000000000000*j2 + 0.968245836552*j3 + 0.000000000000*j4 - 0.968245836552*j5 + 0.968245836552*j6 + 0.000000000000*j7)

	//(+ 3.468245836552*j0 + 0.000000000000*j1 + 0.000000000000*j2 + 0.968245836552*j3 - 0.000000000000*j4 + 0.968245836552*j5 + 0.968245836552*j6 + 0.000000000000*j7)  
	//(+ 3.468245836552*j0 + 0.000000000000*j1 + 0.000000000000*j2 + 0.968245836552*j3 - 0.000000000000*j4 - 0.968245836552*j5 - 0.968245836552*j6 + 0.000000000000*j7)
	//(+ 3.468245836552*j0 + 0.000000000000*j1 + 0.000000000000*j2 - 0.968245836552*j3 + 0.000000000000*j4 + 0.968245836552*j5 - 0.968245836552*j6 + 0.000000000000*j7)
	//(+ 3.468245836552*j0 + 0.000000000000*j1 - 0.000000000000*j2 - 0.968245836552*j3 - 0.000000000000*j4 - 0.968245836552*j5 + 0.968245836552*j6 - 0.000000000000*j7)

	//(+ 2.500000000000*j0 + 0.000000000000*j1 + 0.000000000000*j2 + 0.000000000000*j3 + 0.000000000000*j4 + 0.000000000000*j5 + 1.936491673104*j6 + 0.000000000000*j7)
	//(+ 2.500000000000*j0 + 0.000000000000*j1 + 0.000000000000*j2 + 0.000000000000*j3 + 0.000000000000*j4 - 1.936491673104*j5 + 0.000000000000*j6 + 0.000000000000*j7)
	//(+ 2.500000000000*j0 + 0.000000000000*j1 + 0.000000000000*j2 + 1.936491673104*j3 + 0.000000000000*j4 + 0.000000000000*j5 + 0.000000000000*j6 + 0.000000000000*j7)
	//(+ 2.500000000000*j0 + 0.000000000000*j1 + 0.000000000000*j2 - 1.936491673104*j3 + 0.000000000000*j4 + 0.000000000000*j5 + 0.000000000000*j6 + 0.000000000000*j7)


	//(+ 0.563508326896*j0 + 0.000000000000*j1 + 0.000000000000*j2 + 0.000000000000*j3 + 0.000000000000*j4 + 0.000000000000*j5 + 0.000000000000*j6 + 0.000000000000*j7)
	//(+ 4.436491673104*j0 + 0.000000000000*j1 + 0.000000000000*j2 + 0.000000000000*j3 + 0.000000000000*j4 + 0.000000000000*j5 + 0.000000000000*j6 + 0.000000000000*j7)

	//std::cout << "start x " << sx << "\n" << root(sx) << std::endl << std::endl;
	//[0.00000000000000        12.56637061435917]
	//(-2.3058910256099833e-17,12.566370614359172)

	//component notation
	std::cout.precision(7);

	MX0 a{ +1,-0.5 };
	MX0 b{ 1,0 };
	MX1 c{ {+1, -2}, {8, -1} };
	MX1 cc{ {+1, 2}, {-.5, 0.6} };
	MX2 d{ {{+1,-2},{8, -1}},{{0.8,2},{1,0.7}} };
	MX2 e{ {{+1, 3},{-.5, 0.6}},{{2,0},{1,0}} };

	//(+ 10.400*j0 - 4.300*j1 + 11.700*j2 + 24.600*j3 - 8.000*j4 + 1.000*j5 + 1.000*j6 - 2.000*j7)

	MX0 kj = c.real + c.real - cc.real;
	MX1 pp{ { kj - c.imag * cc.imag},{ c.real * c.real} };
	//std::cout << pp << "\n\n";
	//std::cout << c * cc << "\n\n";

	std::cout << "\nabs(e) " << abs(e);
	std::cout << "\n\n";

	std::vector<std::complex<long double>> D{ {+1,-2},{8, -1},{0.8,2},{1,0.7} };
	std::vector<std::complex<long double>> E{ {+1, 3},{-.5, 0.6},{2,0},{1,0} };

	begin = std::chrono::steady_clock::now();

	D = mul_tx(D, E);

	end = std::chrono::steady_clock::now();

	std::cout.precision(3);
	for (auto& i : D)
		std::cout << i << " ";

	std::cout << "\nDuration D = mul_tx(D,E) " <<
		std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() << " nS\n" << std::endl;

	begin = std::chrono::steady_clock::now();
	d *= e;

	end = std::chrono::steady_clock::now();;

	std::cout << d;
	std::cout << "\nDuration D = D*E " <<
		std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() << " nS\n" << std::endl;

	c = c - 8;

	//std::cout << "\nc = c = c/a " << c << "\n\n";

	c = { {+1, -2}, {8,  -1} };

	a = { +1,0.5 };
	c -= 8;

	//a =  {+1,0.5};
	//d.at_var<0>(2) = {8,9};
	//a = 1;

	//std::cout << "c /=  a " << cc/d << "\n\n";

	//if(a < b)
	//std::cout << "\n true\n" << std::endl;
	//else std::cout << "\n false\n" << std::endl;

	//a.at_var<0>(0) = {5,5};
	//std::cout << a.at_con<0>(0) << "\n\n";
	//c = (c) / 1;

	std::cout.precision(13);


	std::cout << sqrt(x) << std::endl << std::endl; 
	std::cout << "GFG: "
	       << std::is_standard_layout<MX0>::value << '\n'; 

	//std::cout << funStruct<8>::val << std::endl; 
	//std::cout << Factorial<30>::value << "\n";

	std::cout.precision(4);
	VX::vc_eval<REAL>();

	VX::driver_fft<double>();
	return 0;
}

//unrolled 
template<typename T>
std::vector<std::complex<T>> mul_tx
(
	std::vector<std::complex<T>>& v1,
	std::vector<std::complex<T>>& v2
)
{
	//std::vector<std::complex<long double>> D {{+1,-2},{8, -1},{0.8,2},{1,0.7}};
	//std::vector<std::complex<long double>> E {{+1, 3},{-.5, 0.6},{2,0},{1,0}};
	//real * o.real - imag * o.imag | |

	std::vector<std::complex<T>> V(v1.size()), V2(v1.size());

	std::complex<T> a, b;

	//real
	a = v1[0] * v2[0];
	b = v1[1] * v2[1];
	V[0] = a - b;

	//imag
	a = v1[0] * v2[1];
	b = v1[1] * v2[0];
	V[1] = a + b;

	//real * o.real - imag * o.imag,
	//real * o.imag + imag * o.real

	//real
	a = v1[2] * v2[2];
	b = v1[3] * v2[3];
	V[2] = a - b;

	//imag
	a = v1[2] * v2[3];
	b = v1[3] * v2[2];
	V[3] = a + b;

	V[0] -= V[2];
	V[1] -= V[3];

	//////////////

	//std::vector<std::complex<long double>> D {{+1,-2},{8, -1},{0.8,2},{1,0.7}};
	//std::vector<std::complex<long double>> E {{+1, 3},{-.5, 0.6},{2,0},{1,0}};
	//real * o.imag + imag * o.real X X 

	//real
	a = v1[0] * v2[2];
	b = v1[1] * v2[3];
	V2[0] = a - b;

	//imag
	a = v1[0] * v2[3];
	b = v1[1] * v2[2];
	V2[1] = a + b;

	//real
	a = v1[2] * v2[0];
	b = v1[3] * v2[1];
	V2[2] = a - b;

	//imag
	a = v1[2] * v2[1];
	b = v1[3] * v2[0];
	V2[3] = a + b;

	V2[0] += V2[2];
	V2[1] += V2[3];

	V[2] = V2[0];
	V[3] = V2[1];

	return V;
}

int64_t next_dyck_word
(
	int64_t w
)
{
	int64_t const a = w & -w;
	int64_t const b = w + a;
	int64_t       c = w ^ b;
	c = (c / a >> 2) + 1;
	c = ((int64_t(c * c) - 1) & 0xaaaaaaaaaaaaaaaa) | b;
	return c;
}

int ack(int m, int n)
{
	int ans;
	if (m == 0) ans = n + 1;
	else if (n == 0) ans = ack(m - 1, 1);
	else ans = ack(m - 1, ack(m, n - 1));
	return (ans);
}

void Ackermann()
{
	int i, j;
	for (i = 0; i < 6; i++)
		for (j = 0; j < 6; j++)

			printf("ackerman (%d,%d) is: %d\n", i, j, ack(i, j));
}

