
#include <iostream>    
#include <vector>
#include <functional>
#include <limits>
#include <chrono>
#include <complex>
#include <map>
#include <bitset>
#include <random>


struct fail 
{
	fail(char const* const m)
	{
		std::cout << "\n\n *** fail *** " << m << " *** \n\n"; 
		exit(1);
	}
};


//http://tamivox.org/eugene/multicomplex/index.html
template <typename elem, int order>
class multicomplex
{

	static_assert (order >= 0, "multicomplex order must be nonnegative");
	static_assert (order < 25, "sanity_check -- modify if desired");

private:
	std::int64_t i{};

public:
	
	static constexpr int off{ 2 << (order - 1) };

	using one_less = multicomplex<elem, order - 1>;//recursion

	one_less real;
	one_less imag;

	multicomplex() = default;
	virtual ~multicomplex() = default;

	multicomplex(const one_less& r, const one_less& i) : real(r), imag(i) {}
	multicomplex(const elem& y) : real(y) {}
	
	// unary operators
	multicomplex operator+ () const 
	{
		return 
		{ 
			+real, 
			+imag 
		};
	}
	multicomplex operator~ () const 
	{
		return 
		{ 
			+real, 
			-imag 
		};
	}
	multicomplex operator- () const 
	{
		return 
		{ 
			-real, 
			-imag 
		};
	}

	// addition
	template <int o_order>
	typename std::enable_if <(order > o_order), multicomplex<elem, order>>::type
		operator+ (multicomplex<elem, o_order> const& o) const 
	{
		return 
		{ 
			real + o, 
			imag 
		};
	}

	template <int o_order>
	typename std::enable_if <(order < o_order), multicomplex<elem, o_order>>::type
		operator+ (multicomplex<elem, o_order> const& o) const 
	{
		return 
		{ 
			*this + o.real, 
			+o.imag 
		};
	}

	template <int o_order>
	typename std::enable_if <(order == o_order), multicomplex<elem, order>>::type
		operator+ (multicomplex<elem, o_order> const& o) const 
	{
		return 
		{ 
			real + o.real, 
			imag + o.imag 
		};
	}


	template <typename T>
	multicomplex operator+ (T const& o) const
	{
		return *this + multicomplex(o);
	}

	template <int o_order>
	typename std::enable_if <(order == o_order), multicomplex<elem, order>>::type
		operator+ (elem const& o)
	{
		return
		{
			real + o,
			imag + o
		};
	}

	template <int o_order>
	typename std::enable_if <(order > o_order), multicomplex<elem, order>>::type
		operator+= (multicomplex<elem, o_order> const& o) 
	{
		real += o;

		return 
		{ 
			real, 
			imag 
		};
	}

	template <int o_order>
	typename std::enable_if <(order == o_order), multicomplex<elem, order>>::type
		operator+= (multicomplex<elem, o_order> const& o) 
	{
		return *this = *this + o;
	}
	
	multicomplex operator+= (elem const& o)
	{
		return
		{
			real += o,
			imag += o
		};
	}
	
	// prefix increment operator overloading 
	multicomplex operator++ ()
	{
		++i;
		*this += i;
		return *this;
	}

	// postfix increment operator overloading 
	multicomplex operator++ (int)
	{
		i++;
		return *this += i;
	}

	// subtraction
	template <int o_order>
	typename std::enable_if <(order > o_order), multicomplex<elem, order>>::type
		operator- (multicomplex<elem, o_order> const& o) const 
	{
		return 
		{ 
			real - o, 
			imag 
		};
	}

	template <int o_order>
	typename std::enable_if <(order < o_order), multicomplex<elem, o_order>>::type
		operator- (multicomplex<elem, o_order> const& o) const 
	{
		return 
		{ 
			*this - o.real, 
			-o.imag 
		};
	}

	template <int o_order>
	typename std::enable_if <(order == o_order), multicomplex<elem, order>>::type
		operator- (multicomplex<elem, o_order> const& o) const 
	{
		return 
		{ 
			real - o.real, 
			imag - o.imag 
		};
	}

		template <int o_order>
	typename std::enable_if <(order > o_order), multicomplex<elem, order>>::type
		operator-= (multicomplex<elem, o_order> const& o) 
	{
		real -= o;
		return 
		{ 
			real, 
			imag 
		};
	}

	template <int o_order>
	typename std::enable_if <(order == o_order), multicomplex<elem, order>>::type
		operator-= (multicomplex<elem, o_order> const& o) 
	{
		return *this = *this - o;
	}

	multicomplex operator-= (elem const& o) 
	{
		return *this = *this - o;
	}

	template <typename T>
	multicomplex operator- (T const& o) const 
	{
		return *this - multicomplex(o);
	}

	// multiplication
		template <int o_order>
	typename std::enable_if <(order > o_order), multicomplex<elem, order>>::type
		operator* (multicomplex<elem, o_order> const& o) const 
	{
		return 
		{ 
			real * o, 
			imag * o 
		};
	}

	template <int o_order>
	typename std::enable_if <(order < o_order), multicomplex<elem, o_order>>::type
		operator* (multicomplex<elem, o_order> const& o) const 
	{
		return 
		{ 
			*this * o.real, 
			*this * o.imag 
		};
	}

	template <int o_order>
	typename std::enable_if <(order == o_order), multicomplex<elem, order>>::type
		operator* (multicomplex<elem, o_order> const& o) const 
	{
		return 
		{
			real * o.real - imag * o.imag,
			real * o.imag + imag * o.real
		};
	}

	template <int o_order>
	typename std::enable_if <(order > o_order), multicomplex<elem, order>>::type
		operator*= (multicomplex<elem, o_order> const& o) 
	{
		return *this = *this * o;
	}

	template <int o_order>
	typename std::enable_if <(order == o_order), multicomplex<elem, order>>::type
		operator*= (multicomplex<elem, o_order> const& o) 
	{
		return *this = *this * o;
	}

	multicomplex operator*= (elem const& o) 
	{
		return *this = *this * o;
	}

	template <typename T>
	multicomplex operator* (const T& o) const
	{
		return *this * o;
	}

	multicomplex operator* (elem const& o) const 
	{
		return 
		{ 
			real * o, 
			imag * o 
		};
	}

	// division
			template <int o_order>
	typename std::enable_if <(order >= o_order), multicomplex<elem, order>>::type
		operator/ (multicomplex<elem, o_order> const& o) const 
	{
		auto const new_numer { *this * ~o };
		auto const new_denom { o.fold() };
		return new_numer / new_denom;
	}

	template <int o_order>
	typename std::enable_if <(order < o_order), multicomplex<elem, o_order>>::type
		operator/ (multicomplex<elem, o_order> const& o) const 
	{
		auto const new_numer { *this * ~o };
		auto const new_denom { o.fold() };
		return new_numer / new_denom;
	}

	template <typename T>
	multicomplex operator/ (T const& o) const 
	{
		return *this / multicomplex (static_cast<elem>(o));
	}

	multicomplex operator/ (multicomplex<elem, 0> const& o) const 
	{
		return 
		{ 
			real / o, 
			imag / o 
		};
	}

	multicomplex operator/ (elem const& o) const 
	{
		return 
		{ 
			real / o, 
			imag / o 
		};
	}

	multicomplex operator/= (multicomplex<elem, 0> const& o) 
	{
		real /= o;
		imag /= o;
		return *this;
	}

	one_less fold() const {
		return real * real + imag * imag;
	}

	// miscellaneous
	elem norm() const 
	{
		return real.norm() + imag.norm();
	}
	elem signr() const 
	{
		return real.signr();
	}
	elem signi() const 
	{
		return real.signi();
	}
	elem comp() const 
	{
		return real.comp() + imag.comp();
	}

	// selection by position
	private:

	template <int x_order>
	static typename
		std::enable_if <(order > x_order), multicomplex<elem, x_order> const&>::type
		priv_at_con(multicomplex<elem, order> const& source, int const n) {
		static constexpr int half{ 1 << (order - x_order - 1) };
		if (n < half)
			return one_less::template priv_at_con<x_order>(source.real, n);
		else
			return one_less::template priv_at_con<x_order>(source.imag, n - half);
	}

	template <int x_order>
	static typename
		std::enable_if <(order > x_order), multicomplex<elem, x_order>&>::type
		priv_at_var(multicomplex<elem, order>& source, int const n) {
		static constexpr int half{ 1 << (order - x_order - 1) };
		if (n < half)
			return one_less::template priv_at_var<x_order>(source.real, n);
		else
			return one_less::template priv_at_var<x_order>(source.imag, n - half);
	}

	template <int x_order>
	static typename
		std::enable_if <(order == x_order), multicomplex<elem, x_order>&>::type
		priv_at_var(multicomplex<elem, order>& source, int const n) {
		return source;
	}

	template <int x_order>
	static typename
		std::enable_if <(order == x_order), multicomplex<elem, x_order> const&>::type
		priv_at_con(multicomplex<elem, order> const& source, int const n) {
		return source;
	}

public:
	template <int x_order>
	multicomplex<elem, x_order> const& at_con(int const n) const {
		static_assert (x_order >= 0, "x_order < 0");
		static_assert (x_order <= order, "x_order > order");
		static int const full{ 1 << (order - x_order) };
		if (n < 0) throw fail{ "multicomplex::at_con ... n < 0" };
		if (n >= full) throw fail{ "multicomplex::at_con ... n >= full" };
		return priv_at_con<x_order>(*this, n);
	}

	template <int x_order>
	multicomplex<elem, x_order>& at_var(int const n) {
		static_assert (x_order >= 0, "x_order < 0");
		static_assert (x_order <= order, "x_order > order");
		static int const full{ 1 << (order - x_order) };
		if (n < 0) throw fail{ "multicomplex::at_var ... n < 0" };
		if (n >= full) throw fail{ "multicomplex::at_var ... n >= full" };
		return priv_at_var<x_order>(*this, n);
	}

	template<typename elem2, int order2>
	friend class multicomplex;

	// comparison
	template <int x_order>
	bool operator== (multicomplex<elem, x_order>& b) const
	{
		if (comp() == b.comp()) return true;
		return false;
	}

	template <int x_order>
	bool operator!= (multicomplex<elem, x_order>& b) const
	{
		if (comp() == b.comp()) return false;
		return true;
	}

	template <int x_order>
	bool operator> (multicomplex<elem, x_order>& b) const
	{
		if (comp() > b.comp()) return true;
		return false;
	}

	template <int x_order>
	bool operator< (multicomplex<elem, x_order>& b) const
	{
		if (comp() < b.comp()) return true;
		return false;
	}

	template <int x_order>
	bool operator>= (multicomplex<elem, x_order>& b) const
	{
		if (comp() >= b.comp()) return true;
		return false;
	}

	template <int x_order>
	bool operator<= (multicomplex<elem, x_order>& b) const
	{
		if (comp() <= b.comp()) return true;
		return false;
	}

	multicomplex random(elem lower, elem upper);

	void view_i(std::ostream&) const;
	void view_j(std::ostream&) const;
	void view_j(std::ostream&, const int) const;
	void view_k(std::ostream&) const;

};

/////////////

template <typename elem>
class multicomplex<elem, 0>
{

private:
	std::int64_t i{};

public:

	static constexpr int order_value{ 0 };

	elem real{};
	elem imag{};

	multicomplex() = default;
	virtual ~multicomplex() = default;

	multicomplex(const elem& r, const elem& i) : real(r), imag(i) {}
	multicomplex(const elem& y) : real(y) {}

	// unary operators
	multicomplex operator+ () const 
	{
		return 
		{ 
			+real, 
			+imag 
		};
	}
	multicomplex operator~ () const 
	{
		return 
		{ 
			+real, 
			-imag 
		};
	}
	multicomplex operator- () const 
	{
		return 
		{ 
			-real, 
			-imag 
		};
	}

	// addition
		template <int o_order>
	typename std::enable_if <(0 < o_order), multicomplex<elem, o_order>>::type
		operator+ (multicomplex<elem, o_order> const& o) const 
	{
		return 
		{ 
			*this + o.real, 
			+o.imag 
		};
	}

	template <int o_order>
	typename std::enable_if <(0 == o_order), multicomplex<elem, 0>>::type
		operator+ (multicomplex<elem, o_order> const& o) const 
	{
		return 
		{ 
			real + o.real, 
			imag + o.imag 
		};
	}

	template <int o_order>
	typename std::enable_if <(0 == o_order), multicomplex<elem, 0>>::type
		operator+= (multicomplex<elem, o_order> const& o) 
	{
		return *this = *this + o;
	}

	template <typename T>
	multicomplex operator+ (T const& o) const
	{			
		return
		{
			real + o,
			imag
		};		
	}

	// prefix increment operator overloading 
	multicomplex operator++ ()
	{
		++i;
		*this += i;
		return *this;
	}

	// postfix increment operator overloading 
	multicomplex operator++ (int)
	{
		i++;
		*this += i;
		return *this;
	}

	multicomplex operator+= (elem const& o) 
	{
		return
		{
			real += o,
			imag += o
		};
	}

	// subtraction
		template <int o_order>
		typename std::enable_if <(0 < o_order), multicomplex<elem, o_order>>::type
			operator- (multicomplex<elem, o_order> const& o) const 
		{
			return 
			{ 
				*this - o.real, 
				-o.imag 
			};
		}

		template <int o_order>
		typename std::enable_if <(0 == o_order), multicomplex<elem, 0>>::type
			operator- (multicomplex<elem, o_order> const& o) const 
		{
			return 
			{ 
				real - o.real, 
				imag - o.imag 
			};
		}

		template <int o_order>
		typename std::enable_if <(0 == o_order), multicomplex<elem, 0>>::type
			operator-= (multicomplex<elem, o_order> const& o) 
		{
			return *this = *this - o;
		}

		
		template <typename T, int o_order>
		typename std::enable_if <(0 == o_order), multicomplex<elem, 0>>::type
			operator- (T const& o) const 
		{
			return *this - o;
		}

		template <typename T>
		multicomplex operator- (T const& o) const
		{
			return
			{
				real - o,
				imag
			};
		}

		multicomplex operator-= (elem const& o)
		{
			return
			{
				real -= o,
				imag
			};
		}

	// multiplication
		template <int o_order>
	typename std::enable_if <(0 < o_order), multicomplex<elem, o_order>>::type
		operator* (multicomplex<elem, o_order> const& o) const 
	{
		return 
		{ 
			*this * o.real, 
			*this * o.imag 
		};
	}

	template <int o_order>
	typename std::enable_if <(0 == o_order), multicomplex<elem, 0>>::type
		operator* (multicomplex<elem, o_order> const& o) const 
	{
		return 
		{
			real * o.real - imag * o.imag,
			real * o.imag + imag * o.real
		};
	}

	template <int o_order>
	typename std::enable_if <(0 == o_order), multicomplex<elem, 0>>::type
		operator*= (multicomplex<elem, o_order> const& o) 
	{
		return *this = *this * o;
	}

	template <typename T, int o_order>
	typename std::enable_if <(0 == o_order), multicomplex<elem, 0>>::type
		operator* (const T& o) const
	{
		return *this * o;
	}

	multicomplex operator* (elem const& o) const 
	{
		return 
		{ 
			real * o, 
			imag * o 
		};
	}

	multicomplex operator*= (elem const& o) 
	{
		return *this = *this * o;
	}

	// division
	template <int o_order>
	typename std::enable_if <(0 < o_order), multicomplex<elem, o_order>>::type
		operator/ (multicomplex<elem, o_order> const& o) const
	{
		auto const new_numer{ *this * ~o };
		auto const new_denom{ o.fold() };
		return new_numer / new_denom;
	}

	template <int o_order>
	typename std::enable_if <(0 == o_order), multicomplex<elem, 0>>::type
		operator/ (multicomplex<elem, o_order> const& o) const
	{
		return 
		{
			(real * o.real + imag * o.imag) / (o.real * o.real + o.imag * o.imag),
			(imag * o.real - real * o.imag) / (o.real * o.real + o.imag * o.imag)
		};
	}

	template <int o_order>
	typename std::enable_if <(0 == o_order), multicomplex<elem, 0>>::type
		operator/= (multicomplex<elem, o_order> const& o) 
	{
		*this = *this / o;
		return *this;
	}

	template <typename T>
	multicomplex operator/ (T const& o) const 
	{
		return *this / multicomplex<elem, 0>(static_cast<elem>(o));
	}

	multicomplex operator/ (elem const& o) const 
	{
		return 
		{ 
			real / o, 
			imag / o 
		};
	}

	multicomplex operator/= (elem const& o) 
	{
		return *this = *this / o;
	}

	elem fold()
	{
		return real * real + imag * imag;
	}

	// miscellaneous
	elem signr() const 
	{
		return real;
	}
	elem signi() const 
	{
		return imag;
	}
	elem norm() const 
	{
		return real * real + imag * imag;
	}
	elem comp() const 
	{
		return real * imag + real + imag;
	}

	// selection by position
	private:
	template <int x_order>
	static typename
		std::enable_if <(order_value == x_order), multicomplex<elem, x_order> const&>::type
		priv_at_con(multicomplex<elem, 0> const& source, int const n) {
		return source;
	}

	template <int x_order>
	static typename
		std::enable_if <(order_value == x_order), multicomplex<elem, x_order>&>::type
		priv_at_var(multicomplex<elem, 0>& source, int const n) {
		return source;
	}

	public:
		template <int x_order>
		multicomplex<elem, x_order> const& at_con(int const n) const {
			static_assert (x_order >= 0, "x_order < 0");
			static_assert (x_order <= order_value, "x_order > order");
			static int const full{ 1 << (order_value - x_order) };
			if (n < 0) throw fail{ "multicomplex::at_con ... n < 0" };
			if (n >= full) throw fail{ "multicomplex::at_con ... n >= full" };
			return priv_at_con<x_order>(*this, n);
		}

		template <int x_order>
		multicomplex<elem, x_order>& at_var(int const n) {
			static_assert (x_order >= 0, "x_order < 0");
			static_assert (x_order <= order_value, "x_order > order");
			static int const full{ 1 << (order_value - x_order) };
			if (n < 0) throw fail{ "multicomplex::at_var ... n < 0" };
			if (n >= full) throw fail{ "multicomplex::at_var ... n >= full" };
			return priv_at_var<x_order>(*this, n);
		}

		template<typename elem2, int order2>
		friend class multicomplex;

	// comparison
	template <int x_order>
	bool operator== (multicomplex<elem, x_order>& b) const
	{
		if (comp() == b.comp()) return true;
		return false;
	}

	template <int x_order>
	bool operator!= (multicomplex<elem, x_order>& b) const
	{
		if (comp() == b.comp()) return false;
		return true;
	}

	template <int x_order>
	bool operator> (multicomplex<elem, x_order>& b) const
	{
		if (comp() > b.comp()) return true;
		return false;
	}

	template <int x_order>
	bool operator< (multicomplex<elem, x_order>& b) const
	{
		if (comp() < b.comp()) return true;
		return false;
	}

	template <int x_order>
	bool operator>= (multicomplex<elem, x_order>& b) const
	{
		if (comp() >= b.comp()) return true;
		return false;
	}

	template <int x_order>
	bool operator<= (multicomplex<elem, x_order>& b) const
	{
		if (comp() <= b.comp()) return true;
		return false;
	}

	void view_i(std::ostream&) const;
	void view_j(std::ostream&) const;
	void view_j(std::ostream&, const int) const;
	void view_k(std::ostream&) const;

	multicomplex random(elem lower, elem upper);

};

/////////////////////////////////////////

template <typename old_elem, typename new_elem, int order>
multicomplex<new_elem, order> convert_elem(multicomplex<old_elem, order> const& o) 
{
	return multicomplex<new_elem, order> {
		convert_elem<old_elem, new_elem, order - 1>(o.real),
			convert_elem<old_elem, new_elem, order - 1>(o.imag)
	};
}


//cout
	template <typename elem, int order>
void multicomplex<elem, order>::view_i
(
	std::ostream& o
) const
{
	o << "(" << real << " " << "+ " << imag << "*i" << order + 1 << ")";
}

template <typename elem>
void multicomplex<elem, 0>::view_i
(
	std::ostream& o
) const
{
	o << "(";
	if (real >= 0)
		o << "+ " << +real;
	else
		o << "- " << -real;
	if (imag >= 0)
		o << " + " << +imag;
	else
		o << " - " << -imag;
	o << "*i1" << ")";
}

template <typename elem, int order>
void multicomplex<elem, order>::view_j
(
	std::ostream& o
) const
{
	o << "(";
	view_j(o, 0);
	o << ")";
}

template <typename elem, int order>
void multicomplex<elem, order>::view_j
(
	std::ostream& o,
	int const sum
) const
{
	real.view_j(o, sum);
	o << " ";
	imag.view_j(o, sum + off);
}

template <typename elem>
void multicomplex<elem, 0>::view_j
(
	std::ostream& o, int const sum
) const
{
	if (real >= 0)
		o << "+ " << +real;
	else
		o << "- " << -real;
	o << "*j" << sum;

	if (imag >= 0)
		o << " + " << +imag;
	else
		o << " - " << -imag;
	o << "*j" << sum + 1;
}

template <typename elem>
void multicomplex<elem, 0>::view_j
(
	std::ostream& o
) const
{
	o << "(";
	view_j(o, 0);
	o << ")";
}

template <typename elem, int order>
void multicomplex<elem, order>::view_k
(
	std::ostream& o
) const
{
	o << "[";
	real.view_k(o);
	o << " ";
	imag.view_k(o);
	o << "]";
}

template <typename elem>
void multicomplex<elem, 0>::view_k
(
	std::ostream& o
) const
{
	o << "[";
	o << real << " ";
	o << imag;
	o << "]";
}

template <typename elem, int order>
std::ostream& operator <<
(
	std::ostream& o,
	multicomplex<elem, order> const& mc
	)
{
	mc.view_j(o);
	return o;
}

template <typename elem>
std::ostream& operator <<
(
	std::ostream& o,
	multicomplex<elem, 0> const& mc
	)
{
	mc.view_i(o);
	return o;
}

//typedefs
	typedef double REAL;
typedef multicomplex<REAL, 0 > MX0;
typedef multicomplex<REAL, 1 > MX1;
typedef multicomplex<REAL, 2 > MX2;
typedef multicomplex<REAL, 3 > MX3;
typedef multicomplex<REAL, 4 > MX4;
typedef multicomplex<REAL, 5 > MX5;
typedef multicomplex<REAL, 6 > MX6;
typedef multicomplex<REAL, 7 > MX7;
typedef multicomplex<REAL, 8 > MX8;
typedef multicomplex<REAL, 9 > MX9;
typedef multicomplex<REAL, 10 > MX10;
typedef multicomplex<REAL, 11 > MX11;
typedef multicomplex<REAL, 12 > MX12;
typedef multicomplex<REAL, 13 > MX13;
typedef multicomplex<REAL, 14 > MX14;
typedef multicomplex<REAL, 15 > MX15;
typedef multicomplex<REAL, 16 > MX16;


constexpr const REAL pi = 3.141592653589793238462643383279502884197169399375105820974;
constexpr const REAL two_pi = 6.2831853071795864769252867665590057683943387987502116419498891846;
constexpr const REAL half_pi = 1.570796326794896619231321691639751442098584699687552910487;
constexpr const REAL euler = 0.577215664901532860606512090082402431042159335939923598805;
constexpr const REAL half = 0.5;
constexpr const REAL quarter = 0.25;
constexpr const REAL root_two_pi = 1.772453850905516027298167483341145182797549456122387128213;

constexpr const REAL lambda = std::numeric_limits<REAL>::epsilon();


//constexpr const REAL lambda = 1e-5;
typedef std::vector<REAL> Vec;

//////////////////
#include <mxws.hpp>
//////////////////
#include <flat.hpp>
//////////////////
#include <operators.hpp>
//////////////////
#include <basic.hpp>
//////////////////
#include <additional_functions.hpp>
//////////////////
#include <derivative.hpp>
//////////////////
#include <root.hpp>
//////////////////
#include <elementry_functions.hpp>
//////////////////
//#include <integrator.hpp>
//////////////////
#include <test.hpp>
//////////////////
#include <linked_list.hpp>
//////////////////
