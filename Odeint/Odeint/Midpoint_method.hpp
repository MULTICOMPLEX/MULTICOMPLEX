#pragma once

template<typename F, typename T>
void Midpoint_method_implicit(F& f, const T& t, std::vector<T>& y, const T& h, T reset) {

	std::vector<T> k1, k2;
	T k = 0.5 * h;

	k1 = f(t + k, y, reset);

	k2 = f(t + k, 0.5 * (y + (y + h * k1)), reset); // implicit midpoint method

	y += h * k2;

	y *= reset;
}

template<typename F, typename T>
void Midpoint_method_explicit(F& f, const T& t, std::vector<T>& y, const T& h, T reset) {

	std::vector<T> k1, k2;
	T k = 0.5 * h;

	k1 = f(t + k, y, reset);

	k2 = f(t + k, y + k * k1, reset); //explicit midpoint method 

	y += h * k2;

	y *= reset;
}

template<typename F, typename T>
void Midpoint_method_implicit(F& f, const T& t, std::vector<MX0>& y, const T& h, T reset) {

	std::vector<MX0> k1, k2;
	T k = 0.5 * h;

	k1 = f(t + k, y, reset);

	k2 = f(t + k, 0.5 * (y + (y + h * k1)), reset); // implicit midpoint method

	y += h * k2;

	y *= reset;
}

template<typename F, typename T>
void Midpoint_method_explicit(F& f, const T& t, std::vector<MX0>& y, const T& h, T reset) {

	std::vector<MX0> k1, k2;
	T k = 0.5 * h;

	k1 = f(t + k, y, reset);

	k2 = f(t + k, y + 0.5 * h * k1, reset); //explicit midpoint method 

	y += h * k2;

	y *= reset;
}