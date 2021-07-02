template<typename F, typename T>
void Midpoint_method_implicit(F& f, const T& t, std::vector<T>& y, const T& h) {

	std::vector<T> k1, k2;
	T k = 0.5 * h;

	k1 = f(t + k, y);

	k2 = f(t + k, 0.5 * (y + (y + h * k1))); // implicit midpoint method

	y += h * k2;
}

template<typename T>
void EI1(const T& t, std::vector<T>& y) {

	y *= 2 * std::sin(t);
}

template<typename T>
void EI2(const T& t, std::vector<T>& y) {

	y *= 2 * std::cos(t);
}

template<typename F, typename T>
void Midpoint_method_explicit(F& f, const T& t, std::vector<T>& y, const T& h) {

	std::vector<T> k1, k2;
	T k = 0.5 * h;

	k1 = f(t + k, y);

	k2 = f(t + k, y + k * k1); //explicit midpoint method 

	y += h * k2;
}

template<typename F, typename T, int order>
void Midpoint_method_implicit(F& f, const T& t, std::vector<multicomplex<T, order>>& y, const T& h) {

	std::vector<multicomplex<T, order>> k1, k2;
	T k = 0.5 * h;

	k1 = f(t + k, y);

	k2 = f(t + k, 0.5 * (y + (y + h * k1))); // implicit midpoint method

	y += h * k2;
}

template<typename F, typename T, int order>
void Midpoint_method_explicit(F& f, const T& t, std::vector<multicomplex<T,order>>& y, const T& h) {

	std::vector<multicomplex<T, order>> k1, k2;
	T k = 0.5 * h;

	k1 = f(t + k, y);

	k2 = f(t + k, y + 0.5 * h * k1); //explicit midpoint method 

	y += h * k2;
}

