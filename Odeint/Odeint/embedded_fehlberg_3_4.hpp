
template<typename F, typename T>
void Runge_Kutta_3_4(const F& f, const T& t, std::vector<T>& y, const T& h)
{
	static const T a2 = 2.0 / 7.0;
	static const T a3 = 7.0 / 15.0;
	static const T a4 = 35.0 / 38.0;

	static const T b31 = 77.0 / 900.0;
	static const T b32 = 343.0 / 900.0;
	static const T b41 = 805.0 / 1444.0;
	static const T b42 = -77175.0 / 54872.0;
	static const T b43 = 97125.0 / 54872.0;
	static const T b51 = 79.0 / 490.0;
	static const T b53 = 2175.0 / 3626.0;
	static const T b54 = 2166.0 / 9065.0;

	static const T c1 = 229.0 / 1470.0;
	static const T c3 = 1125.0 / 1813.0;
	static const T c4 = 13718.0 / 81585.0;
	static const T c5 = 1.0 / 18.0;

	std::vector<T> k1, k2, k3, k4, k5;
	T h2 = a2 * h;

	k1 = f(t, y);
	k2 = f(t + h2, y + h2 * k1);
	k3 = f(t + a3 * h, y + h * (b31 * k1 + b32 * k2));
	k4 = f(t + a4 * h, y + h * (b41 * k1 + b42 * k2 + b43 * k3));
	k5 = f(t + h, y + h * (b51 * k1 + b53 * k3 + b54 * k4));
	y += h * (c1 * k1 + c3 * k3 + c4 * k4 + c5 * k5);
}

template<typename F, typename T>
void Runge_Kutta_3_4(const F& f, const T& t, std::vector<MX0>& y, const T& h)
{
	static const T a2 = 2.0 / 7.0;
	static const T a3 = 7.0 / 15.0;
	static const T a4 = 35.0 / 38.0;

	static const T b31 = 77.0 / 900.0;
	static const T b32 = 343.0 / 900.0;
	static const T b41 = 805.0 / 1444.0;
	static const T b42 = -77175.0 / 54872.0;
	static const T b43 = 97125.0 / 54872.0;
	static const T b51 = 79.0 / 490.0;
	static const T b53 = 2175.0 / 3626.0;
	static const T b54 = 2166.0 / 9065.0;

	static const T c1 = 229.0 / 1470.0;
	static const T c3 = 1125.0 / 1813.0;
	static const T c4 = 13718.0 / 81585.0;
	static const T c5 = 1.0 / 18.0;

	std::vector<MX0> k1, k2, k3, k4, k5;
	T h2 = a2 * h;

	k1 = f(t, y);
	k2 = f(t + h2, y + h2 * k1);
	k3 = f(t + a3 * h, y + h * (b31 * k1 + b32 * k2));
	k4 = f(t + a4 * h, y + h * (b41 * k1 + b42 * k2 + b43 * k3));
	k5 = f(t + h, y + h * (b51 * k1 + b53 * k3 + b54 * k4));
	y += h * (c1 * k1 + c3 * k3 + c4 * k4 + c5 * k5);
}

template<typename F, typename T>
void Embedded_Fehlberg_3_4(const F& f, const T& t, std::vector<T>& y, const T& h) {

	Runge_Kutta_3_4(f, t, y, h);
}

template<typename F, typename T>
void Embedded_Fehlberg_3_4(const F& f, const T& t, std::vector<MX0>& y, const T& h) {

	Runge_Kutta_3_4(f, t, y, h);
}