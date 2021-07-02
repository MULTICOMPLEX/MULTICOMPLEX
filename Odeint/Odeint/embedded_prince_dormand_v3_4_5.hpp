
template<typename F, typename T>
void Runge_Kutta2(const F& f, const T& t, std::vector<T>& y, const T& h) {
	static const T r_9 = 1.0 / 9.0;
	static const T r_2_9 = 2.0 / 9.0;
	static const T r_12 = 1.0 / 12.0;
	static const T r_324 = 1.0 / 324.0;
	static const T r_330 = 1.0 / 330.0;
	static const T r_28 = 1.0 / 28.0;

	std::vector<T> k1, k2, k3, k4, k5, k6, k7;
	T h29 = r_2_9 * h;
	T h9 = r_9 * h;

	k1 = f(t, y);
	k2 = f(t + h29, y + h29 * k1);
	k3 = f(t + 3.0 * h9, y + r_12 * h * (k1 + 3.0 * k2));
	k4 = f(t + 5.0 * h9, y + r_324 * h * (55.0 * k1 - 75.0 * k2 + 200.0 * k3));
	k5 = f(t + 6.0 * h9, y + r_330 * h * (83.0 * k1 - 195.0 * k2
		+ 305.0 * k3 + 27.0 * k4));
	k6 = f(t + h, y + r_28 * h * (-19.0 * k1 + 63.0 * k2
		+ 4.0 * k3 - 108.0 * k4 + 88.0 * k5));
	k7 = f(t + h, y + 0.0025 * h * (38.0 * k1 + 240.0 * k3 - 243.0 * k4
		+ 330.0 * k5 + 35.0 * k6));
	y += h * (0.0862 * k1 + 0.6660 * k3 - 0.7857 * k4
		+ 0.9570 * k5 + 0.0965 * k6 - 0.0200 * k7);
}

template<typename F, typename T>
void Embedded_Prince_Dormand_v3_4_5(const F& f, const T& t, std::vector<T>& y, const T& h) {

	Runge_Kutta2(f, t, y, h);
}