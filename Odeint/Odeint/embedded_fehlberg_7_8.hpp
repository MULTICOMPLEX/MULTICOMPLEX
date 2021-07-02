
template<typename F, typename T>
void Runge_Kutta_7_8(F& f, const T& t, std::vector<T>& y, const T& h) {

	static const T c_1_11 = 41.0 / 840.0;
	static const T c6 = 34.0 / 105.0;
	static const T c_7_8 = 9.0 / 35.0;
	static const T c_9_10 = 9.0 / 280.0;

	static const T a2 = 2.0 / 27.0;
	static const T a3 = 1.0 / 9.0;
	static const T a4 = 1.0 / 6.0;
	static const T a5 = 5.0 / 12.0;
	static const T a6 = 1.0 / 2.0;
	static const T a7 = 5.0 / 6.0;
	static const T a8 = 1.0 / 6.0;
	static const T a9 = 2.0 / 3.0;
	static const T a10 = 1.0 / 3.0;

	static const T b31 = 1.0 / 36.0;
	static const T b32 = 3.0 / 36.0;
	static const T b41 = 1.0 / 24.0;
	static const T b43 = 3.0 / 24.0;
	static const T b51 = 20.0 / 48.0;
	static const T b53 = -75.0 / 48.0;
	static const T b54 = 75.0 / 48.0;
	static const T b61 = 1.0 / 20.0;
	static const T b64 = 5.0 / 20.0;
	static const T b65 = 4.0 / 20.0;
	static const T b71 = -25.0 / 108.0;
	static const T b74 = 125.0 / 108.0;
	static const T b75 = -260.0 / 108.0;
	static const T b76 = 250.0 / 108.0;
	static const T b81 = 31.0 / 300.0;
	static const T b85 = 61.0 / 225.0;
	static const T b86 = -2.0 / 9.0;
	static const T b87 = 13.0 / 900.0;
	static const T b91 = 2.0;
	static const T b94 = -53.0 / 6.0;
	static const T b95 = 704.0 / 45.0;
	static const T b96 = -107.0 / 9.0;
	static const T b97 = 67.0 / 90.0;
	static const T b98 = 3.0;
	static const T b10_1 = -91.0 / 108.0;
	static const T b10_4 = 23.0 / 108.0;
	static const T b10_5 = -976.0 / 135.0;
	static const T b10_6 = 311.0 / 54.0;
	static const T b10_7 = -19.0 / 60.0;
	static const T b10_8 = 17.0 / 6.0;
	static const T b10_9 = -1.0 / 12.0;
	static const T b11_1 = 2383.0 / 4100.0;
	static const T b11_4 = -341.0 / 164.0;
	static const T b11_5 = 4496.0 / 1025.0;
	static const T b11_6 = -301.0 / 82.0;
	static const T b11_7 = 2133.0 / 4100.0;
	static const T b11_8 = 45.0 / 82.0;
	static const T b11_9 = 45.0 / 164.0;
	static const T b11_10 = 18.0 / 41.0;
	static const T b12_1 = 3.0 / 205.0;
	static const T b12_6 = -6.0 / 41.0;
	static const T b12_7 = -3.0 / 205.0;
	static const T b12_8 = -3.0 / 41.0;
	static const T b12_9 = 3.0 / 41.0;
	static const T b12_10 = 6.0 / 41.0;
	static const T b13_1 = -1777.0 / 4100.0;
	static const T b13_4 = -341.0 / 164.0;
	static const T b13_5 = 4496.0 / 1025.0;
	static const T b13_6 = -289.0 / 82.0;
	static const T b13_7 = 2193.0 / 4100.0;
	static const T b13_8 = 51.0 / 82.0;
	static const T b13_9 = 33.0 / 164.0;
	static const T b13_10 = 12.0 / 41.0;

	static const T err_factor = -41.0 / 840.0;

	std::vector<T> k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13;

	T h2_7 = a2 * h;

	k1 = f(t, y);
	k2 = f(t + h2_7, y + h2_7 * k1);
	k3 = f(t + a3 * h, y + h * (b31 * k1 + b32 * k2));
	k4 = f(t + a4 * h, y + h * (b41 * k1 + b43 * k3));
	k5 = f(t + a5 * h, y + h * (b51 * k1 + b53 * k3 + b54 * k4));
	k6 = f(t + a6 * h, y + h * (b61 * k1 + b64 * k4 + b65 * k5));
	k7 = f(t + a7 * h, y + h * (b71 * k1 + b74 * k4 + b75 * k5 + b76 * k6));
	k8 = f(t + a8 * h, y + h * (b81 * k1 + b85 * k5 + b86 * k6 + b87 * k7));
	k9 = f(t + a9 * h, y + h * (b91 * k1 + b94 * k4 + b95 * k5 + b96 * k6
		+ b97 * k7 + b98 * k8));
	k10 = f(t + a10 * h, y + h * (b10_1 * k1 + b10_4 * k4 + b10_5 * k5 + b10_6 * k6
		+ b10_7 * k7 + b10_8 * k8 + b10_9 * k9));
	k11 = f(t + h, y + h * (b11_1 * k1 + b11_4 * k4 + b11_5 * k5 + b11_6 * k6
		+ b11_7 * k7 + b11_8 * k8 + b11_9 * k9 + b11_10 * k10));
	k12 = f(t, y + h * (b12_1 * k1 + b12_6 * k6 + b12_7 * k7 + b12_8 * k8
		+ b12_9 * k9 + b12_10 * k10));
	k13 = f(t + h, y + h * (b13_1 * k1 + b13_4 * k4 + b13_5 * k5 + b13_6 * k6
		+ b13_7 * k7 + b13_8 * k8 + b13_9 * k9 + b13_10 * k10 + k12));

	y += h * (c_1_11 * (k1 + k11) + c6 * k6 + c_7_8 * (k7 + k8) + c_9_10 * (k9 + k10));
}

template<typename F, typename T>
void Runge_Kutta_7_8(F& f, const T& t, std::vector<MX0>& y, const T& h) {

	static const T c_1_11 = 41.0 / 840.0;
	static const T c6 = 34.0 / 105.0;
	static const T c_7_8 = 9.0 / 35.0;
	static const T c_9_10 = 9.0 / 280.0;

	static const T a2 = 2.0 / 27.0;
	static const T a3 = 1.0 / 9.0;
	static const T a4 = 1.0 / 6.0;
	static const T a5 = 5.0 / 12.0;
	static const T a6 = 1.0 / 2.0;
	static const T a7 = 5.0 / 6.0;
	static const T a8 = 1.0 / 6.0;
	static const T a9 = 2.0 / 3.0;
	static const T a10 = 1.0 / 3.0;

	static const T b31 = 1.0 / 36.0;
	static const T b32 = 3.0 / 36.0;
	static const T b41 = 1.0 / 24.0;
	static const T b43 = 3.0 / 24.0;
	static const T b51 = 20.0 / 48.0;
	static const T b53 = -75.0 / 48.0;
	static const T b54 = 75.0 / 48.0;
	static const T b61 = 1.0 / 20.0;
	static const T b64 = 5.0 / 20.0;
	static const T b65 = 4.0 / 20.0;
	static const T b71 = -25.0 / 108.0;
	static const T b74 = 125.0 / 108.0;
	static const T b75 = -260.0 / 108.0;
	static const T b76 = 250.0 / 108.0;
	static const T b81 = 31.0 / 300.0;
	static const T b85 = 61.0 / 225.0;
	static const T b86 = -2.0 / 9.0;
	static const T b87 = 13.0 / 900.0;
	static const T b91 = 2.0;
	static const T b94 = -53.0 / 6.0;
	static const T b95 = 704.0 / 45.0;
	static const T b96 = -107.0 / 9.0;
	static const T b97 = 67.0 / 90.0;
	static const T b98 = 3.0;
	static const T b10_1 = -91.0 / 108.0;
	static const T b10_4 = 23.0 / 108.0;
	static const T b10_5 = -976.0 / 135.0;
	static const T b10_6 = 311.0 / 54.0;
	static const T b10_7 = -19.0 / 60.0;
	static const T b10_8 = 17.0 / 6.0;
	static const T b10_9 = -1.0 / 12.0;
	static const T b11_1 = 2383.0 / 4100.0;
	static const T b11_4 = -341.0 / 164.0;
	static const T b11_5 = 4496.0 / 1025.0;
	static const T b11_6 = -301.0 / 82.0;
	static const T b11_7 = 2133.0 / 4100.0;
	static const T b11_8 = 45.0 / 82.0;
	static const T b11_9 = 45.0 / 164.0;
	static const T b11_10 = 18.0 / 41.0;
	static const T b12_1 = 3.0 / 205.0;
	static const T b12_6 = -6.0 / 41.0;
	static const T b12_7 = -3.0 / 205.0;
	static const T b12_8 = -3.0 / 41.0;
	static const T b12_9 = 3.0 / 41.0;
	static const T b12_10 = 6.0 / 41.0;
	static const T b13_1 = -1777.0 / 4100.0;
	static const T b13_4 = -341.0 / 164.0;
	static const T b13_5 = 4496.0 / 1025.0;
	static const T b13_6 = -289.0 / 82.0;
	static const T b13_7 = 2193.0 / 4100.0;
	static const T b13_8 = 51.0 / 82.0;
	static const T b13_9 = 33.0 / 164.0;
	static const T b13_10 = 12.0 / 41.0;

	static const T err_factor = -41.0 / 840.0;

	std::vector<MX0> k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13;

	T h2_7 = a2 * h;

	k1 = f(t, y);
	k2 = f(t + h2_7, y + h2_7 * k1);
	k3 = f(t + a3 * h, y + h * (b31 * k1 + b32 * k2));
	k4 = f(t + a4 * h, y + h * (b41 * k1 + b43 * k3));
	k5 = f(t + a5 * h, y + h * (b51 * k1 + b53 * k3 + b54 * k4));
	k6 = f(t + a6 * h, y + h * (b61 * k1 + b64 * k4 + b65 * k5));
	k7 = f(t + a7 * h, y + h * (b71 * k1 + b74 * k4 + b75 * k5 + b76 * k6));
	k8 = f(t + a8 * h, y + h * (b81 * k1 + b85 * k5 + b86 * k6 + b87 * k7));
	k9 = f(t + a9 * h, y + h * (b91 * k1 + b94 * k4 + b95 * k5 + b96 * k6
		+ b97 * k7 + b98 * k8));
	k10 = f(t + a10 * h, y + h * (b10_1 * k1 + b10_4 * k4 + b10_5 * k5 + b10_6 * k6
		+ b10_7 * k7 + b10_8 * k8 + b10_9 * k9));
	k11 = f(t + h, y + h * (b11_1 * k1 + b11_4 * k4 + b11_5 * k5 + b11_6 * k6
		+ b11_7 * k7 + b11_8 * k8 + b11_9 * k9 + b11_10 * k10));
	k12 = f(t, y + h * (b12_1 * k1 + b12_6 * k6 + b12_7 * k7 + b12_8 * k8
		+ b12_9 * k9 + b12_10 * k10));
	k13 = f(t + h, y + h * (b13_1 * k1 + b13_4 * k4 + b13_5 * k5 + b13_6 * k6
		+ b13_7 * k7 + b13_8 * k8 + b13_9 * k9 + b13_10 * k10 + k12));

	y += h * (c_1_11 * (k1 + k11) + c6 * k6 + c_7_8 * (k7 + k8) + c_9_10 * (k9 + k10));
}

template<typename F, typename T>
void Embedded_Fehlberg_7_8(F& f, const T& t, std::vector<MX0>& y, const T& h) {

	Runge_Kutta_7_8(f, t, y, h);

}

template<typename F, typename T>
void Embedded_Fehlberg_7_8(F& f, const T& t, std::vector<T>& y, const T& h) {

	Runge_Kutta_7_8(f, t, y, h);
}
