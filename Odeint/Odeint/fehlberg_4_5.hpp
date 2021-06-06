#pragma once

template<typename F, typename T>
void fehlberg_4_5(F f, T t, std::vector<T>& y, T h) {

  std::vector<T> k1, k2, k3, k4, k5, k6;
  const T a1 = 3 / 8;
  const T a2 = 3 / 32;
  const T a3 = 9 / 32;
  const T a4 = 12 / 13;
  const T a5 = 1932 / 2197;
  const T a6 = 7200 / 2197;
  const T a7 = 7296 / 2197;
  const T a8 = 439 / 216;

  const T a9 = 3680 / 513;
  const T a10 = 845 / 4104;
  const T a11 = 8 / 27;
  const T a12 = 3544 / 2565;
  const T a13 = 1859 / 4104;
  const T a14 = 11 / 40;

  const T a15 = 128 / 4275;
  const T a16 = 2197 / 75240;
  const T a17 = 2 / 55;

  k1 = h * f(t, y);
  k2 = h * f(t + h * 0.25, y + k1 * 0.25);
  k3 = h * f(t + h * 3 / 8, y + k1 * 3 / 32 + k2 * 9 / 32);
  k4 = h * f(t + h * 12 / 13, y + k1 * 1932 / 2197 - k2 * 7200 / 2197 + k3 * 7296 / 2197);
  k5 = h * f(t + h, y + k1 * 439 / 216 - k2 * 8 + k3 * 3680 / 513 - k4 * 845 / 4104);
  k6 = h * f(t + h * .5, y - k1 * 8 / 27 + k2 * 2 - k3 * 3544 / 2565 + k4 * 1859 / 4104 - k5 * 11 / 40);

  //R = abs(k1 / 360 - k3 * 128 / 4275 - k4 * 2197 / 75240 + k5 / 50 + k6 * 2 / 55) / h;

  y += k1 * 25 / 216 + k3 * 1408 / 2565 + k4 * 2197 / 4104 - k5 / 5;
}
