#pragma once

template<typename F, typename T>
void Midpoint_method_implicit(F f, T t, std::vector<T>& y, T h) {

  std::vector<T> k1, k2;

  k1 = f(t + 0.5 * h, y);

  k2 = f(t + 0.5 * h, 0.5 * (y + (y + h * k1))); // implicit midpoint method

  y += h * k2;
}

template<typename F, typename T>
void Midpoint_method_explicit(F f, T t, std::vector<T>& y, T h) {

  std::vector<T> k1, k2;

  k1 = f(t + 0.5 * h, y);

  k2 = f(t + 0.5 * h, y + 0.5 * h * k1); //explicit midpoint method 

  y += h * k2;
}

template<typename F, typename T>
void Midpoint_method_implicit(F f, T t, std::vector<MX0>& y, T h) {

  std::vector<MX0> k1, k2;

  k1 = f(t, y);

  k2 = f(t + 0.5 * h, 0.5 * (y + (y + h * k1))); // implicit midpoint method

  y += h * k2;

}

template<typename F, typename T>
void Midpoint_method_explicit(F f, T t, std::vector<MX0>& y, T h) {

  std::vector<MX0> k1, k2;

  k1 = f(t, y);

  k2 = f(t + 0.5 * h, y + 0.5 * h * k1); //explicit midpoint method 

  y += h * k2;

}