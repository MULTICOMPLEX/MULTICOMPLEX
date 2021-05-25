#pragma once

template<typename F, typename T>
void Euler_method(F f, T t, std::vector<T>& y, T h);

template<typename F, typename T>
int Embedded_Euler(F f, T t, std::vector<T>& y, T h) {

  Euler_method(f, t, y, h);

  return 0;
}

template<typename F, typename T>
void Euler_method(F f, T t, std::vector<T>& y, T h)
{
  std::vector<T> k1;

  k1 = f(t, y);
  y += h * k1;

}
