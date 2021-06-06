#pragma once

template<typename F, typename T>
void Euler_method(F f, T t, std::vector<T>& y, T h)
{
  std::vector<T> k1;

  k1 = f(t, y);
  y += h * k1;
}
