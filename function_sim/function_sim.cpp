

#include <iostream>
#include <chrono>
#include <iomanip>
#include <string>
#include "mxws.hpp"
#include <numeric>
 
void coprime_driver();

//https://www.jstor.org/stable/2685243?seq=1


template <class T>
double EXP_sim(T& rng, size_t& TEL, double x = 0.9, size_t n_samples = 10000000000)
{
  if (x == 0) return 1;
  
  //mxms rng_32;
  auto e = 2.7182818284590452353602874713526624977572470936999595749669676277;
 
  double h = 0;
  double xi = 1;
  size_t tel = 0;

  if (x < -1 || x > 1)
  {
    x = std::modf(x, &xi);
    //The integer part is stored in the object pointed by intpart, 
    //and the fractional part is returned by the function.
    xi = double(std::pow(e, xi));
  }

  if (x == 0)return xi;

  if (x > 0)
  {
    std::random_device r;

    auto k = 1 / x;
    for (size_t i = 0; i < n_samples; i++)
    {
      while (h < 1)
      {
        h += rng(k);
        //h += rnb(mt);
        tel++;
      }
      h = 0;
    }

    TEL = tel;
    
    return double(tel) / n_samples * xi;
  }
  
  else {
    x = abs(x);
    for (size_t i = 0; i < n_samples; i++)
    {
      while (h < x)
      {
        h += rng(1.);
        tel++;
      }
      h = 0;
    }

    TEL = tel;
    return n_samples / double(tel) * xi;
  }
}

template <class T>
double exp_sim(size_t& TEL, double x, size_t throws, T& rng)
{ 
  std::cout << "x                     = 0x" << std::hex << std::uppercase << rng.x << std::dec << std::endl << std::endl;

  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

  auto r = EXP_sim(rng, TEL, x, throws);

  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

  std::cout << "number of throws      = " << TEL << std::endl << std::endl;
  std::cout << "simulated exp(" << x << ")  = ";
  std::cout << std::fixed << std::setprecision(12) << r << std::endl;
  
  double r2;
  if (std::exp(x) >= r)r2 = exp(x) / r;
  else r2 = r / exp(x);

  std::cout << std::defaultfloat << 
    "exact     exp(" << x << ")  = " << std::setprecision(12) << std::fixed << std::exp(x) << std::endl << std::endl << 
    "difference            = " << abs(1-r2) << std::endl <<
    "central limit         = " << 1./sqrt(TEL) << std::endl << std::endl <<
    "Time difference       = " <<
    std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

  return r;
}

template <class T>
double sqrt_sim(T& rng, size_t& TEL, double z = 2, size_t throws = 10000000000)
{
  size_t tel = 0, i = 0;
  double h;

  std::random_device r;
  
  while (i < throws)
  {
    h = rng(z);
    h *= h;
    if (h < z)tel++;
    i++;
  }

  TEL = i;

  if (z < 1) return (1.0 / z) * double(tel) / throws;
  else return z * double(tel) / throws;
}


int main(int argc, char** argv)
{
  coprime_driver();
  
  //mxws rng;
 
  mxws rng;
  
  size_t throws = 100000000;
  size_t tel = 0;

  if (argc == 1) {
    exp_sim(tel, 9.823475, throws, rng);
    std::cout << std::endl << std::defaultfloat;
    exp_sim(tel, -1.875, throws, rng);
  }

  if (argc == 2)
    exp_sim(tel, std::stof(argv[0]), throws, rng);

  if (argc == 3)
    exp_sim(tel, std::stof(argv[0]), std::stoi(argv[1]), rng);

  std::cout << std::endl << std::endl;

  //sqrt sim
  auto begin = std::chrono::steady_clock::now();
    auto r = sqrt_sim(rng, tel, 8, throws);
  auto end = std::chrono::steady_clock::now();
  double st = 8;

  double r2;
  if (std::sqrt(st) >= r)r2 = sqrt(st) / r;
  else r2 = r / sqrt(st);
  
  std::cout << 
    "number of throws      = " << tel << std::endl << std::endl << 
    "sqrt_sim("<< std::defaultfloat << st <<")           = " << std::setprecision(12) << std::fixed << r << std::endl <<
    "exact sqrt("<< std::defaultfloat << st <<")         = " << std::fixed << std::setprecision(12) << sqrt(st) << std::endl <<
    "difference            = " << abs(1-r2) << std::endl <<
    "central limit         = " << 1. / sqrt(throws) << std::endl << std::endl <<
    "Time difference       = " <<
    std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

  return 0;
}

template<class T>
inline T gcd(T a, T b) {
  while (b) {
    auto t = a % b;
    a = b;
    b = t;
  }
  return a;
}

template<class T>
inline bool are_coprimes(T a, T b) {
  if (!((a | b) & 1))
    return false; // Both are even numbers, divisible by at least 2.
  return 1 == gcd(a, b);
}

void coprime_driver() 
{ 
  mxws_64 rng;
  
  size_t samples = 0;
  size_t successes = 0;

  uint64_t x=0, y=0;

  uint64_t throws = 1000000;

  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

  for (size_t i = 0; i < throws; i++)
  {
    x = rng();
    y = rng();
    samples ++;
    if (are_coprimes(x, y))
      successes ++;
  }

  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

  const double pi = 3.14159265359;
  //const double e =  2.71828182845;

  double d = std::sqrt(6. * ((double)samples) / ((double)successes));

  double r2;
  if (pi >= d)r2 = pi / d;
  else r2 = d / pi;

  std::cout <<
    "coprime sim pi        = " << std::setprecision(12) << std::fixed << d << std::endl <<
    "exact pi              = " << pi << std::endl <<
    "difference            = " << abs(1-r2) << std::endl <<
    "central limit         = " << 1. / sqrt(throws) << std::endl << std::endl << std::defaultfloat <<
    "Time difference       = " <<
    std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl << std::endl;
}