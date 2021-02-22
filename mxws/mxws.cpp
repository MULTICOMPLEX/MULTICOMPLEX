
#include <mxws.hpp>
#include <iostream>
#include <chrono>
#include <iomanip>

void eval();
void timing(size_t n = pow(2, 30));

class mxws mxws;
class mxws_64 mxws_64;

int main()
{

  eval();

  mxws.x = 1;
  mxws.w = 1;

  eval();

  timing();

  return 0;
}

void eval()
{
  std::cout << "x = 0x" << mxws.x << std::endl;
  std::cout << "w = 0x" << mxws.w << std::endl;


  for (int t = 0; t < 32; t++)
    std::cout << "0x" << std::setfill('0') << std::setw(8) << std::hex << std::uppercase << mxws() <<
    std::dec << std::endl;

  std::cout << std::endl << std::endl;
}

//https://en.wikipedia.org/wiki/Middle-square_method
//https://crypto.stackexchange.com/questions/62750/attack-of-the-middle-square-weyl-sequence-prng

inline static uint32_t msws()
{ 
  static uint64_t x = 0, w = 0, s = 0xb5ad4eceda1ce2a9;
  x *= x; x += (w += s); return uint32_t(x = (x >> 32) | (x << 32));
}

inline uint64_t msws64()
{
  static uint64_t x1 = 0, w1 = 0, s1 = 0xeeeeeeeeeeeeeeee;
  static uint64_t x2 = 0, w2 = 0, s2 = 0xeeeeeeeeeeeeeeef;
  
  x1 *= x1;
  x1 += (w1 += s1);
  x1 = (x1 >> 32) | (x1 << 32);

  x2 *= x2;
  x2 += (w2 += s2);
  x2 = (x2 >> 32) | (x2 << 32);

  return (x1<<32) | uint32_t(x2);
}

inline std::chrono::steady_clock::time_point now()
{
  return std::chrono::steady_clock::now();
}

template <typename T>
inline size_t duration
(
  const std::chrono::steady_clock::time_point& begin,
  const std::chrono::steady_clock::time_point& end
)
{
  return std::chrono::duration_cast<T>(end - begin).count();
}


void timing(size_t n)
{ 
  std::mt19937 rng1;
  std::mt19937_64 rng2;
  std::uniform_real_distribution<float> dist1;
  std::uniform_real_distribution<double> dist2;
  
  //std::vector<uint32_t> rv(n/8);

  using std::chrono::milliseconds;
  
  auto begin = now();
  for (size_t t = 0; t < n; t++)
    (void)rng1();
  auto end = now();
  std::cout << "\nDuration Mersenne Twister, " << n << " random unsigned 32bit integers = " << std::endl <<
    duration<milliseconds>(begin, end) << "[ms]" << std::endl << std::endl;

  begin = now();
  for (size_t t = 0; t < n; t++)
    (void)dist1(rng1);
  end = now();
  std::cout << "\nDuration Mersenne Twister, " << n << " random floats                  = " << std::endl <<
    duration<milliseconds>(begin, end) << "[ms]" << std::endl << std::endl;

  begin = now();
  for (size_t t = 0; t < n; t++)
    (void)rng2();
  end = now();
  std::cout << "\nDuration Mersenne Twister, " << n << " random unsigned 64bit integers = " << std::endl <<
    duration<milliseconds>(begin, end) << "[ms]" << std::endl << std::endl;

  begin = now();
  for (size_t t = 0; t < n; t++)
    (void)dist2(rng2);
  end = now();
  std::cout << "\nDuration Mersenne Twister, " << n << " random doubles                 = " << std::endl <<
    duration<milliseconds>(begin, end) << "[ms]" << std::endl << std::endl;

  begin = now();
  for (size_t t = 0; t < n; t++)
    msws();
  end = now();
  std::cout << "\nDuration msws(),           " << n << " random unsigned 32bit integers = " << std::endl <<
    duration<milliseconds>(begin, end) << "[ms]" << std::endl << std::endl;

  begin = now();
  for (size_t t = 0; t < n; t++)
    msws64();
  end = now();
  std::cout << "\nDuration msws64(),         " << n << " random unsigned 64bit integers = " << std::endl <<
    duration<milliseconds>(begin, end) << "[ms]" << std::endl << std::endl;

  begin = now();
  for (size_t t = 0; t < n; t++)
    mxws();

  //std::transform(rv.begin(), rv.end(), rv.begin(), [](uint32_t& a) {return a = mxws(); });
  end = now();
  std::cout << "\nDuration mxws(),           " << n << " random unsigned 32bit integers = " << std::endl <<
    duration<milliseconds>(begin, end) << "[ms]" << std::endl << std::endl;

  begin = now();
  for (size_t t = 0; t < n; t++)
    mxws_64();
  end = now();
  std::cout << "\nDuration mxws_64(),        " << n << " random unsigned 64bit integers = " << std::endl <<
    duration<milliseconds>(begin, end) << "[ms]" << std::endl << std::endl;

  begin = now();
  for (size_t t = 0; t < n; t++)
    mxws(1.);
  end = now();
  std::cout << "\nDuration mxws(1.),         " << n << " random floats                  = " << std::endl <<
    duration<milliseconds>(begin, end) << "[ms]" << std::endl << std::endl;

  begin = now();
  for (size_t t = 0; t < n; t++)
    mxws_64(1.);
  end = now();
  std::cout << "\nDuration mxws_64(1.),      " << n << " random doubles                 = " << std::endl <<
    duration<milliseconds>(begin, end) << "[ms]" << std::endl << std::endl;

}