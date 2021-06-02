#pragma once

template<typename F, typename T>
  class Laplace
  {
  public:

    const F& mf;
    T mEpsilon;

    //LaplaceTransform[sin(t),t,0.5] = 0.8
    Laplace(const T epsilon, const F& f) : mEpsilon(epsilon), mf(f) {
      InitStehfest(DefaultStehfestN);
    }
 
    std::vector<T> V;       //  Stehfest coefficients
    const T ln2 = log(2.0);       //  log of 2
    const int DefaultStehfestN = 24;

    Laplace() = default;
    virtual ~Laplace() = default;

    T Transform(double s)
    {
      const int DefaultIntegralN = 50000;
      T du = 0.5 / (T)DefaultIntegralN;
      T y = -mf(0) / 2.0;
      T u = 0;
      T limit = 1.0 - mEpsilon;
      while (u < limit)
      {
        u += du;
        y += 2.0 * pow(u, s - 1) * mf(-log(u));
        u += du;
        y += pow(u, s - 1) * mf(-log(u));
      }
      return 2.0 * y * du / 3.0;
    }

    //Gaver/Stehfest Algorithm
    T InverseTransform(T t)
    {
      T ln2t = ln2 / t;
      T x = 0;
      T y = 0;
      for (size_t i = 0; i < V.size(); i++)
      {
        x += ln2t;
        y += V[i] * mf(x);
      }
      return ln2t * y;
    }

    T Factorial(int N)
    {
      T x = 1;
      if (N > 1)
      {
        for (int i = 2; i <= N; i++)
          x = i * x;
      }
      return x;
    }

    //public static double Integrate(FunctionDelegate f, double Min, double Max)
    //{
    //    return Integrate(f, Min, Max, 100);
    //}

    //public static double Integrate(FunctionDelegate f, double XMin, double XMax, int N)
    //{
    //    double dx = (XMax - XMin) / (double)N / 2.0;
    //    double y = (f(XMin) - f(XMax))/2.0;
    //    double x = XMin;
    //    double limit = XMax - 1e-10;
    //    while (x < limit)
    //    {
    //        x += dx;
    //        y += 2.0*f(x);
    //        x += dx;
    //        y += f(x);
    //    }
    //    return 2.0 * y * dx / 3.0;
    //}

    void InitStehfest()
    {
      InitStehfest(DefaultStehfestN);
    }

 
    void InitStehfest(int N)
    {
      
      int N2 = N / 2;
      int NV = 2 * N2;
      V.resize(NV);
      int sign = 1;
      if ((N2 % 2) != 0)
        sign = -1;
      for (int i = 0; i < NV; i++)
      {
        int kmin = (i + 2) / 2;
        int kmax = i + 1;
        if (kmax > N2)
          kmax = N2;
        V[i] = 0;
        sign = -sign;
        for (int k = kmin; k <= kmax; k++)
        {
          V[i] = V[i] + (pow(k, N2) / Factorial(k)) * (Factorial(2 * k)
            / Factorial(2 * k - i - 1)) / Factorial(N2 - k) / Factorial(k - 1)
            / Factorial(i + 1 - k);
        }
        V[i] = sign * V[i];
      }
    }
  };

  inline void Laplace_driver()
  {
    auto V = [&](const auto& x)
    {
      return sin(x);
    };

    std::cout.setf(std::ios::fixed, std::ios::floatfield);
    std::cout.precision(15);
    double mEpsilon = 1e-12;
    const auto laplace = new Laplace(mEpsilon, V);

    auto t = laplace->Transform(0.5);
    std::cout << "laplace->Transform(0.5) " << t << std::endl;

    //std::cout << "laplace->InverseTransform(0.5) " << laplace->InverseTransform(t) << std::endl;
  }