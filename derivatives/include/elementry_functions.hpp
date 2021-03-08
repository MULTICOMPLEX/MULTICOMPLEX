/// primary functions

template <typename elem, int order>
class multicomplex;

/// sqrt
template <typename elem, int order> 
multicomplex<elem,order> sqrt
(
  const multicomplex<elem,order>& x
) 
{
 // if(x<0 && x > -lambda)return 0;
 // if(x>0 && x <  lambda)return 0;
  
  //multicomplex<elem,order> x_0{half,half};//initial guess
  
  multicomplex<elem,order> x_0;//initial guess 
  
  elem a = 1;
  elem k = 1;
  
  //start in the right quadrant
  if(x.signr()<0){k = -k;}
  if(x.signi()<0){a = -a;}
  
  multicomplex<elem,0> ig{k,a};//initial guess
  x_0 += ig; 
  
  unroll<14>([&](int iter){
    
  x_0 *= (x_0*x_0 + 3*x)/(3*x_0*x_0+x); //Halley 
  //  x_0 = half*(x_0 + (x/x_0)); //Newton 
  });
    
  return x_0;
}

//---------------------------------------------------

template <typename elem, int order> 
inline multicomplex<elem,order> sqrto 
(
  const multicomplex<elem,order>& x
) //Newton method complex
{ 
  multicomplex<elem,order> x_0{half,half};//initial guess
  
  for(int t = 0; t < 8; t++){
    x_0 *= (x_0*x_0 + 3*x)/(3*x_0*x_0+x); //Halley
    //x_0 = half*(x_0 + (x/x_0)); //Newton 
  }
  return x_0;
}

//---------------------------------------------------

template <typename elem, int order> 
multicomplex<elem,order> sqrts 
(
  const multicomplex<elem,order>& x
) //Newton method complex
{ 
  multicomplex<elem,order> x_0{half,half};//initial guess
 
  for(int t = 0; t < 4; t++){
    x_0 *= (x_0*x_0 + 3*x)/(3*x_0*x_0+x); //Halley 
   // x_0 = half*(x_0 + (x/x_0)); //Newton 
  }
  return x_0;
}

//---------------------------------------------------

template <typename elem, int order> 
multicomplex<elem,order> sqrtpow 
(
  const multicomplex<elem,order>& x
) 
{ 
  return pow(x,{half,0});
}

//---------------------------------------------------

/// log
template <typename elem, int order>
multicomplex<elem,order> log 
(
  const multicomplex<elem,order> x
)// log Newton
{   
 // multicomplex<elem,order> x_0(sqrts(x));
 // multicomplex<elem,order> x_0{1,1};
  multicomplex<elem,order> x_0;//initial guess 
  
  elem a = 1;
  elem k = 1;
  
  //start in the right quadrant
  if(x.signr()<0){k = -k;}
  if(x.signi()<0){a = -a;}
  
  multicomplex<elem,0> ig{k,a};//initial guess
  x_0 += ig; 
   
  for(int t = 0; t < 15; t++)
    x_0 -= 2*((exp(x_0)-x)/(exp(x_0)+x));
    
  return x_0;
  //log derivative
  //MX1 ax; 
  //return {ax.real,x.imag/x.real};
}

//---------------------------------------------------

///pow ^complex
template <typename elem, int order>
multicomplex<elem,order> pow 
(
  const multicomplex<elem,order>& b, 
  const multicomplex<elem,order>& e
)
{
	return exp(e * log(b));
}

//---------------------------------------------------

///pow ^REAL
template <typename elem, int order>
multicomplex<elem, order> pow
(
  const multicomplex<elem, order>& b,
  double exp
)
{
  return pow(b, static_cast<multicomplex<elem, order>>(exp));
}

//---------------------------------------------------

// pow ^int
template <typename elem, int order>
inline multicomplex<elem, order> pow
(
  const multicomplex<elem, order>& b, 
  size_t exp
) 
{
  if( exp == 0)
    return 1;
  
  multicomplex<elem,order> temp = pow(b, exp/2);       
  
  if (exp%2 == 0)
      return temp*temp;
  else 
  {
    if(exp > 0)
        return b*temp*temp;
    else
      return (temp*temp)/b; //negative exponent computation 
  }
}

template <typename elem, int order>
inline multicomplex<elem, order> pow
(
  const multicomplex<elem, order>& b,
  int exp
)
{
  if (exp == 0)
    return 1;

  multicomplex<elem, order> temp = pow(b, exp / 2);

  if (exp % 2 == 0)
    return temp * temp;
  else
  {
    if (exp > 0)
      return b * temp * temp;
    else
      return (temp * temp) / b; //negative exponent computation 
  }
}


template <typename elem, int order>
inline multicomplex<elem, order> pow
(
  const multicomplex<elem, order>& b,
  elem exp
)
{
  if (exp == 0)
    return 1;

  multicomplex<elem, order> temp = pow(b, exp / 2);

  if (exp % 2 == 0)
    return temp * temp;
  else
  {
    if (exp > 0)
      return b * temp * temp;
    else
      return (temp * temp) / b; //negative exponent computation 
  }
}
//---------------------------------------------------

/// trigonometric functions

// sine
template <typename elem, int order>
inline multicomplex<elem,order> sin 
(
  multicomplex<elem,order> const & z
) 
{
  return 
  {
    sin (z.real) * cosh (z.imag),
    cos (z.real) * sinh (z.imag)
  };
}

//---------------------------------------------------

template <typename elem>
inline multicomplex<elem, 0> sin 
(
  multicomplex<elem,0> const & z
) 
{
  return 
  {
    std::sin (z.real) * std::cosh (z.imag),
    std::cos (z.real) * std::sinh (z.imag)
  };
}

//---------------------------------------------------

// sine hyperbolic
template <typename elem, int order>
inline multicomplex<elem,order> sinh 
(
  multicomplex<elem,order> const & z) 
{
  return 
  {
    sinh (z.real) * cos (z.imag),
    cosh (z.real) * sin (z.imag)
  };
}

//---------------------------------------------------

template <typename elem>
inline multicomplex<elem, 0> sinh 
(
  multicomplex<elem,0> const & z
)
{
  return 
  {
    std::sinh (z.real) * std::cos (z.imag),
    std::cosh (z.real) * std::sin (z.imag)
  };
}

//---------------------------------------------------

// cosine
template <typename elem, int order>
inline multicomplex<elem,order> cos 
(
  multicomplex<elem,order> const & z
) 
{
 return 
  {
    + cos (z.real) * cosh (z.imag),
    - sin (z.real) * sinh (z.imag)
  };
}

//---------------------------------------------------

template <typename elem>
inline multicomplex<elem, 0> cos 
(
  multicomplex<elem,0> const & z) 
{
  return 
  {
    + std::cos (z.real) * std::cosh (z.imag),
    - std::sin (z.real) * std::sinh (z.imag)
  };
}

//---------------------------------------------------

// tan
template <typename elem, int order>
inline multicomplex<elem,order> tan 
(
  multicomplex<elem,order> const & z
) 
{
  return 
  {
    sin(z) / cos(z)
  };
}

//---------------------------------------------------

// cosine hyperbolic
template <typename elem, int order>
inline multicomplex<elem,order> cosh 
(
  multicomplex<elem,order> const & z
) 
{
  return 
  {
    cosh (z.real) * cos (z.imag),
    sinh (z.real) * sin (z.imag)
  };
}

//---------------------------------------------------

template <typename elem>
inline multicomplex<elem, 0> cosh 
(
  multicomplex<elem,0> const & z
) 
{
  return 
  {
    std::cosh (z.real) * std::cos (z.imag),
    std::sinh (z.real) * std::sin (z.imag)
  };
}

//---------------------------------------------------

// tanh
template <typename elem, int order>
inline multicomplex<elem,order> tanh 
(
  multicomplex<elem,order> const & z
) 
{
  return 
  {
    sinh(z) / cosh(z)
  };
}

//---------------------------------------------------

// exponential
template <typename elem, int order>
inline multicomplex<elem,order> exp 
(
  multicomplex<elem,order> const & z
) 
{
  multicomplex<elem,order-1> const r {exp (z.real)};
  return  
  {
    r * cos (z.imag),
    r * sin (z.imag)
  };
}

//---------------------------------------------------

template <typename elem>
inline multicomplex<elem, 0> exp 
(
  multicomplex<elem,0> const & z
) 
{
  return 
  {
    std::exp (z.real) * std::cos (z.imag),
    std::exp (z.real) * std::sin (z.imag)
  };
}

//---------------------------------------------------

template <typename elem, int order>
inline multicomplex<elem, order> ldexp 
(
  multicomplex<elem,0> const & z,
  int const & e
) 
{
  return 
  {
    z * elem(std::pow(2,e))
  };
}

//---------------------------------------------------

template <typename elem>
inline multicomplex<elem, 0> ldexp 
(
  multicomplex<elem,0> const & z,
  int const & e
) 
{
  return 
  {
    z * elem(std::pow(2,e))
  };
}

//---------------------------------------------------

template <typename elem, int order>
inline multicomplex<elem,order> expl 
(
  multicomplex<elem,order> const & z
) 
{
  multicomplex<elem,order-1> const r {expl (z.real)};
  return  
  {
    r * cos (z.imag),
    r * sin (z.imag)
  };
}

//---------------------------------------------------

template <typename elem>
inline multicomplex<elem, 0> expl 
(
  multicomplex<elem,0> const & z
) 
{
  return 
  {
    elem(expl (z.real)) * std::cos (z.imag),
    elem(expl (z.real)) * std::sin (z.imag)
  };
}

//---------------------------------------------------

template <typename T>
T factorial
(
  std::size_t number
) 
{ 
  T num = T(1); 
  for (size_t i = 1; i <= number; i++) 
      num *= i; 
  return num; 
}

//---------------------------------------------------

const int N = 15;

template <typename elem, int order> 
inline multicomplex<elem,order> Sin 
(
  const multicomplex<elem,order>& x
)
{ 
  multicomplex<elem,order> result;
  elem j;

  unroll<N>([&](size_t n){
 // for(int n = 0; n < N; n++)
 // { 
    j = (1.0/factorial<elem>(2*n+1)) * std::pow(-1, n);
    result += pow(x, 2*n+1) * j;
   // std::cout << result << std::endl;
 // }
  
    });
    
  return result;
}

//---------------------------------------------------

template <typename elem, int order> 
multicomplex<elem,order> Cos 
(
  const multicomplex<elem,order>& x
) 
{ 
  multicomplex<elem,order> result;
  elem j;
  
  for(int n = 0; n < N; n++)
  {
    j = (1.0/factorial<elem>(2*n)) * pow(-1, n);
    result += pow(x, 2*n) * j;
  }

  return result;
}

//---------------------------------------------------

template <typename elem, int order> 
multicomplex<elem,order> Cosh 
(
  const multicomplex<elem,order>& x
)
{ 
  multicomplex<elem,order> result;
  elem j;
  
  for(int n = 0; n < N; n++)
  {
    j = 1.0/factorial<elem>(2*n);
    result += pow(x, 2*n) * j;
  }

  return result;
}

//---------------------------------------------------

template <typename elem, int order> 
multicomplex<elem,order> Sinh 
(
  const multicomplex<elem,order>& x
)
{ 
  multicomplex<elem,order> result;
  elem j;
  
  for(int n = 0; n < N; n++)
  {
    j = 1.0/factorial<elem>(2*n+1);
    result += pow(x, 2*n+1) * j;
  }

  return result;
}

//---------------------------------------------------

template <typename elem, int order> 
multicomplex<elem,order> Exp 
(
  const multicomplex<elem,order>& x
)
{ 
  multicomplex<elem,order> result;
  elem j;  
  elem k=1;
  for(int n = 0; n < 2*N; n++)
  {
    if(n>1) k *= n;
    j = 1.0/k;
    result += pow(x, n) * j;
  }

  return result;
}

//---------------------------------------------------

template <typename elem, int order>
inline const multicomplex<elem, order-1> abs 
(
  const multicomplex<elem, order>& a 
) 
{
  return 
  {
    sqrt(a.real*a.real + a.imag*a.imag)
  };
}

//---------------------------------------------------

template <typename elem>
inline const elem abs 
(
  const multicomplex<elem, 0>& a 
) 
{
  return 
  {
    std::sqrt(a.real*a.real + a.imag*a.imag)
  };
}

//---------------------------------------------------

template<typename elem, int order>
inline multicomplex<elem,order> arg 
(
  const multicomplex<elem,order> & z
) 
{
  std::complex a(z.real.real,z.real.imag);
  auto x = std::arg(a);
  
  multicomplex<elem,order> b{x};
  return 
  {
    b
  };
}

//---------------------------------------------------

template<typename elem>
inline multicomplex<elem,0> arg 
(
  const multicomplex<elem,0> & z
) 
{
  std::complex a(z.real,z.imag);
  auto x = std::arg(a);
  
  multicomplex<elem,0> b{x};
  return 
  {
    b
  };
}

//---------------------------------------------------

template <typename T>
T Fac
(
  std::size_t number
) 
{ 
  T num = T(1); 
  for (size_t i = 1; i <= number; i++) 
      num *= i; 
  return num; 
}

//---------------------------------------------------

template<typename A, typename T>
A Binomial_Coefficient(const T& n, const T& k)
{
  return Fac<A>(n) / (Fac<A>(n-k) * Fac<A>(k));
}

//---------------------------------------------------

template<typename elem, int order>
multicomplex<elem,order> Riemann_Zeta
(
  const multicomplex<elem, order>& s
)
{
  multicomplex<elem, order> sum1 = 0, sum2 = 0;
  
  for(size_t n=0; n <= 34; n++){
  
    sum1 = 0;
    for(size_t k=0; k <= n; k++){
    
    sum1 += elem(std::pow(-1, k)) * 
    Binomial_Coefficient<size_t>(n,k) * pow(multicomplex<elem,order>(elem(k)+1), -s);
  }
  
  sum1 *= 1 / elem(pow(2,n+1));

  sum2 += sum1;
  }
  
  return (1 / (1-pow(multicomplex<elem,order>(2),1-s))) * sum2;
}

//---------------------------------------------------

//RiemannSiegelTheta()
template<typename elem, int order>
multicomplex<elem, order> Riemann_Siegel_theta
(
  const multicomplex<elem, order>& t
)
{
  multicomplex<elem,0> i = {0,1};
  return arg(gamma(quarter + half * i * t.real)) - half * t * std::log(pi);
}

//---------------------------------------------------

template <typename elem>
multicomplex<elem, 0> Polygamma
(
  int n,
  const multicomplex<elem, 0>& z
) 
{   
  mcdv mcdv; 
  
  multicomplex<elem, 0> r;
  
  MX1 x1;
  MX2 x2;
  MX3 x3;
  MX4 x4;
  MX5 x5;
  
  n +=1; 

  if(n == 1){ mcdv.sh<0>(x1, z);  r = mcdv.dv<0>(log(gamma(x1)));}
  if(n == 2){ mcdv.sh<0>(x2, z);  r = mcdv.dv<0>(log(gamma(x2)));} 
  if(n == 3){ mcdv.sh<0>(x3, z);  r = mcdv.dv<0>(log(gamma(x3)));} 
  if(n == 4){ mcdv.sh<0>(x4, z);  r = mcdv.dv<0>(log(gamma(x4)));} 
  if(n == 5){ mcdv.sh<0>(x5, z);  r = mcdv.dv<0>(log(gamma(x5)));}  
  
  return r;
} 

//---------------------------------------------------

template <typename elem>
elem Kronecker_Delta
(
  const elem & k
) 
{   
  if(k==0)return 1;
  return 0;
}

//---------------------------------------------------

template<typename elem>
multicomplex<elem, 0> log_gamma
(
  const multicomplex<elem, 0>& t
)
{
  return log(gamma(t));
}

//---------------------------------------------------

template<typename elem>
multicomplex<elem, 0> Riemann_Siegel_theta
(
  const multicomplex<elem, 0>& t
)
{
  multicomplex<elem,0> i = {0,1};
  
  return -(i/2) * (log_gamma(0.25L+i*t/2) - log_gamma(0.25L-i*t/2)) - ((std::log(pi)*t)/2);
}

//---------------------------------------------------

//RiemannSiegelZ()
template<typename elem, int order>
multicomplex<elem, order> Riemann_Siegel_Z
(
  const multicomplex<elem,order>& t
)
{
  multicomplex<elem,0> i = {0,1};

  return exp(i * Riemann_Siegel_theta(t)) * Riemann_Zeta(half + i * t);
}

//---------------------------------------------------

template<typename elem, int order>
multicomplex<elem, order> Riemann_Siegel_Z2
(
  const multicomplex<elem,order>& t
)
{
  multicomplex<elem,0> i = {0,1};

  return  Riemann_Zeta(half + i * t);
}

//---------------------------------------------------

template <typename elem, int order> 
multicomplex<elem, order> erf
(
  multicomplex<elem, order> const & z
)
{
  multicomplex<elem, order> s = 0;
  for(int n=0; n < 40; n++)
    s += elem(std::pow(-1,n)) * pow(z,2*n+1) / (factorial<elem> (n) * (2*n+1));
  
  return (2./std::sqrt(pi)) * s;
}

//---------------------------------------------------

template <typename T>
T erf
(
  T const & z
)
{
  T s = 0;
  for(int n=0; n < 40; n++)
    s += T(std::pow(-1,n)) * pow(z,2*n+1) / (factorial<T> (n) * (2*n+1));
  
  return (2./std::sqrt(pi)) * s;
}