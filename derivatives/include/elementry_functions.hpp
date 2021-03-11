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


template <typename elem, int order> 
multicomplex<elem,order> asin 
(
	const multicomplex<elem,order> & z 
) 
{
	MX0 i;
	i.real = 0;
	i.imag = 1;
		
	return -i * log(sqrt(1 - z*z) + i * z);	//-i log(sqrt(1 - z^2) + i z)	
}

template <typename elem, int order> 
multicomplex<elem,order> acos 
(
	const multicomplex<elem,order> & z 
) 
{
	MX0 i;
	i.real = 0;
	i.imag = 1;
		
	return half_pi + i * log(sqrt(1 - z*z) + i * z);//π/2 + i log(sqrt(1 - z^2) + i z)		
}

template <typename elem, int order> 
multicomplex<elem,order> atan
(
  const multicomplex<elem,order>& z
) 
{
	MX0 i;
	i.real = 0;
	i.imag = 1;
	return half * i * log(1 - i * z) - half * i * log(1 + i * z);//1/2 i log(1 - i z) - 1/2 i log(1 + i z)		
}


template <typename elem, int order> 
multicomplex<elem,order> asinh 
(
	const multicomplex<elem,order> & z 
) 
{		
	return log(sqrt(z*z + 1) + z);	//log(sqrt(z^2 + 1) + z)
}

template <typename elem, int order> 
multicomplex<elem,order> acosh 
(
	const multicomplex<elem,order> & z 
) 
{		
	return log(z + sqrt(z - 1) * sqrt(z + 1));//log(z + sqrt(z - 1) sqrt(z + 1))
}

template <typename elem, int order> 
multicomplex<elem,order> atanh
(
  const multicomplex<elem,order>& z
) 
{
	return half * log(z + 1) - half * log(1 - z);//1/2 log(z + 1) - 1/2 log(1 - z)		
}

template <typename elem, int order> 
multicomplex<elem,order> function
(
  const multicomplex<elem,order>& z0,
  const int formula
) 
{ 
			
			MX0 i;
			i.real = 0;
			i.imag = 1;
			
			MX1 z1;
			MX2 z2;
			MX3 z3;
			MX4 z4;
			
				 if(formula>=19 && formula <=30) sh(z1, z0);
			else if(formula>=31 && formula <=42) sh(z2, z0);
			else if(formula>=43 && formula <=53) sh(z3, z0);
			else if(formula>=54 && formula <=64) sh(z4, z0);
			else if(formula>=65 && formula <=74) sh(z4, z0);
			else if(formula>74) sh(z2, z0);
			
			
			if(formula == 1) return z0;
			
			else if(formula == 2) return log(z0);
			else if(formula == 3) return sqrt(z0);
			else if(formula == 4) return sin(z0);
			else if(formula == 5) return cos(z0);
			else if(formula == 6) return tan(z0);
			else if(formula == 7) return sinh(z0);
			else if(formula == 8) return cosh(z0);
			else if(formula == 9) return tanh(z0);
			else if(formula == 10) return asin(z0);
			else if(formula == 11) return acos(z0);
			else if(formula == 12) return atan(z0);
			else if(formula == 13) return asinh(z0);
			else if(formula == 14) return acosh(z0);
			else if(formula == 15) return atanh(z0);
			else if(formula == 16) return exp(z0);
			else if(formula == 17) return gamma(z0);
			else if(formula == 18) return LambertW(0,z0);
			
			else if(formula == 19) return dv((z1));
			else if(formula == 20) return dv(log(z1));
			else if(formula == 21) return dv(sqrt(z1));
			else if(formula == 22) return dv(sin(z1));
			else if(formula == 23) return dv(cos(z1));
			else if(formula == 24) return dv(tan(z1));
			else if(formula == 25) return dv(sinh(z1));
			else if(formula == 26) return dv(cosh(z1));
			else if(formula == 27) return dv(tanh(z1));
			else if(formula == 28) return dv(exp(z1));
			else if(formula == 29) return dv(gamma(z1));
			else if(formula == 30) return dv(LambertW(0,z1));
			
			else if(formula == 31) return dv((z2));
			else if(formula == 32) return dv(log(z2));
			else if(formula == 33) return dv(sqrt(z2));
			else if(formula == 34) return dv(sin(z2));
			else if(formula == 35) return dv(cos(z2));
			else if(formula == 36) return dv(tan(z2));
			else if(formula == 37) return dv(sinh(z2));
			else if(formula == 38) return dv(cosh(z2));
			else if(formula == 39) return dv(tanh(z2));
			else if(formula == 40) return dv(exp(z2));
			else if(formula == 41) return dv(gamma(z2));
			else if(formula == 42) return dv(LambertW(0,z2));
			
			else if(formula == 43) return dv(log(z3));
			else if(formula == 44) return dv(sqrt(z3));
			else if(formula == 45) return dv(sin(z3));
			else if(formula == 46) return dv(cos(z3));
			else if(formula == 47) return dv(tan(z3));
			else if(formula == 48) return dv(sinh(z3));
			else if(formula == 49) return dv(cosh(z3));
			else if(formula == 50) return dv(tanh(z3));
			else if(formula == 51) return dv(exp(z3));
			else if(formula == 52) return dv(gamma(z3));
			else if(formula == 53) return dv(LambertW(0,z3));
			
			else if(formula == 54) return dv(log(z4));
			else if(formula == 55) return dv(sqrt(z4));
			else if(formula == 56) return dv(sin(z4));
			else if(formula == 57) return dv(cos(z4));
			else if(formula == 58) return dv(tan(z4));
			else if(formula == 59) return dv(sinh(z4));
			else if(formula == 60) return dv(cosh(z4));
			else if(formula == 61) return dv(tanh(z4));
			else if(formula == 62) return dv(exp(z4));
			else if(formula == 63) return dv(gamma(z4));
			else if(formula == 64) return dv(LambertW(0,z4));
			
			else if(formula == 65) return dv(gamma(log(z2))) * i;
			else if(formula == 66) return dv(gamma(sqrt(z2))) * i;
			else if(formula == 67) return dv(gamma(sin(z2))) * i;
			else if(formula == 68) return dv(gamma(cos(z2))) * i;
			else if(formula == 69) return dv(gamma(tan(z2))) * i;
			else if(formula == 70) return dv(gamma(sinh(z2))) * i;
			else if(formula == 71) return dv(gamma(cosh(z2))) * i;
			else if(formula == 72) return dv(gamma(tanh(z2))) * i;
			else if(formula == 73) return dv(gamma(exp(z2))) * i;
			else if(formula == 74) return dv(gamma(LambertW(0,z2))) * i;
			
			
			else return dv(gamma(log(z2)));
		};






