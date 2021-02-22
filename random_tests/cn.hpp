#include "stdafx.hpp"

template <typename T>
T standardDeviation
(
  std::vector<T> const & data
)
{
    T sum = 0.0, mean, standardDeviation = 0.0;

    auto n = data.size();
    size_t i;

    for(i = 0; i < n; ++i)
    {
        sum += data[i];
    }

    mean = sum/n;

    for(i = 0; i < n; ++i)
        standardDeviation += pow(data[i] - mean, 2);

    return sqrt(standardDeviation / n);
    //normalized data = (X - mean) / standardDeviation
}

template <typename T>
std::vector<T> Normalize_normal_distribution
(
  std::vector<T> & data
)
{
    T sum = 0.0, mean, standardDeviation = 0.0;

    auto n = data.size();
    size_t i;

    for(i = 0; i < n; ++i)
    {
        sum += data[i];
    }

    mean = sum/n;

    for(i = 0; i < n; ++i)
        standardDeviation += pow(data[i] - mean, 2);

    auto sd = sqrt(standardDeviation / n);
    for(i = 0; i < n; ++i)
      data[i] = (data[i] - mean) / sd;
    
    return data;
}

template <typename T> 
T ShannonEntropy
(
  const std::vector<T>& data
)
{
    T entropy = 0;
    std::map<T, long> counts;

    for (size_t dataIndex = 0; dataIndex < data.size(); ++dataIndex) {
        counts[floor(data[dataIndex] * 50)]++;
        //counts[data[dataIndex]]++;
    }
    
    auto it = counts.begin();
    while(it != counts.end()){
        T p_x = (T)it->second/data.size();
        if (p_x > 0) entropy -= p_x*log(p_x)/log(2);
        it++;
    }
    return entropy;
}

template <typename T> 
T Entropy
(
  const std::vector<T>& data
)
{
    T entropy = 0;
    std::map<T, long> counts;
        
    for (size_t dataIndex = 0; dataIndex < data.size(); ++dataIndex) {
        counts[floor(data[dataIndex] * 50)]++;
        //counts[data[dataIndex]]++;
    }
    
    auto it = counts.begin();
    while(it != counts.end()){
        T p_x = (T)it->second/data.size();
        if (p_x > 0) entropy -= p_x*log(p_x);
        it++;
    }
    return entropy;
}

//https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform
//uniform to normal
template <typename T>
std::vector<T> Box_Muller_transform
(
  const T& U1, 
  const T& U2
)
{   
  std::vector<T> v(2);
  v[0] = sqrt(-2. * log(U1)) * cos(2. * pi * U2);
  v[1] = sqrt(-2. * log(U1)) * sin(2. * pi * U2);
    
  return v; 
}

//https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform
//uniform to normal
template <typename T>
std::vector<T> Box_Muller_transform_uni
(
  const T& mean, 
  const T& stdDev
)
{   
  mxws_64 rng;

  auto u = rng(1.);
  auto v = rng(1.);
  
  std::vector<T> vec(2);
  
  T s = T(stdDev * std::sqrt(-2.0 * std::log(u)));
  
  vec[0] = T(std::cos ( two_pi * v )) * s + mean;
  vec[1] = T(std::sin ( two_pi * v )) * s + mean;

  return vec; 
}

template <typename T>
std::vector<T> Marsaglia_polar
(
  const T& mean, 
  const T& stdDev
) 
{
  mxws_64 rng;

  std::vector<T> vec(2);
  T u, v, s;
    do {
        u = T(rng(-1.,1.));
        v = T(rng(-1.,1.));
        s = u * u + v * v;
    } while (s >= 1 || s == 0);
        
    s = stdDev * T(sqrt(-2.0 * log(s) / s));
      
   vec[0] = u * s + mean;
   vec[1] = v * s + mean;
 
   return vec;
}

//uniform to normal
template <typename elem, int order> 
multicomplex<elem, order> Box_Muller_transform
(
  multicomplex<elem,order> const & z,
  elem const & mu, 
  elem const & sigma
)
{   
  multicomplex<elem,order-1> const r {Box_Muller_transform (z.real, mu, sigma)};
  multicomplex<elem,order-1> const i {Box_Muller_transform (z.imag, mu, sigma)};

  return 
  { 
    r,
    i
  };
}

template <typename elem>
multicomplex<elem, 0> Box_Muller_transform
(
  multicomplex<elem, 0> const& z,
  elem const& mu,
  elem const& sigma
)
{
  elem v1, v2;

  mxws_64 rng;

  v1 = rng(1);
  v2 = rng(1);

  auto mag = sigma * std::sqrt(-2.0 * std::log(v1));

  return
  {
    mag * std::cos(two_pi * v2) + mu,
    mag * std::sin(two_pi * v2) + mu
  };

}

//https://en.wikipedia.org/wiki/Normal_distribution#Cumulative_distribution_function
//normal to uniform
template <typename T>
T Inverse_Box_Muller_transform
(
  const T& x,
  const T& mu,
  const T& sigma
)
{   
  return  (1/2.)*(1 + erf((x-mu)/(sigma * std::sqrt(2)))); 
}

template <typename elem, int order>
multicomplex<elem, order> Inverse_Box_Muller_transform
(
  multicomplex<elem, order> const& x,
  const elem& mu,
  const elem& sigma
)
{   
  return  (1/2.)*(1 + erf((x-mu)/(sigma * std::sqrt(2)))); 
}

template <typename T>
void Inverse_Box_Muller_transform_driver()
{
  std::map<T, T> hist;
  std::queue <int> gquiz;
  
  gquiz.push(10);
  if(!gquiz.empty()){auto p = gquiz.front();std::cout << p << std::endl;}
  if(gquiz.empty()){std::cout << "p" << std::endl;}
  
  T mu = 0.9;
  T sigma = 3.87;
  
  mxws_64 rng;

  std::normal_distribution<> dist(mu, sigma);
  
  MX1 d;
  d = Box_Muller_transform(d, mu, sigma);
   
  //sh(d,{1,2});
  //d = sinh(d);//d/dx(sinh(x)),x=1+2i
  //std::cout << dv(d) << std::endl << std::endl;
  
  std::cout << d << std::endl;
  std::cout << Inverse_Box_Muller_transform(d.real.real, mu, sigma) << std::endl;
  
  std::cout << std::endl;

  MX0 g{1,2};
  g = Box_Muller_transform(g, mu, sigma);
  std::cout << g << std::endl;
  std::cout << Inverse_Box_Muller_transform(g.real, mu, sigma) << std::endl;
  
  std::cout << std::endl;
  
  size_t N = 10000000;
  std::vector<T> data(N);
  
  std::vector<T> vv {0, 1.86, 4.75, 4.75, 1.7554, 1};//1.56071
  std::cout << "Entropy " << Entropy(vv) << std::endl << std::endl;
  
  for (size_t n = 0; n < N; ++n)
  {
    d = Box_Muller_transform(d, mu, sigma);
    ++hist[ floor( Inverse_Box_Muller_transform(d.real.real, mu, sigma) * 20) ];
    
    data[n] = dist(rng);
    
  }
  
  auto k = standardDeviation(data);
  std::cout << "SD " << k << std::endl << std::endl;
  
  std::cout << "Entropy " << Entropy(data) << std::endl << std::endl;
  
  for (auto& p : hist) 
  {
    std::cout << std::internal << std::setw(3) << int(p.first)

    << " : " << std::right << std::setw(6) << std::fixed << std::setprecision(2) << p.second/5000
     
    << " : " << std::string(int(p.second/10000), '*') << std::endl;
  }

  std::cout << std::endl;

}

