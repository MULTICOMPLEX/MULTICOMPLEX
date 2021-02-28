/// flatten class

#include <vector>

template <typename elem, int order>
class multicomplex;

template <typename elem, int order>
class flat 
{
  
public:

  std::vector<elem> vec;

  const int raw_size {2 << order};

  template <int x_order>
  void help_lower (multicomplex<elem,x_order> const & mc) 
  {
    help_lower<x_order-1> (mc.real);//recursion
    help_lower<x_order-1> (mc.imag);//recursion
  }

  template <int x_order>
  void help_lower (multicomplex<elem,0> const & mc) 
  {
    vec.push_back (mc.real);
    vec.push_back (mc.imag);
  }

  template <int x_order> 
  typename std::enable_if <(x_order > 0), multicomplex<elem,x_order>>::type
  help_raise (size_t const lo, size_t const hi)
  {
    size_t const md {lo + (hi - lo) / 2};

    return 
    {
      help_raise<x_order-1> (lo, md),
      help_raise<x_order-1> (md, hi)
    };
  }

  template <int x_order> 
  typename std::enable_if <(x_order == 0), multicomplex<elem,x_order>>::type
  help_raise (size_t const lo, size_t const hi)
  {
    return {vec.at(lo), vec.at(lo+1)};
  }

  public:
  std::vector<elem> v;
  flat () = delete;
  flat (std::initializer_list<elem> l) : v(l){};
  flat (flat const &) = default;
  flat (flat &&) = default;
  flat & operator= (flat const &) = default;
  flat & operator= (flat &&) = default;
  ~flat () = default;

  flat (multicomplex<elem,order> const & mc) 
  {
    vec.reserve (raw_size);
    help_lower<order> (mc);
  }

  elem & at (int const j)
  { return vec.at(j); }

  int size () const
  { return raw_size; }

  operator multicomplex<elem,order> () 
  { return help_raise<order> (0, raw_size); }

};


// prime for derivative complex

//https://en.wikipedia.org/wiki/Numerical_differentiation




template <typename elem, int order>
multicomplex<elem, order+1> di
(
  const multicomplex<elem, order>& x
)
{
  multicomplex<elem, order+1> z;
  
  flat<elem, order+1> a = z;
  
  a.at(0) = x.real;
  a.at(1) = x.imag;
  
  for(int j=2; j < a.size(); j*=2)
    a.at(j) = lambda;// 2 4 8 16 32 64...
    
  z = a;
  
  return z;
}

template <typename elem, int order>
void sh
(
  multicomplex<elem, order>& z,
  const multicomplex<elem, 0>& x
)
{
  flat<elem, order> a = z;
  
  a.at(0) = x.real;
  a.at(1) = x.imag;
  
  for(int j=2; j < a.size(); j*=2)
    a.at(j) = lambda;// 2 4 8 16 32 64...
    
  z = a;
}

/// retrieve derivative complex

template <typename elem, int order>
multicomplex<elem, 0> dv
(
  const multicomplex<elem, order>& z
)
{ 
  flat<elem, order> a = z;
  auto s = a.size();

  multicomplex<elem, 0> k = {a.at(s-2), a.at(s-1)};
  return k / std::pow(lambda,order);
}

//multicomplex derivatives
class mcdv
{
public :
  
  mcdv()=default;
  virtual ~mcdv()=default;
  
  template <int dv_order, typename elem, int order>  
  void sh
  (
    multicomplex<elem, order>& z,
    const multicomplex<elem, dv_order>& x
  )
  {
     
    flat<elem, order> a = z;
    flat<elem, dv_order> b = x;
  
    int e = int(std::pow(2,dv_order+1));
  
    for(int i = 0; i < e; i++)
      a.at(i) = b.at(i);
  
    for(int j=e; j < a.size(); j*=2)
      a.at(j) = lambda;// 2 4 8 16 32 64...
    
    z = a;
  }
  
  template <int dv_order, typename elem, int order>  
  multicomplex<elem, dv_order> dv
  (
    const multicomplex<elem, order>& z
  )
  { 
    flat<elem, order> a = z;
    auto s = a.size();
  
    int e = int(std::pow(2,dv_order+1));
  
    multicomplex<elem, dv_order> r;
    
    flat<elem, dv_order> k = r;
  
    for(int i = 0; i < e; i++)  
      k.at(i) = a.at(s-i-1);
  
    flat<elem, dv_order> j = k;
    for(int i = 0; i < e; i++)  
      k.at(i) = j.at(e-i-1);
  
    r = k;
    
    return r / std::pow(lambda,order-dv_order);
  }
   
};

template <typename elem, int order>  
multicomplex<elem, order> multicomplex<elem, order>::random(elem lower, elem upper)
{
  class mxws_64 rng;
  flat<elem, order> a = *this;

  for(int j=0; j < a.size(); j++)
    a.at(j) = rng(lower, upper);
    
  *this = a;
  return *this;
}

/// get first 

template <typename elem, int order>
elem gr
(
  const multicomplex<elem, order>& z
)
{
  flat<elem, order> a = z;
  return a.at(0);
}

