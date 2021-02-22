
#include <random>

class mxws
{

private:

	std::random_device r;

public:

	uint64_t x, w;

	typedef uint32_t result_type;

	void seed()
	{
		init();
	}
	
	mxws(const std::seed_seq& seq)
	{
		if (seq.size() == 2)
		{
			std::vector<uint32_t> seeds(seq.size());
			seq.param(seeds.rbegin());
			w = (uint64_t(seeds[1]) << 32) | seeds[0];
			x = 1;
		}

		else init();
	}

	mxws()
	{
		w = (uint64_t(r()) << 32) | r();
		x = 1;
	}

	mxws(const uint64_t& seed)
	{
		init();
	}

	void init()
	{
		w = (uint64_t(r()) << 32) | r();
		x = 1;
	}

	void init(const uint64_t& seed)
	{
		w = seed;
		x = 1;
	}

	virtual ~mxws() = default;

	static uint32_t min() { return std::numeric_limits<uint32_t>::min(); }
	static uint32_t max() { return std::numeric_limits<uint32_t>::max(); }
	
	inline uint32_t operator()()
	{
		x *= w;
		x = (x >> 32) | (x << 32);
		w += x;
		return uint32_t(x);
	}
	 
	inline double operator()(const double& f)
	{
		return double((*this)()) / 4294967296. * f;
	} 

	inline long double operator()(const long double& min, const long double& max)
	{
		return (*this)(1.) * (max - min) + min;
	}
	
	inline double operator()(const double& min, const double& max)
	{
		return (*this)(1.) * (max - min) + min;
	}

	inline uint32_t operator()(const uint32_t& max)
	{
		return (*this)() % (max + 1);
	}

	inline uint32_t operator()(const uint32_t& min, const uint32_t& max)
	{
		return min + ((*this)() % (max - min + 1));
	}

};


class mxws_64
{

private:

	std::random_device r;

public:

	uint64_t x1, x2, w1, w2;

	typedef uint64_t result_type;

	void seed()
	{
		init();
		x1 = x2 = 1;
	}

	mxws_64(const std::seed_seq& seq)
	{
		if (seq.size() == 2)
		{
			std::vector<uint32_t> seeds(seq.size());
			seq.param(seeds.rbegin());
			w1 = (uint64_t(seeds[1]) << 32) | seeds[0];
			w2 = w1 + 1;
			x1 = x2 = 1;
		}
		 
		else init();
	}

	mxws_64(const uint64_t& seed)
	{
		w1 = (uint64_t(r()) << 32) | r();
		w2 = seed + w1 + 1;
		x1 = x2 = 1;
	}

	void init(const uint64_t& seed)
	{
		w1 = (uint64_t(r()) << 32) | r();
		w2 = seed + w1 + 1;
		x1 = x2 = 1;
	}

	mxws_64()
	{
		init();
	}

	void init()
	{
		w1 = (uint64_t(r()) << 32) | r();
		w2 = w1 + 1;
		x1 = x2 = 1;
	}

	virtual ~mxws_64() = default;

	static uint64_t min() { return std::numeric_limits<uint64_t>::min(); }
	static uint64_t max() { return std::numeric_limits<uint64_t>::max(); }


	inline uint64_t operator()()
	{
		x1 *= w1;
		x1 = (x1 >> 32) | (x1 << 32);
		w1 += x1;

		x2 *= w2;
		x2 = (x2 >> 32) | (x2 << 32);
		w2 += x2;

		return (x1 << 32) | uint32_t(x2);
	}

	inline double operator()(const double& f)
	{
		return double((*this)() >> 11) / 9007199254740992. * f;
	}

	inline double operator()(const double& min, const double& max)
	{
		return (*this)(1.) * (max - min) + min;
	}

	inline long double operator()(const long double& min, const long double& max)
	{
		return (*this)(1.) * (max - min) + min;
	}

	inline uint64_t operator()(const uint64_t& max)
	{
		return (*this)() % (max + 1);
	}
	
	inline uint64_t operator()(const uint64_t& min, const uint64_t& max)
	{
		return min + ((*this)() % (max - min + 1));
	}

};


#ifdef __clang__

template<unsigned T>
class mxws_t_LLVM
{

private:

	std::random_device r;

public:

	static const unsigned n = T * 2;
	typedef unsigned _ExtInt(n) uint;

	uint x, w;

	typedef unsigned _ExtInt(T) result_type;

	void seed()
	{
		init();
		x = 1;
	}

	mxws_t_LLVM(const std::seed_seq& seq)
	{
		if (seq.size() * 32 == n)
		{
			std::vector<uint32_t> seeds(seq.size());
			seq.param(seeds.rbegin());
			for (unsigned x = 0; x < n; x += 32)
				w |= (uint(seeds[x / 32]) << x);
			x = 1;
		}

		else init();
	}

	mxws_t_LLVM(const uint& seed)
	{
		for (unsigned x = 0; x < n; x += 32)
			w |= (uint(r()) << x);
		x = 1;
	}

	void init(const uint64_t& seed)
	{
		for (unsigned x = 0; x < n; x += 32)
			w |= (uint(r()) << x);
		x = 1;
	}

	mxws_t_LLVM()
	{
		init();
	}

	void init()
	{
		for (unsigned x = 0; x < n; x += 32)
			w |= (uint(r()) << x);
		x = 1;
	}

	virtual ~mxws_t_LLVM() = default;

	result_type min() { return std::numeric_limits<result_type>::min(); }
	result_type max() { return std::numeric_limits<result_type>::max(); }

	inline result_type operator()()
	{
		x *= w;
		x = (x >> (n / 2)) | (x << (n / 2));
		w += x;

		return result_type(x);
	}

	inline result_type operator()(const result_type& max)
	{
		return (*this)() % (max + 1);
	}

	inline result_type operator()(const result_type& min, const result_type& max)
	{
		return min + ((*this)() % (max - min + 1));
	}

};

template <unsigned T>
std::ostream& operator <<(std::ostream& os, const unsigned _ExtInt(T)& ExtInt)
{
	os << std::hex;
	for (unsigned i = 0; i < T; i += 4) {
		int b = int(ExtInt >> (T - i - 4)) & 0xf;
		os << b;
	}
	return os;
}

template <unsigned T>
std::ostream& operator <<(std::ostream& os, const _ExtInt(T)& ExtInt)
{
	os << std::hex;
	for (unsigned i = 0; i < T; i += 4) {
		int b = int(ExtInt >> (T - i - 4)) & 0xf;
		os << b;
	}
	return os;
}

template<unsigned T>
inline unsigned _ExtInt(T) rand_t_LLVM()
{
	static thread_local mxws_t_LLVM<T> mxws_t_LLVM;
	return mxws_t_LLVM();
}

#endif
