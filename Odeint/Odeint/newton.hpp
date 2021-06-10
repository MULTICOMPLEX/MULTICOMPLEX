
template<typename T, typename F, typename FP>
class Newton {

public:
    Newton(double epsilon, const F &f, const FP &fPrime) : mEpsilon(epsilon), mf(f), mfPrime(fPrime) {}
 
    Newton() = default;
    virtual ~Newton() = default;

    T solve(T x) {
        resetNumberOfIterations();

        T fx = mf(x);
        T fxPrime = mfPrime(x);
        incrementNumberOfIterations();

        while(abs(fx) >= epsilon()) {
            x = calculateX(x, fx, fxPrime);

            fx = mf(x);
            fxPrime = mfPrime(x);

            incrementNumberOfIterations();
        }

        return x;
    }

private:
    T calculateX(T x, T fx, T fxPrime) {
        assert(abs(fxPrime) >= (std::numeric_limits<T>::min)());

        return x - fx/fxPrime;
    }

    const F& mf;
    const FP& mfPrime;

   
    int numberOfIterations() const { return mNumberOfIterations; }

    void resetNumberOfIterations() { mNumberOfIterations = 0; }
    int incrementNumberOfIterations() { return mNumberOfIterations++; }
    T epsilon() const { return mEpsilon; }

    const T mEpsilon;

    int mNumberOfIterations = 0;

};


