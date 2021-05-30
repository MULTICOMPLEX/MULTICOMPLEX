
//https://thoughts-on-coding.com/2019/06/06/numerical-methods-with-cpp-part-3-root-approximation-algorithms/
template<typename T, typename F>
class Dekker {
public:
    Dekker(const T epsilon, const F& f) : mEpsilon(epsilon), mf(f) {}
    
    Dekker() = default;
    virtual ~Dekker() = default;

    double solve(T a, T b) {
        resetNumberOfIterations();

        T fa = mf(a);
        T fb = mf(b);

        checkAndFixAlgorithmCriteria(a, b, fa, fb);

        //incrementNumberOfIterations();

        T lastB = a;
        T lastFb = fa;

        while (abs(fb) > epsilon() && abs(b - a) > epsilon()) {
          const T s = calculateSecant(b, fb, lastB, lastFb);
          const T m = calculateBisection(a, b);

          lastB = b;

          b = useSecantMethod(b, s, m) ? s : m;

          lastFb = fb;
          fb = mf(b);

          if (fa * fb > 0 && fb * lastFb < 0) {
            a = lastB;
          }

          fa = mf(a);
          checkAndFixAlgorithmCriteria(a, b, fa, fb);

          incrementNumberOfIterations();
        }

        return b;
    }

private:
    void checkAndFixAlgorithmCriteria(T&a, T&b, T &fa, T &fb) {
        //Algorithm works in range [a,b] if criteria f(a)*f(b) < 0 and f(a) > f(b) is fulfilled
        //assert(fa*fb < 0);
        if (abs(fa) < abs(fb)) {
            std::swap(a, b);
            std::swap(fa, fb);
        }
    }

    T calculateSecant(T b, T fb, T lastB, T lastFb) {
        //No need to check division by 0, in this case the method returns NAN which is taken care by useSecantMethod method
        return b-fb*(b-lastB)/(fb-lastFb);
    }

    T calculateBisection(T a, T b) {
        return 0.5*(a+b);
    }

    bool useSecantMethod(T b, T s, T m) {
        //Value s calculated by secant method has to be between m and b
        return (b > m && s > m && s < b) ||
               (b < m && s > b && s < m);
    }

    const F& mf;

    int numberOfIterations() const { return mNumberOfIterations; }

    void resetNumberOfIterations() { mNumberOfIterations = 0; }
    int incrementNumberOfIterations() { return mNumberOfIterations++; }
    T epsilon() const { return mEpsilon; }

    const T mEpsilon;

    int mNumberOfIterations = 0;
};
