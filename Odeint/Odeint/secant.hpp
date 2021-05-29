
template<typename T, typename F>
class Secant {
public:
    Secant(const T epsilon, const F& f) : mEpsilon(epsilon), mf(f) {}

    Secant() = default;
    virtual ~Secant() = default;

    T solve(T a, T b) {
        resetNumberOfIterations();

        if(mf(a) > mf(b)) {
            std::swap(a,b);
        }

        T x = b;
        T lastX = a;
        T fx = mf(b);
        T lastFx = mf(a);

        incrementNumberOfIterations();

        while(abs(fx) >= epsilon()) {
            const T x_tmp = calculateX(x, lastX, fx, lastFx);

            lastFx = fx;
            lastX = x;
            x = x_tmp;

            fx = mf(x);

            incrementNumberOfIterations();
        }


        return x;
    }

private:
    T calculateX(T x, T lastX, T fx, T lastFx) {
        const T functionDifference = fx - lastFx;
        assert(abs(functionDifference) >= (std::numeric_limits<T>::min)());

        return x - fx*(x-lastX)/functionDifference;
    }

    const F& mf;

    int numberOfIterations() const { return mNumberOfIterations; }

    void resetNumberOfIterations() { mNumberOfIterations = 0; }
    int incrementNumberOfIterations() { return mNumberOfIterations++; }
    T epsilon() const { return mEpsilon; }

    const T mEpsilon;

    int mNumberOfIterations = 0;
};


