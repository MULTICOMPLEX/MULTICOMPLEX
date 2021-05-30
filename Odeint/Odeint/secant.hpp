
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


template<typename T, typename F>
class Secant2 {
public:
  Secant2(const T epsilon, const F& f) : mEpsilon(epsilon), mf(f) {}

  Secant2() = default;
  virtual ~Secant2() = default;

  T solve(T x1, T x2) {
    T xm, x0, c;

    do {
      // calculate the intermediate value
      x0 = (x1 * mf(x2) - x2 * mf(x1)) / (mf(x2) - mf(x1));

      // check if x0 is root of equation or not
      c = mf(x1) * mf(x0);

      // update the value of interval
      x1 = x2;
      x2 = x0;

      // update number of iteration
      incrementNumberOfIterations();

      // if x0 is the root of equation then break the loop
      if (c == 0)
        break;
      xm = (x1 * mf(x2) - x2 * mf(x1)) / (mf(x2) - mf(x1));
    } while (fabs(xm - x0) >= mEpsilon); // repeat the loop
                            // until the convergence

    return x0;
  }

private:

  const F& mf;

  int numberOfIterations() const { return mNumberOfIterations; }

  void resetNumberOfIterations() { mNumberOfIterations = 0; }
  int incrementNumberOfIterations() { return mNumberOfIterations++; }
  T epsilon() const { return mEpsilon; }

  const T mEpsilon;

  int mNumberOfIterations = 0;
};