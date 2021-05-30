

template<typename T, typename F>
class Secant {
public:
  Secant(const T epsilon, const F& f) : mEpsilon(epsilon), mf(f) {}

  Secant() = default;
  virtual ~Secant() = default;

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