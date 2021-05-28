
#include <functional>
//#include <fmt/format.h>

class Iteration {
public:
    inline virtual double solve(double x) { return {}; }
    inline virtual double solve(double a, double b) { return {}; }

    int numberOfIterations() const { return mNumberOfIterations; }
    
    Iteration() = default;
    virtual ~Iteration() = default;

protected:
    explicit Iteration(double epsilon) : mEpsilon(epsilon), mNumberOfIterations(0) {}

    inline void resetNumberOfIterations() { mNumberOfIterations = 0; }
    inline int incrementNumberOfIterations() { return mNumberOfIterations++; }
    double epsilon() const { return mEpsilon; }

private:
    const double mEpsilon;
    int mNumberOfIterations;
};

