
template<typename T, typename F>
class Brent {

public:

    Brent(const T epsilon, const F& f) : mEpsilon(epsilon), mf(f) {}

    Brent() = default;
    virtual ~Brent() = default;
   
    T solve(T a, T b) {
        
        T fa = mf(a);
        T fb = mf(b);
       
        checkAndFixAlgorithmCriteria(a, b, fa, fb);
      
        //incrementNumberOfIterations();
        T lastB = a; // b_{k-1}
        T lastFb = fa;
      
        T s = (std::numeric_limits<T>::max)();
        T fs = (std::numeric_limits<T>::max)();
        T penultimateB = a; // b_{k-2}

        bool bisection = true;
        while(abs(fb) > epsilon() && abs(fs) > epsilon() && abs(b-a) > epsilon()) {
            if(useInverseQuadraticInterpolation(fa, fb, lastFb)) {
                s = calculateInverseQuadraticInterpolation(a, b, lastB, fa, fb, lastFb);
            }
            else {
                s = calculateSecant(a, b, fa, fb);
            }

            if(useBisection(bisection, a, b, lastB, penultimateB, s)) {
                s = calculateBisection(a, b);
                bisection = true;
            }
            else {
                bisection = false;
            }

            fs = mf(s);
            penultimateB = lastB;
            lastB = b;

            if(fa*fs < 0) {
                b = s;
            }
            else {
                a = s;
            }

            fa = mf(a);
            lastFb = fb;
            fb = mf(b);
            checkAndFixAlgorithmCriteria(a, b, fa, fb);

            incrementNumberOfIterations();
        }

        return fb < fs ? b : s;
    }

private:
  T calculateBisection(T a, T b) {
        return 0.5*(a+b);
    }

  T calculateSecant(T a, T b, T fa, T fb) {
        //No need to check division by 0, in this case the method returns NAN which is taken care by useSecantMethod method
        return b-fb*(b-a)/(fb-fa);
    }

  T calculateInverseQuadraticInterpolation(T a, T b, T lastB, T fa, T fb, T lastFb) {
        return a*fb*lastFb/((fa-fb)*(fa-lastFb)) +
               b*fa*lastFb/((fb-fa)*(fb-lastFb)) +
               lastB*fa*fb/((lastFb-fa)*(lastFb-fb));
    }

    inline bool useInverseQuadraticInterpolation(T fa, T fb, T lastFb) {
        return fa != lastFb && fb != lastFb;
    }

    void checkAndFixAlgorithmCriteria(T &a, T &b, T &fa, T &fb) {
        //Algorithm works in range [a,b] if criteria f(a)*f(b) < 0 and f(a) > f(b) is fulfilled
        //assert(fa*fb < 0);
        if (abs(fa) < abs(fb)) {
            std::swap(a, b);
            std::swap(fa, fb);
        }
    }

    bool useBisection(bool bisection, T a, T b, T lastB, T penultimateB, T s) const {
        T DELTA = epsilon() + (std::numeric_limits<T>::min)();

        return (bisection && abs(s-b) >= 0.5*fabs(b-lastB)) ||                 //Bisection was used in last step but |s-b|>=|b-lastB|/2 <- Interpolation step would be to rough, so still use bisection
               (!bisection && abs(s-b) >= 0.5*fabs(lastB-penultimateB)) ||     //Interpolation was used in last step but |s-b|>=|lastB-penultimateB|/2 <- Interpolation step would be to small
               (bisection  && abs(b-lastB) < DELTA) ||                         //If last iteration was using bisection and difference between b and lastB is < delta use bisection for next iteration
               (!bisection && abs(lastB-penultimateB) < DELTA);                //If last iteration was using interpolation but difference between lastB ond penultimateB is < delta use biscetion for next iteration
    }

  const F& mf;

  int numberOfIterations() const { return mNumberOfIterations; }

  void resetNumberOfIterations() { mNumberOfIterations = 0; }
  int incrementNumberOfIterations() { return mNumberOfIterations++; }
  T epsilon() const { return mEpsilon; }

  const T mEpsilon;

  int mNumberOfIterations = 0;

};

