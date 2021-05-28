
#include "iteration.h"

template <typename T>
class Brent : public Iteration {
public:
    Brent(T epsilon, const std::function<T (T)> &f) : Iteration(epsilon), mf(f) {}

    Brent() = default;
    virtual ~Brent() = default;
   
    T solve(T a, T b) override {
        resetNumberOfIterations();
        
        auto fa = mf(a);
        auto fb = mf(b);
       
        checkAndFixAlgorithmCriteria(a, b, fa, fb);
      
        //fmt::print("Brent -> [{:}, {:}]\n", a, b);
        //fmt::print("{:<5}|{:<20}|{:<20}|{:<20}|{:<20}|{:<20}\n", "K", "a", "b", "f(a)", "f(b)", "f(s)");
        //fmt::print("---------------------------------------------------------------------------------------------------------- \n");
        //fmt::print("{:<5}|{:<20.15f}|{:<20.15f}|{:<20.15f}|{:<20.15f}\n", incrementNumberOfIterations(), a, b, fa, fb, "");
        incrementNumberOfIterations();
        T lastB = a; // b_{k-1}
        T lastFb = fa;
      

        T s = (std::numeric_limits<double>::max)();
        T fs = (std::numeric_limits<double>::max)();
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

            //fmt::print("{:<5}|{:<20.15f}|{:<20.15f}|{:<20.15f}|{:<20.15f}|{:<20.15f}\n", incrementNumberOfIterations(), a, b, fa, fb, fs);
            incrementNumberOfIterations();
        }

        //fmt::print("\n");

        return fb < fs ? b : s;
    }

private:
  inline T calculateBisection(T a, T b) {
        return 0.5*(a+b);
    }

  inline T calculateSecant(T a, T b, T fa, T fb) {
        //No need to check division by 0, in this case the method returns NAN which is taken care by useSecantMethod method
        return b-fb*(b-a)/(fb-fa);
    }

  inline T calculateInverseQuadraticInterpolation(T a, T b, T lastB, T fa, T fb, T lastFb) {
        return a*fb*lastFb/((fa-fb)*(fa-lastFb)) +
               b*fa*lastFb/((fb-fa)*(fb-lastFb)) +
               lastB*fa*fb/((lastFb-fa)*(lastFb-fb));
    }

    inline bool useInverseQuadraticInterpolation(T fa, T fb, T lastFb) {
        return fa != lastFb && fb != lastFb;
    }

    inline void checkAndFixAlgorithmCriteria(T &a, T &b, T &fa, T &fb) {
        //Algorithm works in range [a,b] if criteria f(a)*f(b) < 0 and f(a) > f(b) is fulfilled
        //assert(fa*fb < 0);
        if (abs(fa) < abs(fb)) {
            std::swap(a, b);
            std::swap(fa, fb);
        }
    }

    inline bool useBisection(bool bisection, T a, T b, T lastB, T penultimateB, T s) const {
        auto DELTA = epsilon() + (std::numeric_limits<T>::min)();

        return (bisection && abs(s-b) >= 0.5*fabs(b-lastB)) ||                 //Bisection was used in last step but |s-b|>=|b-lastB|/2 <- Interpolation step would be to rough, so still use bisection
               (!bisection && abs(s-b) >= 0.5*fabs(lastB-penultimateB)) ||     //Interpolation was used in last step but |s-b|>=|lastB-penultimateB|/2 <- Interpolation step would be to small
               (bisection  && abs(b-lastB) < DELTA) ||                         //If last iteration was using bisection and difference between b and lastB is < delta use bisection for next iteration
               (!bisection && abs(lastB-penultimateB) < DELTA);                //If last iteration was using interpolation but difference between lastB ond penultimateB is < delta use biscetion for next iteration
    }

    const std::function<T (T)> mf;
};

