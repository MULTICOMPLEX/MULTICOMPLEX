#pragma once

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Description:                                                              //
//     The Runge-Kutta-Fehlberg method is an adaptive procedure for approxi-  //
//     mating the solution of the differential equation y'(x) = f(x,y) with   //
//     initial condition y(x0) = c.  This implementation evaluates f(x,y) five//
//     times per step using embedded third order and fourth order Runge-Kutta //
//     estimates to estimate the not only the solution but also the error.    //
//     The next step size is then calculated using the preassigned tolerance  //
//     and error estimate.                                                    //
//     For step i+1,                                                          //
//        y[i+1] = y[i] +  h * ( 229 / 1470 * k1 + 1125 / 1813 * k3           //
//                              + 13718 / 81585 * k4 + 1 / 18 * k5 )          //
//     where                                                                  //
//     k1 = f( x[i],y[i] ),                                                   //
//     k2 = f( x[i]+2h/7, y[i] + 2/7h*k1 ),                                   //
//     k3 = f( x[i]+7h/15, y[i]+h*(77/900 k1 + 343/900 k2) ),                 //
//     k4 = f( x[i]+35h/38, y[i]+h*(805/1444 k1 - 77175/54872 k2              //
//                                                     + 97125/54872 k3) ),   //
//     k5 = f( x[i]+h, y[i]+h*(79/490 k1 + 2175/3626 k3 + 2166/9065 k4)),     //
//     x[i+1] = x[i] + h.                                                     //
//                                                                            //
//     The error is estimated to be                                           //
//        err = h * ( - 888 k1 + 3375 k3 - 11552 k4 + 9065 k5 ) / 163170      //
//     The step size h is then scaled by the scale factor                     //
//         scale = 0.8 * | epsilon * y[i] / [err * (xmax - x[0])] | ^ 1/3     //
//     The scale factor is further constrained 0.125 < scale < 4.0.           //
//     The new step size is h := scale * h.                                   //
////////////////////////////////////////////////////////////////////////////////


template<typename F, typename T>
void Runge_Kutta_3_4(F& f, const T& t, std::vector<T>& y, const T& h, T reset);

template<typename F, typename T>
void Runge_Kutta_3_4(F& f, const T& t, std::vector<MX0>& y, const T& h, T reset);

template<typename F, typename T>
int Embedded_Fehlberg_3_4(F& f, const T& t, std::vector<T>& y, const T& h, T reset) {

  Runge_Kutta_3_4(f, t, y, h, reset);

  return 0;
}

template<typename F, typename T>
int Embedded_Fehlberg_3_4(F& f, const T& t, std::vector<MX0>& y, const T& h, T reset) {

  Runge_Kutta_3_4(f, t, y, h, reset);

  return 0;
}

template<typename F, typename T>
void Runge_Kutta_3_4(F& f, const T& t, std::vector<T>& y, const T& h, T reset)
{
  static const T a2 = 2.0 / 7.0;
  static const T a3 = 7.0 / 15.0;
  static const T a4 = 35.0 / 38.0;

  static const T b31 = 77.0 / 900.0;
  static const T b32 = 343.0 / 900.0;
  static const T b41 = 805.0 / 1444.0;
  static const T b42 = -77175.0 / 54872.0;
  static const T b43 = 97125.0 / 54872.0;
  static const T b51 = 79.0 / 490.0;
  static const T b53 = 2175.0 / 3626.0;
  static const T b54 = 2166.0 / 9065.0;

  static const T c1 = 229.0 / 1470.0;
  static const T c3 = 1125.0 / 1813.0;
  static const T c4 = 13718.0 / 81585.0;
  static const T c5 = 1.0 / 18.0;

  static const T d1 = -888.0 / 163170.0;
  static const T d3 = 3375.0 / 163170.0;
  static const T d4 = -11552.0 / 163170.0;
  static const T d5 = 9065.0 / 163170.0;
  
  std::vector<T> k1, k2, k3, k4, k5;
  T h2 = a2 * h;

  k1 = f(t, y, reset);
  k2 = f(t + h2, y + h2 * k1, reset);
  k3 = f(t + a3 * h, y + h * (b31 * k1 + b32 * k2), reset);
  k4 = f(t + a4 * h, y + h * (b41 * k1 + b42 * k2 + b43 * k3), reset);
  k5 = f(t + h, y + h * (b51 * k1 + b53 * k3 + b54 * k4), reset);
  y += h * (c1 * k1 + c3 * k3 + c4 * k4 + c5 * k5);
  //return  d1 * k1 + d3 * k3 + d4 * k4 + d5 * k5;
}

template<typename F, typename T>
void Runge_Kutta_3_4(F& f, const T& t, std::vector<MX0>& y, const T& h, T reset)
{
  static const T a2 = 2.0 / 7.0;
  static const T a3 = 7.0 / 15.0;
  static const T a4 = 35.0 / 38.0;

  static const T b31 = 77.0 / 900.0;
  static const T b32 = 343.0 / 900.0;
  static const T b41 = 805.0 / 1444.0;
  static const T b42 = -77175.0 / 54872.0;
  static const T b43 = 97125.0 / 54872.0;
  static const T b51 = 79.0 / 490.0;
  static const T b53 = 2175.0 / 3626.0;
  static const T b54 = 2166.0 / 9065.0;

  static const T c1 = 229.0 / 1470.0;
  static const T c3 = 1125.0 / 1813.0;
  static const T c4 = 13718.0 / 81585.0;
  static const T c5 = 1.0 / 18.0;

  static const T d1 = -888.0 / 163170.0;
  static const T d3 = 3375.0 / 163170.0;
  static const T d4 = -11552.0 / 163170.0;
  static const T d5 = 9065.0 / 163170.0;

  std::vector<MX0> k1, k2, k3, k4, k5;
  T h2 = a2 * h;

  k1 = f(t, y, reset);
  k2 = f(t + h2, y + h2 * k1, reset);
  k3 = f(t + a3 * h, y + h * (b31 * k1 + b32 * k2), reset);
  k4 = f(t + a4 * h, y + h * (b41 * k1 + b42 * k2 + b43 * k3), reset);
  k5 = f(t + h, y + h * (b51 * k1 + b53 * k3 + b54 * k4), reset);
  y += h * (c1 * k1 + c3 * k3 + c4 * k4 + c5 * k5);
  y *= reset;
  //return  d1 * k1 + d3 * k3 + d4 * k4 + d5 * k5;
}
