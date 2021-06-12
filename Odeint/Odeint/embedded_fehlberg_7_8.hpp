#pragma once


////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Description:                                                              //
//     The Runge-Kutta-Fehlberg method is an adaptive procedure for approxi-  //
//     mating the solution of the differential equation y'(x) = f(x,y) with   //
//     initial condition y(x0) = c.  This implementation evaluates f(x,y)     //
//     thirteen times per step using embedded seventh order and eight order   //
//     Runge-Kutta estimates to estimate the not only the solution but also   //
//     the error.                                                             //
//     The next step size is then calculated using the preassigned tolerance  //
//     and error estimate.                                                    //
//     For step i+1,                                                          //
//        y[i+1] = y[i] +  h * (41/840 * k1 + 34/105 * k6 + 9/35 * k7         //
//                        + 9/35 * k8 + 9/280 * k9 + 9/280 k10 + 41/840 k11 ) //
//     where                                                                  //
//     k1 = f( x[i],y[i] ),                                                   //
//     k2 = f( x[i]+2h/27, y[i] + 2h*k1/27),                                  //
//     k3 = f( x[i]+h/9, y[i]+h/36*( k1 + 3 k2) ),                            //
//     k4 = f( x[i]+h/6, y[i]+h/24*( k1 + 3 k3) ),                            //
//     k5 = f( x[i]+5h/12, y[i]+h/48*(20 k1 - 75 k3 + 75 k4)),                //
//     k6 = f( x[i]+h/2, y[i]+h/20*( k1 + 5 k4 + 4 k5 ) ),                    //
//     k7 = f( x[i]+5h/6, y[i]+h/108*( -25 k1 + 125 k4 - 260 k5 + 250 k6 ) ), //
//     k8 = f( x[i]+h/6, y[i]+h*( 31/300 k1 + 61/225 k5 - 2/9 k6              //
//                                                            + 13/900 K7) )  //
//     k9 = f( x[i]+2h/3, y[i]+h*( 2 k1 - 53/6 k4 + 704/45 k5 - 107/9 k6      //
//                                                      + 67/90 k7 + 3 k8) ), //
//     k10 = f( x[i]+h/3, y[i]+h*( -91/108 k1 + 23/108 k4 - 976/135 k5        //
//                             + 311/54 k6 - 19/60 k7 + 17/6 K8 - 1/12 k9) ), //
//     k11 = f( x[i]+h, y[i]+h*( 2383/4100 k1 - 341/164 k4 + 4496/1025 k5     //
//          - 301/82 k6 + 2133/4100 k7 + 45/82 K8 + 45/164 k9 + 18/41 k10) )  //
//     k12 = f( x[i], y[i]+h*( 3/205 k1 - 6/41 k6 - 3/205 k7 - 3/41 K8        //
//                                                   + 3/41 k9 + 6/41 k10) )  //
//     k13 = f( x[i]+h, y[i]+h*( -1777/4100 k1 - 341/164 k4 + 4496/1025 k5    //
//                      - 289/82 k6 + 2193/4100 k7 + 51/82 K8 + 33/164 k9 +   //
//                                                        12/41 k10 + k12) )  //
//     x[i+1] = x[i] + h.                                                     //
//                                                                            //
//     The error is estimated to be                                           //
//        err = -41/840 * h * ( k1 + k11 - k12 - k13)                         //
//     The step size h is then scaled by the scale factor                     //
//         scale = 0.8 * | epsilon * y[i] / [err * (xmax - x[0])] | ^ 1/7     //
//     The scale factor is further constrained 0.125 < scale < 4.0.           //
//     The new step size is h := scale * h.                                   //
////////////////////////////////////////////////////////////////////////////////


template<typename F, typename T>
void Runge_Kutta_7_8(F& f, const T& t, std::vector<T>& y, const T& h, T reset);

template<typename F, typename T>
void Runge_Kutta_7_8(F& f, const T& t, std::vector<MX0>& y, const T& h, T reset);

template<typename F, typename T>
int Embedded_Fehlberg_7_8(F& f, const T& t, std::vector<MX0>& y, const T& h, T reset) {

  Runge_Kutta_7_8(f, t, y, h, reset);

  return 0;
}

template<typename F, typename T>
int Embedded_Fehlberg_7_8(F& f, const T& t, std::vector<T>& y, const T& h, T reset) {

  Runge_Kutta_7_8(f, t, y, h, reset);

  return 0;
}

template<typename F, typename T>
void Runge_Kutta_7_8(F& f, const T& t, std::vector<T>& y, const T& h, T reset) {

  static const T c_1_11 = 41.0 / 840.0;
  static const T c6 = 34.0 / 105.0;
  static const T c_7_8 = 9.0 / 35.0;
  static const T c_9_10 = 9.0 / 280.0;
               
  static const T a2 = 2.0 / 27.0;
  static const T a3 = 1.0 / 9.0;
  static const T a4 = 1.0 / 6.0;
  static const T a5 = 5.0 / 12.0;
  static const T a6 = 1.0 / 2.0;
  static const T a7 = 5.0 / 6.0;
  static const T a8 = 1.0 / 6.0;
  static const T a9 = 2.0 / 3.0;
  static const T a10 = 1.0 / 3.0;
  
  static const T b31 = 1.0 / 36.0;
  static const T b32 = 3.0 / 36.0;
  static const T b41 = 1.0 / 24.0;
  static const T b43 = 3.0 / 24.0;
  static const T b51 = 20.0 / 48.0;
  static const T b53 = -75.0 / 48.0;
  static const T b54 = 75.0 / 48.0;
  static const T b61 = 1.0 / 20.0;
  static const T b64 = 5.0 / 20.0;
  static const T b65 = 4.0 / 20.0;
  static const T b71 = -25.0 / 108.0;
  static const T b74 = 125.0 / 108.0;
  static const T b75 = -260.0 / 108.0;
  static const T b76 = 250.0 / 108.0;
  static const T b81 = 31.0 / 300.0;
  static const T b85 = 61.0 / 225.0;
  static const T b86 = -2.0 / 9.0;
  static const T b87 = 13.0 / 900.0;
  static const T b91 = 2.0;
  static const T b94 = -53.0 / 6.0;
  static const T b95 = 704.0 / 45.0;
  static const T b96 = -107.0 / 9.0;
  static const T b97 = 67.0 / 90.0;
  static const T b98 = 3.0;
  static const T b10_1 = -91.0 / 108.0;
  static const T b10_4 = 23.0 / 108.0;
  static const T b10_5 = -976.0 / 135.0;
  static const T b10_6 = 311.0 / 54.0;
  static const T b10_7 = -19.0 / 60.0;
  static const T b10_8 = 17.0 / 6.0;
  static const T b10_9 = -1.0 / 12.0;
  static const T b11_1 = 2383.0 / 4100.0;
  static const T b11_4 = -341.0 / 164.0;
  static const T b11_5 = 4496.0 / 1025.0;
  static const T b11_6 = -301.0 / 82.0;
  static const T b11_7 = 2133.0 / 4100.0;
  static const T b11_8 = 45.0 / 82.0;
  static const T b11_9 = 45.0 / 164.0;
  static const T b11_10 = 18.0 / 41.0;
  static const T b12_1 = 3.0 / 205.0;
  static const T b12_6 = -6.0 / 41.0;
  static const T b12_7 = -3.0 / 205.0;
  static const T b12_8 = -3.0 / 41.0;
  static const T b12_9 = 3.0 / 41.0;
  static const T b12_10 = 6.0 / 41.0;
  static const T b13_1 = -1777.0 / 4100.0;
  static const T b13_4 = -341.0 / 164.0;
  static const T b13_5 = 4496.0 / 1025.0;
  static const T b13_6 = -289.0 / 82.0;
  static const T b13_7 = 2193.0 / 4100.0;
  static const T b13_8 = 51.0 / 82.0;
  static const T b13_9 = 33.0 / 164.0;
  static const T b13_10 = 12.0 / 41.0;
               
  static const T err_factor = -41.0 / 840.0;
  
  std::vector<T> k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13;
  
  T h2_7 = a2 * h;

  k1 = f(t, y, reset);
  k2 = f(t + h2_7, y + h2_7 * k1, reset);
  k3 = f(t + a3 * h, y + h * (b31 * k1 + b32 * k2), reset);
  k4 = f(t + a4 * h, y + h * (b41 * k1 + b43 * k3), reset);
  k5 = f(t + a5 * h, y + h * (b51 * k1 + b53 * k3 + b54 * k4), reset);
  k6 = f(t + a6 * h, y + h * (b61 * k1 + b64 * k4 + b65 * k5), reset);
  k7 = f(t + a7 * h, y + h * (b71 * k1 + b74 * k4 + b75 * k5 + b76 * k6), reset);
  k8 = f(t + a8 * h, y + h * (b81 * k1 + b85 * k5 + b86 * k6 + b87 * k7), reset);
  k9 = f(t + a9 * h, y + h * (b91 * k1 + b94 * k4 + b95 * k5 + b96 * k6
    + b97 * k7 + b98 * k8), reset);
  k10 = f(t + a10 * h, y + h * (b10_1 * k1 + b10_4 * k4 + b10_5 * k5 + b10_6 * k6
    + b10_7 * k7 + b10_8 * k8 + b10_9 * k9), reset);
  k11 = f(t + h, y + h * (b11_1 * k1 + b11_4 * k4 + b11_5 * k5 + b11_6 * k6
    + b11_7 * k7 + b11_8 * k8 + b11_9 * k9 + b11_10 * k10), reset);
  k12 = f(t, y + h * (b12_1 * k1 + b12_6 * k6 + b12_7 * k7 + b12_8 * k8
    + b12_9 * k9 + b12_10 * k10), reset);
  k13 = f(t + h, y + h * (b13_1 * k1 + b13_4 * k4 + b13_5 * k5 + b13_6 * k6
    + b13_7 * k7 + b13_8 * k8 + b13_9 * k9 + b13_10 * k10 + k12), reset);
  
  y += h * (c_1_11 * (k1 + k11) + c6 * k6 + c_7_8 * (k7 + k8) + c_9_10 * (k9 + k10));
  //y += h * k1;
 // return err_factor * (k1 + k11 - k12 - k13);
  y *= reset;
}


template<typename F, typename T>
void Runge_Kutta_7_8(F& f, const T& t, std::vector<MX0>& y, const T& h, T reset) {

  static const T c_1_11 = 41.0 / 840.0;
  static const T c6 = 34.0 / 105.0;
  static const T c_7_8 = 9.0 / 35.0;
  static const T c_9_10 = 9.0 / 280.0;

  static const T a2 = 2.0 / 27.0;
  static const T a3 = 1.0 / 9.0;
  static const T a4 = 1.0 / 6.0;
  static const T a5 = 5.0 / 12.0;
  static const T a6 = 1.0 / 2.0;
  static const T a7 = 5.0 / 6.0;
  static const T a8 = 1.0 / 6.0;
  static const T a9 = 2.0 / 3.0;
  static const T a10 = 1.0 / 3.0;

  static const T b31 = 1.0 / 36.0;
  static const T b32 = 3.0 / 36.0;
  static const T b41 = 1.0 / 24.0;
  static const T b43 = 3.0 / 24.0;
  static const T b51 = 20.0 / 48.0;
  static const T b53 = -75.0 / 48.0;
  static const T b54 = 75.0 / 48.0;
  static const T b61 = 1.0 / 20.0;
  static const T b64 = 5.0 / 20.0;
  static const T b65 = 4.0 / 20.0;
  static const T b71 = -25.0 / 108.0;
  static const T b74 = 125.0 / 108.0;
  static const T b75 = -260.0 / 108.0;
  static const T b76 = 250.0 / 108.0;
  static const T b81 = 31.0 / 300.0;
  static const T b85 = 61.0 / 225.0;
  static const T b86 = -2.0 / 9.0;
  static const T b87 = 13.0 / 900.0;
  static const T b91 = 2.0;
  static const T b94 = -53.0 / 6.0;
  static const T b95 = 704.0 / 45.0;
  static const T b96 = -107.0 / 9.0;
  static const T b97 = 67.0 / 90.0;
  static const T b98 = 3.0;
  static const T b10_1 = -91.0 / 108.0;
  static const T b10_4 = 23.0 / 108.0;
  static const T b10_5 = -976.0 / 135.0;
  static const T b10_6 = 311.0 / 54.0;
  static const T b10_7 = -19.0 / 60.0;
  static const T b10_8 = 17.0 / 6.0;
  static const T b10_9 = -1.0 / 12.0;
  static const T b11_1 = 2383.0 / 4100.0;
  static const T b11_4 = -341.0 / 164.0;
  static const T b11_5 = 4496.0 / 1025.0;
  static const T b11_6 = -301.0 / 82.0;
  static const T b11_7 = 2133.0 / 4100.0;
  static const T b11_8 = 45.0 / 82.0;
  static const T b11_9 = 45.0 / 164.0;
  static const T b11_10 = 18.0 / 41.0;
  static const T b12_1 = 3.0 / 205.0;
  static const T b12_6 = -6.0 / 41.0;
  static const T b12_7 = -3.0 / 205.0;
  static const T b12_8 = -3.0 / 41.0;
  static const T b12_9 = 3.0 / 41.0;
  static const T b12_10 = 6.0 / 41.0;
  static const T b13_1 = -1777.0 / 4100.0;
  static const T b13_4 = -341.0 / 164.0;
  static const T b13_5 = 4496.0 / 1025.0;
  static const T b13_6 = -289.0 / 82.0;
  static const T b13_7 = 2193.0 / 4100.0;
  static const T b13_8 = 51.0 / 82.0;
  static const T b13_9 = 33.0 / 164.0;
  static const T b13_10 = 12.0 / 41.0;

  static const T err_factor = -41.0 / 840.0;

  std::vector<MX0> k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13;

  T h2_7 = a2 * h;

  k1 = f(t, y, reset);
  k2 = f(t + h2_7, y + h2_7 * k1, reset);
  k3 = f(t + a3 * h, y + h * (b31 * k1 + b32 * k2), reset);
  k4 = f(t + a4 * h, y + h * (b41 * k1 + b43 * k3), reset);
  k5 = f(t + a5 * h, y + h * (b51 * k1 + b53 * k3 + b54 * k4), reset);
  k6 = f(t + a6 * h, y + h * (b61 * k1 + b64 * k4 + b65 * k5), reset);
  k7 = f(t + a7 * h, y + h * (b71 * k1 + b74 * k4 + b75 * k5 + b76 * k6), reset);
  k8 = f(t + a8 * h, y + h * (b81 * k1 + b85 * k5 + b86 * k6 + b87 * k7), reset);
  k9 = f(t + a9 * h, y + h * (b91 * k1 + b94 * k4 + b95 * k5 + b96 * k6
    + b97 * k7 + b98 * k8), reset);
  k10 = f(t + a10 * h, y + h * (b10_1 * k1 + b10_4 * k4 + b10_5 * k5 + b10_6 * k6
    + b10_7 * k7 + b10_8 * k8 + b10_9 * k9), reset);
  k11 = f(t + h, y + h * (b11_1 * k1 + b11_4 * k4 + b11_5 * k5 + b11_6 * k6
    + b11_7 * k7 + b11_8 * k8 + b11_9 * k9 + b11_10 * k10), reset);
  k12 = f(t, y + h * (b12_1 * k1 + b12_6 * k6 + b12_7 * k7 + b12_8 * k8
    + b12_9 * k9 + b12_10 * k10), reset);
  k13 = f(t + h, y + h * (b13_1 * k1 + b13_4 * k4 + b13_5 * k5 + b13_6 * k6
    + b13_7 * k7 + b13_8 * k8 + b13_9 * k9 + b13_10 * k10 + k12), reset);

  y += h * (c_1_11 * (k1 + k11) + c6 * k6 + c_7_8 * (k7 + k8) + c_9_10 * (k9 + k10));
  //y += h * k1;
 // return err_factor * (k1 + k11 - k12 - k13);
  y *= reset;
}

