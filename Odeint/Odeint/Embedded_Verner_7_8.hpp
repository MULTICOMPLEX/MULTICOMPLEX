#pragma once

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Description:                                                              //
//     The Runge-Kutta-Verner method is an adaptive procedure for approxi-    //
//     mating the solution of the differential equation y'(x) = f(x,y) with   //
//     initial condition y(x0) = c.  This implementation evaluates f(x,y)     //
//     thirteen times per step using embedded seventh order and eight order   //
//     Runge-Kutta estimates to estimate the not only the solution but also   //
//     the error.                                                             //
//     The next step size is then calculated using the preassigned tolerance  //
//     and error estimate.                                                    //
//     For step i+1,                                                          //
//        y[i+1] = y[i] +  h * (13/288 * k1 + 32/125 * k6 + 31213/144000 * k7 //
//        + 2401/12375 * k8 + 1701/14080 * k9 + 2401/19200 k10 + 19/450 k11 ) //
//     where                                                                  //
//     k1 = f( x[i],y[i] ),                                                   //
//     k2 = f( x[i]+h/4, y[i] + h*k1/4),                                      //
//     k3 = f( x[i]+h/12, y[i]+h/72*(5 k1 + k2) ),                            //
//     k4 = f( x[i]+h/8, y[i]+h/32*( k1 + 3 k3) ),                            //
//     k5 = f( x[i]+2h/5, y[i]+h/125*(106 k1 - 408 k3 + 352 k4)),             //
//     k6 = f( x[i]+h/2, y[i]+h*( k1/48 + 8 k4 / 33 + 125 k5/528 ) ),         //
//     k7 = f( x[i]+6h/7, y[i]+h/26411*( -13893 k1 + 39936 k4 - 64125 k5      //
//                                                            + 60720 k6 ) ), //
//     k8 = f( x[i]+h/7, y[i]+h*( 37/392 k1 + 1625/9408 k5 - 2/15 k6          //
//                                                           + 61/6720 K7) )  //
//     k9 = f( x[i]+2h/3, y[i]+h*( 17176/25515 k1 - 47104/25515 k4            //
//        + 1325/504 k5 - 41792/25515 k6 + 20237/145800 k7 + 4312/6075 k8) ), //
//     k10 = f( x[i]+2h/7, y[i]+h*( -23834/180075 k1 - 77824/1980825 k4       //
//             - 636635/633864 k5 + 254048/300125 k6 - 183/7000 k7 + 8/11 K8  //
//                                                 - 324/3773 k9) )           //
//     k11 = f( x[i]+h, y[i]+h*( 12733/7600 k1 - 20032/5225 k4                //
//                + 456485/80256 k5 - 42599/7125 k6 + 339227/912000 k7        //
//                          - 1029/4108 K8 + 1701/1408 k9 + 5145/2432 k10) )  //
//     k12 = f( x[i]+h/3, y[i]+h*( -27061/204120 k1 + 40448/280665 k4         //
//                   - 1353775/1197504 k5 + 17662/25515 k6 - 71687/1166400 k7 //
//                        + 98/225 K8 + 1/16 k9 + 3773/11664 k10) )           //
//     k13 = f( x[i]+h, y[i]+h*( 11203/8680 k1 - 38144/11935 k4               //
//             + 2354425/458304 k5 - 84046/16275 k6 + 673309/1636800 k7       //
//             + 4704/8525 K8 + 9477/10912 k9 - 1029/992 k10 + 19/341 k12) )  //
//     x[i+1] = x[i] + h.                                                     //
//                                                                            //
//     The error is estimated to be                                           //
//        err = h*( 6600 k1 + 135168 k6 + 14406 k7 - 57624 k8 - 54675 k9      //
//         + 396165 k10 + 133760 k11 - 437400 k12 - 136400 k13) / 3168000     //
//     The step size h is then scaled by the scale factor                     //
//         scale = 0.8 * | epsilon * y[i] / [err * (xmax - x[0])] | ^ 1/7     //
//     The scale factor is further constrained 0.125 < scale < 4.0.           //
//     The new step size is h := scale * h.                                   //
////////////////////////////////////////////////////////////////////////////////

template<typename F, typename T>
void Runge_Kutta_7_8_v(F f, T t, std::vector<T>& y, T h);

  template<typename F, typename T>
  int Embedded_Verner_7_8(F f, T t, std::vector<T>&y, T h) {

  Runge_Kutta_7_8_v(f, t, y, h);

  return 0;
}

template<typename F, typename T>
void Runge_Kutta_7_8_v(F f, T t, std::vector<T>& y, T h) {

  static const T c1 = 13.0 / 288.0;
  static const T c6 = 32.0 / 125.0;
  static const T c7 = 31213.0 / 144000.0;
  static const T c8 = 2401.0 / 12375.0;
  static const T c9 = 1701.0 / 14080.0;
  static const T c10 = 2401.0 / 19200.0;
  static const T c11 = 19.0 / 450.0;
               
  static const T a2 = 1.0 / 4.0;
  static const T a3 = 1.0 / 12.0;
  static const T a4 = 1.0 / 8.0;
  static const T a5 = 2.0 / 5.0;
  static const T a6 = 1.0 / 2.0;
  static const T a7 = 6.0 / 7.0;
  static const T a8 = 1.0 / 7.0;
  static const T a9 = 2.0 / 3.0;
  static const T a10 = 2.0 / 7.0;
  static const T a12 = 1.0 / 3.0;
               
  static const T b31 = 5.0 / 72.0;
  static const T b32 = 1.0 / 72.0;
  static const T b41 = 1.0 / 32.0;
  static const T b43 = 3.0 / 32.0;
  static const T b51 = 106.0 / 125.0;
  static const T b53 = -408.0 / 125.0;
  static const T b54 = 352.0 / 125.0;
  static const T b61 = 1.0 / 48.0;
  static const T b64 = 8.0 / 33.0;
  static const T b65 = 125.0 / 528.0;
  static const T b71 = -13893.0 / 26411.0;
  static const T b74 = 39936.0 / 26411.0;
  static const T b75 = -64125.0 / 26411.0;
  static const T b76 = 60720.0 / 26411.0;
  static const T b81 = 37.0 / 392.0;
  static const T b85 = 1625.0 / 9408.0;
  static const T b86 = -2.0 / 15.0;
  static const T b87 = 61.0 / 6720.0;
  static const T b91 = 17176.0 / 25515.0;
  static const T b94 = -47104.0 / 25515.0;
  static const T b95 = 1325.0 / 504.0;
  static const T b96 = -41792.0 / 25515.0;
  static const T b97 = 20237.0 / 145800.0;
  static const T b98 = 4312.0 / 6075.0;
  static const T b10_1 = -23834.0 / 180075.0;
  static const T b10_4 = -77824.0 / 1980825.0;
  static const T b10_5 = -636635.0 / 633864.0;
  static const T b10_6 = 254048.0 / 300125.0;
  static const T b10_7 = -183.0 / 7000.0;
  static const T b10_8 = 8.0 / 11.0;
  static const T b10_9 = -324.0 / 3773.0;
  static const T b11_1 = 12733.0 / 7600.0;
  static const T b11_4 = -20032.0 / 5225.0;
  static const T b11_5 = 456485.0 / 80256.0;
  static const T b11_6 = -42599.0 / 7125.0;
  static const T b11_7 = 339227.0 / 912000.0;
  static const T b11_8 = -1029.0 / 4180.0;
  static const T b11_9 = 1701.0 / 1408.0;
  static const T b11_10 = 5145.0 / 2432.0;
  static const T b12_1 = -27061.0 / 204120.0;
  static const T b12_4 = 40448.0 / 280665.0;
  static const T b12_5 = -1353775.0 / 1197504.0;
  static const T b12_6 = 17662.0 / 25515.0;
  static const T b12_7 = -71687.0 / 1166400.0;
  static const T b12_8 = 98.0 / 225.0;
  static const T b12_9 = 1.0 / 16.0;
  static const T b12_10 = 3773.0 / 11664.0;
  static const T b13_1 = 11203.0 / 8680.0;
  static const T b13_4 = -38144.0 / 11935.0;
  static const T b13_5 = 2354425.0 / 458304.0;
  static const T b13_6 = -84046.0 / 16275.0;
  static const T b13_7 = 673309.0 / 1636800.0;
  static const T b13_8 = 4704.0 / 8525.0;
  static const T b13_9 = 9477.0 / 10912.0;
  static const T b13_10 = -1029.0 / 992.0;
  static const T b13_12 = 729.0 / 341.0;
               
  static const T e1 = -6600.0 / 3168000.0;
  static const T e6 = -135168.0 / 3168000.0;
  static const T e7 = -14406.0 / 3168000.0;
  static const T e8 = 57624.0 / 3168000.0;
  static const T e9 = 54675.0 / 3168000.0;
  static const T e10 = -396165.0 / 3168000.0;
  static const T e11 = -133760.0 / 3168000.0;
  static const T e12 = 437400.0 / 3168000.0;
  static const T e13 = 136400.0 / 3168000.0;

  std::vector<T>k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13;

  T h4 = a2 * h;
  T h6 = a3 * h;

  k1 = f(t, y);
  k2 = f(t + h4, y + h4 * k1);
  k3 = f(t + a3 * h, y + h * (b31 * k1 + b32 * k2));
  k4 = f(t + a4 * h, y + h * (b41 * k1 + b43 * k3));
  k5 = f(t + a5 * h, y + h * (b51 * k1 + b53 * k3 + b54 * k4));
  k6 = f(t + a6 * h, y + h * (b61 * k1 + b64 * k4 + b65 * k5));
  k7 = f(t + a7 * h, y + h * (b71 * k1 + b74 * k4 + b75 * k5 + b76 * k6));
  k8 = f(t + a8 * h, y + h * (b81 * k1 + b85 * k5 + b86 * k6 + b87 * k7));
  k9 = f(t + a9 * h, y + h * (b91 * k1 + b94 * k4 + b95 * k5 + b96 * k6
    + b97 * k7 + b98 * k8));
  k10 = f(t + a10 * h, y + h * (b10_1 * k1 + b10_4 * k4 + b10_5 * k5 + b10_6 * k6
    + b10_7 * k7 + b10_8 * k8 + b10_9 * k9));
  k11 = f(t + h, y + h * (b11_1 * k1 + b11_4 * k4 + b11_5 * k5 + b11_6 * k6
    + b11_7 * k7 + b11_8 * k8 + b11_9 * k9 + b11_10 * k10));
  k12 = f(t + a12 * h, y + h * (b12_1 * k1 + b12_4 * k4 + b12_5 * k5 + b12_6 * k6
    + b12_7 * k7 + b12_8 * k8 + b12_9 * k9 + b12_10 * k10));
  k13 = f(t + h, y + h * (b13_1 * k1 + b13_4 * k4 + b13_5 * k5 + b13_6 * k6
    + b13_7 * k7 + b13_8 * k8 + b13_9 * k9 + b13_10 * k10 + b13_12 * k12));
    y += h * (c1 * k1 + c6 * k6 + c7 * k7 + c8 * k8 + c9 * k9
    + c10 * k10 + c11 * k11);
  
  //return e1 * k1 + e6 * k6 + e7 * k7 + e8 * k8 + e9 * k9 + e10 * k10 + e11 * k11
  //  + e12 * k12 + e13 * k13;
}
