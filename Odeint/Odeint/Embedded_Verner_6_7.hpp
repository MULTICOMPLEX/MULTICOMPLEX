#pragma once

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Description:                                                              //
//     The Runge-Kutta-Verner method is an adaptive procedure for approxi-    //
//     mating the solution of the differential equation y'(x) = f(x,y) with   //
//     initial condition y(x0) = c.  This implementation evaluates f(x,y)     //
//     ten times per step using embedded sixth order and seventh order        //
//     Runge-Kutta estimates to estimate the not only the solution but also   //
//     the error.                                                             //
//     The next step size is then calculated using the preassigned tolerance  //
//     and error estimate.                                                    //
//     For step i+1,                                                          //
//        y[i+1] = y[i] +  h / 90 * ( 7 * k1 + 32 * k4 + 32 * k5              //
//                                                 + 12 * k7 + 7 * k8 )       //
//     where                                                                  //
//     k1 = f( x[i],y[i] ),                                                   //
//     k2 = f( x[i]+h/12, y[i] + h*k1/12),                                    //
//     k3 = f( x[i]+h/6, y[i]+h/6*( k2) ),                                    //
//     k4 = f( x[i]+h/4, y[i]+h/16*( k1 + 3 k3) ),                            //
//     k5 = f( x[i]+3h/4, y[i]+h/16*(21 k1 - 81 k3 + 72 k4)),                 //
//     k6 = f( x[i]+16h/17, y[i]+h/250563*(1344688 k1 - 5127552 k3            //
//                                               + 4096896 k4 - 78208 k5 ) ), //
//     k7 = f( x[i]+h/2, y[i]+h/234624*( -341549 k1 + 1407744 k3 - 1018368 k4 //
//                                                 + 84224 k5 - 14739 k6 ) ), //
//     k8 = f( x[i]+h, y[i]+h/136864*( -381875 k1 + 1642368 k3 - 1327872 k4   //
//                                      + 72192 k5 + 14739 k6 + 117312 K7) )  //
//     k9 = f( x[i]+2h/3, y[i]+h/16755336*( -2070757 k1 + 9929088 k3          //
//                       + 584064 k4 + 3023488 k5 - 447083 k6 + 151424 k7) ), //
//     k10 = f( x[i]+h, y[i]+h/10743824*( 130521209 k1 - 499279872 k3         //
//                                  - 391267968 k4 + 13012608 k5 - 3522621 k6 //
//                                    + 9033024 K7 - 30288492 k9) )           //
//     x[i+1] = x[i] + h.                                                     //
//                                                                            //
//     The error is estimated to be                                           //
//        err = h*( - 1090635 k1 + 9504768 k4 - 171816960 k5 + 72412707 k6    //
//                  - 55840512 k7 - 13412672 k8 + 181730952  k9               //
//                  - 21487648 k10) / 172448640                               //
//     The step size h is then scaled by the scale factor                     //
//         scale = 0.8 * | epsilon * y[i] / [err * (xmax - x[0])] | ^ 1/6     //
//     The scale factor is further constrained 0.125 < scale < 4.0.           //
//     The new step size is h := scale * h.                                   //
////////////////////////////////////////////////////////////////////////////////

template<typename F, typename T>
void Runge_Kutta(F f, T t, Vec& y, T h);


template<typename F, typename T>
int Embedded_Verner_6_7(F f, T t, std::vector<T>& y, T h) {

  Runge_Kutta_6_7(f, t, y, h);

  return 0;
}


template<typename F, typename T>
void Runge_Kutta(F f, T t, Vec& y, T h) {

  static const T c1 = 7.0 / 90.0;
  static const T c4 = 32.0 / 90.0;
  static const T c5 = 32.0 / 90.0;
  static const T c7 = 12.0 / 90.0;
  static const T c8 = 7.0 / 90.0;
               
  static const T a2 = 1.0 / 12.0;
  static const T a3 = 1.0 / 6.0;
  static const T a4 = 1.0 / 4.0;
  static const T a5 = 3.0 / 4.0;
  static const T a6 = 16.0 / 17.0;
  static const T a7 = 1.0 / 2.0;
  static const T a9 = 2.0 / 3.0;
               
  static const T b41 = 1.0 / 16.0;
  static const T b43 = 3.0 / 16.0;
  static const T b51 = 21.0 / 16.0;
  static const T b53 = -81.0 / 16.0;
  static const T b54 = 72.0 / 16.0;
  static const T b61 = 1344688.0 / 250563.0;
  static const T b63 = -5127552.0 / 250563.0;
  static const T b64 = 4096896.0 / 250563.0;
  static const T b65 = -78208.0 / 250563.0;
  static const T b71 = -341549.0 / 234624.0;
  static const T b73 = 1407744.0 / 234624.0;
  static const T b74 = -1018368.0 / 234624.0;
  static const T b75 = 84224.0 / 234624.0;
  static const T b76 = -14739.0 / 234624.0;
  static const T b81 = -381875.0 / 136864.0;
  static const T b83 = 1642368.0 / 136864.0;
  static const T b84 = -1327872.0 / 136864.0;
  static const T b85 = 72192.0 / 136864.0;
  static const T b86 = 14739.0 / 136864.0;
  static const T b87 = 117312.0 / 136864.0;
  static const T b91 = -2070757.0 / 16755336.0;
  static const T b93 = 9929088.0 / 16755336.0;
  static const T b94 = 584064.0 / 16755336.0;
  static const T b95 = 3023488.0 / 16755336.0;
  static const T b96 = -447083.0 / 16755336.0;
  static const T b97 = 151424.0 / 16755336.0;
  static const T b10_1 = 130521209.0 / 10743824.0;
  static const T b10_3 = -499279872.0 / 10743824.0;
  static const T b10_4 = 391267968.0 / 10743824.0;
  static const T b10_5 = 13012608.0 / 10743824.0;
  static const T b10_6 = -3522621.0 / 10743824.0;
  static const T b10_7 = 9033024.0 / 10743824.0;
  static const T b10_9 = -30288492.0 / 10743824.0;

  std::vector<T> k1, k2, k3, k4, k5, k6, k7, k8, k9, k10;

  const T h12 = a2 * h;
  const T h6 = a3 * h;

    k1 = f(t, y);
    k2 = f(t + h12, y + h12 * k1);
    k3 = f(t + h6, y + h6 * k2);
    k4 = f(t + a4 * h, y + h * (b41 * k1 + b43 * k3));
    k5 = f(t + a5 * h, y + h * (b51 * k1 + b53 * k3 + b54 * k4));
    k6 = f(t + a6 * h, y + h * (b61 * k1 + b63 * k3 + b64 * k4 + b65 * k5));
    k7 = f(t + a7 * h, y + h * (b71 * k1 + b73 * k3 + b74 * k4 + b75 * k5 + b76 * k6));
    k8 = f(t + h, y + h * (b81 * k1 + b83 * k3 + b84 * k4 + b85 * k5 + b86 * k6
      + b87 * k7));
    k9 = f(t + a9 * h, y + h * (b91 * k1 + b93 * k3 + b94 * k4 + b95 * k5 + b96 * k6
      + b97 * k7));
    k10 = f(t + h, y + h * (b10_1 * k1 + b10_3 * k3 + b10_4 * k4 + b10_5 * k5
      + b10_6 * k6 + b10_7 * k7 + b10_9 * k9));
   
    y += h * (c1 * k1 + c4 * k4 + c5 * k5 + c7 * k7 + c8 * k8);
}
