#pragma once

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Description:                                                              //
//     The Runge-Kutta-Verner method is an adaptive procedure for approxi-    //
//     mating the solution of the differential equation y'(x) = f(x,y) with   //
//     initial condition y(x0) = c.  This implementation evaluates f(x,y)     //
//     eight times per step using embedded fifth order and sixth order        //
//     Runge-Kutta estimates to estimate the not only the solution but also   //
//     the error.                                                             //
//     The next step size is then calculated using the preassigned tolerance  //
//     and error estimate.                                                    //
//     For step i+1,                                                          //
//        y[i+1] = y[i] +  h / 5600 * ( 210 * k1 + 896 * k3 + 1215 * k4       //
//                                                 + 2695 * k5 + 584 * k6 )   //
//     where                                                                  //
//     k1 = f( x[i],y[i] ),                                                   //
//     k2 = f( x[i]+h/18, y[i] + h*k1/18),                                    //
//     k3 = f( x[i]+h/6, y[i]+h/12*(- k1 + 3 k2) ),                           //
//     k4 = f( x[i]+2h/9, y[i]+h/81*(-2 k1 + 12 k2 + 8 k3) ),                 //
//     k5 = f( x[i]+2h/3, y[i]+h/33*(40 k1 - 12 k2 - 168 k3 + 162 k4))        //
//     k6 = f( x[i]+h, y[i]+h/1752*(-8856 k1 + 1728 k2 + 43040 k3             //
//                                                   - 36855 k4 + 2695 k5) )  //
//     k7 = f( x[i]+8h/9, y[i]+h/891*( -8716 k1 + 1968 k2 + 39520 k3          //
//                                                   - 33696 k4 + 1716 k5) )  //
//     k8 = f( x[i]+h, y[i]+h/9984*( 117585 k1 - 22464 k2 - 540032 k3         //
//                                + 466830 k4 - 14014 k5 + 2079 K7) )          //
//     x[i+1] = x[i] + h.                                                     //
//                                                                            //
//     The error is estimated to be                                           //
//        err = h*( 15015 k1 - 118272 k3 + 115830 k4 - 30030 k5 - 30368 k6    //
//                             + 31185 k7 + 16640 k8 ) / 291200               //
//     The step size h is then scaled by the scale factor                     //
//         scale = 0.8 * | epsilon * y[i] / [err * (xmax - x[0])] | ^ 1/5     //
//     The scale factor is further constrained 0.125 < scale < 4.0.           //
//     The new step size is h := scale * h.                                   //
////////////////////////////////////////////////////////////////////////////////

template<typename F, typename T>
void Runge_Kutta_5_6(F f, T t, std::vector<T>& y, T h);

template<typename F, typename T>
int Embedded_Verner_5_6(F f, T t, std::vector<T>& y, T h) 
{

  Runge_Kutta_5_6(f, y, x, h);

  return 0;
}


////////////////////////////////////////////////////////////////////////////////
//  static double Runge_Kutta(double (*f)(double,double), double *y,          //
//                                                       double x0, double h) //
//                                                                            //
//  Description:                                                              //
//     This routine uses Verner's embedded 5th and 6th order methods to       //
//     approximate the solution of the differential equation y'=f(x,y) with   //
//     the initial condition y = y[0] at x = x0.  The value at x + h is       //
//     returned in y[1].  The function returns err / h ( the absolute error   //
//     per step size ).                                                       //
//                                                                            //
//  Arguments:                                                                //
//     double *f  Pointer to the function which returns the slope at (x,y) of //
//                integral curve of the differential equation y' = f(x,y)     //
//                which passes through the point (x0,y[0]).                   //
//     double y[] On input y[0] is the initial value of y at x, on output     //
//                y[1] is the solution at x + h.                              //
//     double x   Initial value of x.                                         //
//     double h   Step size                                                   //
//                                                                            //
//  Return Values:                                                            //
//     This routine returns the err / h.  The solution of y(x) at x + h is    //
//     returned in y[1].                                                      //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

template<typename F, typename T>
void Runge_Kutta_5_6(F f, T t, std::vector<T>& y, T h)
{
  static const T c1 = 210.0 / 5600.0;
  static const T c3 = 896.0 / 5600.0;
  static const T c4 = 1215.0 / 5600.0;
  static const T c5 = 2695.0 / 5600.0;
  static const T c6 = 584.0 / 5600.0;
               
  static const T a2 = 1.0 / 18.0;
  static const T a3 = 1.0 / 6.0;
  static const T a4 = 2.0 / 9.0;
  static const T a5 = 2.0 / 3.0;
  static const T a7 = 8.0 / 9.0;
               
  static const T b31 = -1.0 / 12.0;
  static const T b32 = 3.0 / 12.0;
  static const T b41 = -2.0 / 81.0;
  static const T b42 = 12.0 / 81.0;
  static const T b43 = 8.0 / 81.0;
  static const T b51 = 40.0 / 33.0;
  static const T b52 = -12.0 / 33.0;
  static const T b53 = -168.0 / 33.0;
  static const T b54 = 162.0 / 33.0;
  static const T b61 = -8856.0 / 1752.0;
  static const T b62 = 1728.0 / 1752.0;
  static const T b63 = 43040.0 / 1752.0;
  static const T b64 = -36855.0 / 1752.0;
  static const T b65 = 2695.0 / 1752.0;
  static const T b71 = -8716.0 / 891.0;
  static const T b72 = 1968.0 / 891.0;
  static const T b73 = 39520.0 / 891.0;
  static const T b74 = -33696.0 / 891.0;
  static const T b75 = 1716.0 / 891.0;
  static const T b81 = 117585.0 / 9984.0;
  static const T b82 = -22464.0 / 9984.0;
  static const T b83 = -540032.0 / 9984.0;
  static const T b84 = 466830.0 / 9984.0;
  static const T b85 = -14014.0 / 9984.0;
  static const T b87 = 2079.0 / 9984.0;
               
  static const T e1 = 15015.0 / 291200.0;
  static const T e3 = -118272.0 / 291200.0;
  static const T e4 = 115830.0 / 291200.0;
  static const T e5 = -30030.0 / 291200.0;
  static const T e6 = -30368.0 / 291200.0;
  static const T e7 = 31185.0 / 291200.0;
  static const T e8 = 16640.0 / 291200.0;

  std::vector<T> k1, k2, k3, k4, k5, k6, k7, k8;

  T h18 = a2 * h;

  k1 = f(t, y);
  k2 = f(t + h18, y + h18 * k1);
  k3 = f(t + a3 * h, y + h * (b31 * k1 + b32 * k2));
  k4 = f(t + a4 * h, y + h * (b41 * k1 + b42 * k2 + b43 * k3));
  k5 = f(t + a5 * h, y + h * (b51 * k1 + b52 * k2 + b53 * k3 + b54 * k4));
  k6 = f(t + h, y + h * (b61 * k1 + b62 * k2 + b63 * k3 + b64 * k4 + b65 * k5));
  k7 = f(t + a7 * h, y + h * (b71 * k1 + b72 * k2 + b73 * k3 + b74 * k4 + b75 * k5));
  k8 = f(t + h, y + h * (b81 * k1 + b82 * k2 + b83 * k3 + b84 * k4
    + b85 * k5 + b87 * k7));
  y += y + h * (c1 * k1 + c3 * k3 + c4 * k4 + c5 * k5 + c6 * k6);
  //return e1 * k1 + e3 * k3 + e4 * k4 + e5 * k5 + e6 * k6 + e7 * k7 + e8 * k8;
}
