#pragma once

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Description:                                                              //
//     The Runge-Kutta-Fehlberg method is an adaptive procedure for approxi-  //
//     mating the solution of the differential equation y'(x) = f(x,y) with   //
//     initial condition y(x0) = c.  This implementation evaluates f(x,y)     //
//     eight times per step using embedded fifth order and sixth order        //
//     Runge-Kutta estimates to estimate the not only the solution but also   //
//     the error.                                                             //
//     The next step size is then calculated using the preassigned tolerance  //
//     and error estimate.                                                    //
//     For step i+1,                                                          //
//        y[i+1] = y[i] +  h / 8448 * ( 682 * k1 + 3375 * k3 + 2376 * k4      //
//                                                 + 1375 * k5 + 640 * k6 )   //
//     where                                                                  //
//     k1 = f( x[i],y[i] ),                                                   //
//     k2 = f( x[i]+h/6, y[i] + h*k1/6 ),                                     //
//     k3 = f( x[i]+4h/15, y[i]+h/75*(4 k1 + 16 k2) ),                        //
//     k4 = f( x[i]+2h/3, y[i]+h/6*(5 k1 - 16 k2 + 15 k3) ),                  //
//     k5 = f( x[i]+4h/5, y[i]+h/25*(-40 k1 + 144 k2 - 100 k3 + 16 k4))       //
//     k6 = f( x[i]+h, y[i]+h/640*(722 k1 - 2304 k2 + 2035 k3                 //
//                                                       - 88 k4 + 275 k5) )  //
//     k7 = f( x[i], y[i]+h/1280*( -22 k1 + 55 k3 - 88 k4 + 55 k5) )          //
//     k8 = f( x[i]+h, y[i]+h/1280*( 186 k1 - 4608 k2 + 4015 k3 - 88 k4       //
//                                             + 495 k5 + 1280 K7) )          //
//     x[i+1] = x[i] + h.                                                     //
//                                                                            //
//     The error is estimated to be                                           //
//        err = - 5*h*( k1 + k6 - k7 - k8 ) / 66                              //
//     The step size h is then scaled by the scale factor                     //
//         scale = 0.8 * | epsilon * y[i] / [err * (xmax - x[0])] | ^ 1/5     //
//     The scale factor is further constrained 0.125 < scale < 4.0.           //
//     The new step size is h := scale * h.                                   //
////////////////////////////////////////////////////////////////////////////////

template<typename F, typename T>
std::vector<T> Runge_Kutta_5_6(F f, T t, std::vector<T>& y, T h, T reset);

////////////////////////////////////////////////////////////////////////////////
// int Embedded_Fehlberg_4_5( double (*f)(double, double), double y[],        //
//       double x, double h, double xmax, double *h_next, double tolerance )  //
//                                                                            //
//  Description:                                                              //
//     This function solves the differential equation y'=f(x,y) with the      //
//     initial condition y(x) = y[0].  The value at xmax is returned in y[1]. //
//     The function returns 0 if successful or -1 if it fails.                //
//                                                                            //
//  Arguments:                                                                //
//     double *f  Pointer to the function which returns the slope at (x,y) of //
//                integral curve of the differential equation y' = f(x,y)     //
//                which passes through the point (x0,y0) corresponding to the //
//                initial condition y(x0) = y0.                               //
//     double y[] On input y[0] is the initial value of y at x, on output     //
//                y[1] is the solution at xmax.                               //
//     double x   The initial value of x.                                     //
//     double h   Initial step size.                                          //
//     double xmax The endpoint of x.                                         //
//     double *h_next   A pointer to the estimated step size for successive   //
//                      calls to Embedded_Fehlberg_5_6.                       //
//     double tolerance The tolerance of y(xmax), i.e. a solution is sought   //
//                so that the relative error < tolerance.                     //
//                                                                            //
//  Return Values:                                                            //
//     0   The solution of y' = f(x,y) from x to xmax is stored y[1] and      //
//         h_next has the value to the next size to try.                      //
//    -1   The solution of y' = f(x,y) from x to xmax failed.                 //
//    -2   Failed because either xmax < x or the step size h <= 0.            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
// 

template<typename F, typename T>
int Embedded_Fehlberg_5_6(F f, T t, std::vector<T>& y, T h, T reset)
{
  Runge_Kutta_5_6(f, t, y, h, reset);

  return 0;
}


////////////////////////////////////////////////////////////////////////////////
//  static double Runge_Kutta(double (*f)(double,double), double *y,          //
//                                                       double x0, double h) //
//                                                                            //
//  Description:                                                              //
//     This routine uses Fehlberg's embedded 5th and 6th order methods to     //
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
std::vector<T> Runge_Kutta_5_6(F f, T t, std::vector<T>& y, T h, T reset)
{
  static const T c1 = 682.0 / 8448.0;
  static const T c3 = 3375.0 / 8448.0;
  static const T c4 = 2376.0 / 8448.0;
  static const T c5 = 1375.0 / 8448.0;
  static const T c6 = 640.0 / 8448.0;
               
  static const T a2 = 1.0 / 6.0;
  static const T a3 = 4.0 / 15.0;
  static const T a4 = 2.0 / 3.0;
  static const T a5 = 0.8;
               
  static const T b31 = 4.0 / 75.0;
  static const T b32 = 16.0 / 75.0;
  static const T b41 = 5.0 / 6.0;
  static const T b42 = -16.0 / 6.0;
  static const T b43 = 15.0 / 6.0;
  static const T b51 = -40.0 / 25.0;
  static const T b52 = 144.0 / 25.0;
  static const T b53 = -100.0 / 25.0;
  static const T b54 = 16.0 / 25.0;
  static const T b61 = 722.0 / 640.0;
  static const T b62 = -2304.0 / 640.0;
  static const T b63 = 2035.0 / 640.0;
  static const T b64 = -88.0 / 640.0;
  static const T b65 = 275.0 / 640.0;
  static const T b71 = -22.0 / 1280.0;
  static const T b73 = 55.0 / 1280.0;
  static const T b74 = -88.0 / 1280.0;
  static const T b75 = 55.0 / 1280.0;
  static const T b81 = 186.0 / 1280.0;
  static const T b82 = -4608.0 / 1280.0;
  static const T b83 = 4015.0 / 1280.0;
  static const T b84 = -88.0 / 1280.0;
  static const T b85 = 495.0 / 1280.0;
               
  static const T err_factor = -5.0 / 66.0;

  std::vector<T>  k1, k2, k3, k4, k5, k6, k7, k8;
  T h6 = a2 * h;

  k1 = f(t, y, reset);
  k2 = f(t + h6, y + h6 * k1, reset);
  k3 = f(t + a3 * h, y + h * (b31 * k1 + b32 * k2), reset);
  k4 = f(t + a4 * h, y + h * (b41 * k1 + b42 * k2 + b43 * k3), reset);
  k5 = f(t + a5 * h, y + h * (b51 * k1 + b52 * k2 + b53 * k3 + b54 * k4), reset);
  k6 = f(t + h, y + h * (b61 * k1 + b62 * k2 + b63 * k3 + b64 * k4 + b65 * k5), reset);
  k7 = f(t, y + h * (b71 * k1 + b73 * k3 + b74 * k4 + b75 * k5), reset);
  k8 = f(t + h, y + h * (b81 * k1 + b82 * k2 + b83 * k3 + b84 * k4
    + b85 * k5 + k7), reset);
  y += h * (c1 * k1 + c3 * k3 + c4 * k4 + c5 * k5 + c6 * k6);
  y *= reset;
  return err_factor * (k1 + k6 - k7 - k8);
}
