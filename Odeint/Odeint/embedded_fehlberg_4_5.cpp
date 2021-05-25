
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Description:                                                              //
//     The Runge-Kutta-Fehlberg method is an adaptive procedure for approxi-  //
//     mating the solution of the differential equation y'(x) = f(x,y) with   //
//     initial condition y(x0) = c.  This implementation evaluates f(x,y) six //
//     times per step using embedded fourth order and fifth order Runge-Kutta //
//     estimates to estimate the not only the solution but also the error.    //
//     The next step size is then calculated using the preassigned tolerance  //
//     and error estimate.                                                    //
//     For step i+1,                                                          //
//        y[i+1] = y[i] +  h * ( 25 / 216 * k1 + 1408 / 2565 * k3             //
//                              + 2197 / 4104 * k4 - 1 / 5 * k5 )             //
//     where                                                                  //
//     k1 = f( x[i],y[i] ),                                                   //
//     k2 = f( x[i]+h/4, y[i] + h*k1/4 ),                                     //
//     k3 = f( x[i]+3h/8, y[i]+h*(3/32 k1 + 9/32 k2) ),                       //
//     k4 = f( x[i]+12h/13, y[i]+h*(1932/2197 k1 - 7200/2197 k2               //
//                                                       + 7296/2197 k3) ),   //
//     k5 = f( x[i]+h, y[i]+h*(439/216 k1 - 8 k2 + 3680/513 k3 - 845/4104 k4))//
//     k6 = f( x[i]+h/2, y[i]+h*(-8/27 k1 + 2 k2 - 3544/2565 k3               //
//                                              + 1859/4104 k4 - 11/40 k5) )  //
//     x[i+1] = x[i] + h.                                                     //
//                                                                            //
//     The error is estimated to be                                           //
//        err = h*( k1 / 360 - 128 k3 / 4275 - 2197 k4 / 75240 + k5 / 50      //
//              + 2 k6 / 55 )                                                 //
//     The step size h is then scaled by the scale factor                     //
//         scale = 0.8 * | epsilon * y[i] / [err * (xmax - x[0])] | ^ 1/4     //
//     The scale factor is further constrained 0.125 < scale < 4.0.           //
//     The new step size is h := scale * h.                                   //
////////////////////////////////////////////////////////////////////////////////

static double Runge_Kutta_4_5(double (*f)(double,double), double *y, double x,
                                                                   double h);

int Embedded_Fehlberg_4_5( double (*f)(double, double), double y[], double x,
                   double h ) {

   
         Runge_Kutta_4_5(f, y, x, h);
         
   return 0;
}


static const double a2 = 0.25;
static const double a3 = 0.375;
static const double a4 = 12.0 / 13.0;
static const double a6 = 0.5;
static const double b21 = 0.25;
static const double b31 = 3.0 / 32.0;
static const double b32 = 9.0 / 32.0;
static const double b41 = 1932.0 / 2197.0;
static const double b42 = -7200.0 / 2197.0;
static const double b43 = 7296.0 / 2197.0;
static const double b51 = 439.0 / 216.0;
static const double b52 = -8.0;
static const double b53 = 3680.0 / 513.0;
static const double b54 = -845.0 / 4104.0;
static const double b61 = -8.0 / 27.0;
static const double b62 = 2.0;
static const double b63 = -3544.0 / 2565.0;
static const double b64 = 1859.0 / 4104.0;
static const double b65 = -11.0 / 40.0;
static const double c1 = 25.0 / 216.0;
static const double c3 = 1408.0 / 2565.0;
static const double c4 = 2197.0 / 4104.0;
static const double c5 = -0.20;
static const double d1 = 1.0 / 360.0;
static const double d3 = -128.0 / 4275.0;
static const double d4 = -2197.0 / 75240.0;
static const double d5 = 0.02;
static const double d6 = 2.0 / 55.0;


static double Runge_Kutta_4_5(double (*f)(double,double), double *y, double x0,
                                                                   double h) {
   
   double k1, k2, k3, k4, k5, k6;
   double h2 = a2 * h, h3 = a3 * h, h4 = a4 * h, h6 = a6 * h;

   k1 = (*f)(x0, *y);
   k2 = (*f)(x0+h2, *y + h * b21 * k1);
   k3 = (*f)(x0+h3, *y + h * ( b31*k1 + b32*k2) );
   k4 = (*f)(x0+h4, *y + h * ( b41*k1 + b42*k2 + b43*k3) );
   k5 = (*f)(x0+h,  *y + h * ( b51*k1 + b52*k2 + b53*k3 + b54*k4) );
   k6 = (*f)(x0+h6, *y + h * ( b61*k1 + b62*k2 + b63*k3 + b64*k4 + b65*k5) );
   *(y+1) = *y +  h * (c1*k1 + c3*k3 + c4*k4 + c5*k5);
   return d1*k1 + d3*k3 + d4*k4 + d5*k5 + d6*k6;
}
