////////////////////////////////////////////////////////////////////////////////
// File: embedded_verner_6_7.c                                                //
// Routines:                                                                  //
//    Embedded_Verner_6_7                                                     //
////////////////////////////////////////////////////////////////////////////////

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

#include <math.h>

#define max(x,y) ( (x) < (y) ? (y) : (x) )
#define min(x,y) ( (x) < (y) ? (x) : (y) )

#define ATTEMPTS 12
#define MIN_SCALE_FACTOR 0.125
#define MAX_SCALE_FACTOR 4.0

static double Runge_Kutta(double (*f)(double,double), double *y, double x,
                                                                   double h);

////////////////////////////////////////////////////////////////////////////////
// int Embedded_Verner_6_7( double (*f)(double, double), double y[],          //
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
//                      calls to Embedded_Verner_6_7.                         //
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
//                                                                            //
int Embedded_Verner_6_7( double (*f)(double, double), double y[], double x,
                   double h, double xmax, double *h_next, double tolerance ) {

   static const double err_exponent = 1.0 / 6.0;

   double scale;
   double temp_y[2];
   double err;
   double yy;
   int i;
   int last_interval = 0;
   
      // Verify that the step size is positive and that the upper endpoint //
      // of integration is greater than the initial enpoint.               //

   if (xmax < x || h <= 0.0) return -2;
   
       // If the upper endpoint of the independent variable agrees with the //
       // initial value of the independent variable.  Set the value of the  //
       // dependent variable and return success.                            //

   *h_next = h;
   y[1] = y[0];
   if (xmax == x) return 0; 

       // Insure that the step size h is not larger than the length of the //
       // integration interval.                                            //
  
   h = min(h, xmax - x);

        // Redefine the error tolerance to an error tolerance per unit    //
        // length of the integration interval.                            //

   tolerance /= (xmax - x);

        // Integrate the diff eq y'=f(x,y) from x=x to x=xmax trying to  //
        // maintain an error less than tolerance * (xmax-x) using an     //
        // initial step size of h and initial value: y = y[0]            //

   temp_y[0] = y[0];
   while ( x < xmax ) {
      scale = 1.0;
      for (i = 0; i < ATTEMPTS; i++) {
         err = fabs( Runge_Kutta(f, temp_y, x, h) );
         if (err == 0.0) { scale = MAX_SCALE_FACTOR; break; }
         yy = (temp_y[0] == 0.0) ? tolerance : fabs(temp_y[0]);
         scale = 0.8 * pow( tolerance * yy /  err , err_exponent );
         scale = min( max(scale,MIN_SCALE_FACTOR), MAX_SCALE_FACTOR);
         if ( err < ( tolerance * yy ) ) break;
         h *= scale;
         if ( x + h > xmax ) h = xmax - x;
         else if ( x + h + 0.5 * h > xmax ) h = 0.5 * h;
      }
      if ( i >= ATTEMPTS ) { *h_next = h * scale; return -1; };
      temp_y[0] = temp_y[1];         
      x += h;
      h *= scale;
      *h_next = h;
      if ( last_interval ) break;
      if (  x + h > xmax ) { last_interval = 1; h = xmax - x; }
      else if ( x + h + 0.5 * h > xmax ) h = 0.5 * h;
   }
   y[1] = temp_y[1];
   return 0;
}


////////////////////////////////////////////////////////////////////////////////
//  static double Runge_Kutta(double (*f)(double,double), double *y,          //
//                                                       double x0, double h) //
//                                                                            //
//  Description:                                                              //
//     This routine uses Verner's embedded 6th and 7th order methods to       //
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

static double Runge_Kutta(double (*f)(double,double), double *y, double x0,
                                                                   double h) {
   
   static const double c1 = 7.0 / 90.0;
   static const double c4 = 32.0 / 90.0;
   static const double c5 = 32.0 / 90.0;
   static const double c7 = 12.0 / 90.0;
   static const double c8 = 7.0 / 90.0;

   static const double a2 = 1.0 / 12.0;
   static const double a3 = 1.0 / 6.0;
   static const double a4 = 1.0 / 4.0;
   static const double a5 = 3.0 / 4.0;
   static const double a6 = 16.0 / 17.0;
   static const double a7 = 1.0 / 2.0;
   static const double a9 = 2.0 / 3.0;

   static const double b41 = 1.0 / 16.0;
   static const double b43 = 3.0 / 16.0;
   static const double b51 = 21.0 / 16.0;
   static const double b53 = -81.0 / 16.0;
   static const double b54 = 72.0 / 16.0;
   static const double b61 = 1344688.0 / 250563.0;
   static const double b63 = -5127552.0 / 250563.0;
   static const double b64 = 4096896.0 / 250563.0;
   static const double b65 = -78208.0 / 250563.0;
   static const double b71 = -341549.0 / 234624.0;
   static const double b73 = 1407744.0 / 234624.0;
   static const double b74 = -1018368.0 / 234624.0;
   static const double b75 = 84224.0 / 234624.0;
   static const double b76 = -14739.0 / 234624.0;
   static const double b81 = -381875.0 / 136864.0;
   static const double b83 = 1642368.0 / 136864.0;
   static const double b84 = -1327872.0 / 136864.0;
   static const double b85 = 72192.0 / 136864.0;
   static const double b86 = 14739.0 / 136864.0;
   static const double b87 = 117312.0 / 136864.0;
   static const double b91 = -2070757.0 / 16755336.0;
   static const double b93 = 9929088.0 / 16755336.0;
   static const double b94 = 584064.0 / 16755336.0;
   static const double b95 = 3023488.0 / 16755336.0;
   static const double b96 = -447083.0 / 16755336.0;
   static const double b97 = 151424.0 / 16755336.0;
   static const double b10_1 = 130521209.0 / 10743824.0;
   static const double b10_3 = -499279872.0 / 10743824.0;
   static const double b10_4 = 391267968.0 / 10743824.0;
   static const double b10_5 = 13012608.0 / 10743824.0;
   static const double b10_6 = -3522621.0 / 10743824.0;
   static const double b10_7 = 9033024.0 / 10743824.0;
   static const double b10_9 = -30288492.0 / 10743824.0;
   
   static const double e1 = -1090635.0 / 172448640.0;
   static const double e4 = 9504768.0 / 172448640.0;
   static const double e5 = - 171816960.0 / 172448640.0;
   static const double e6 =  72412707.0 / 172448640.0;
   static const double e7 = - 55840512.0 / 172448640.0;
   static const double e8 = - 13412672.0 / 172448640.0;
   static const double e9 =  181730952.0 / 172448640.0;
   static const double e10 = - 21487648.0 / 172448640.0;

   double k1, k2, k3, k4, k5, k6, k7, k8, k9, k10;
   double h12 = a2 * h;
   double h6 = a3 * h;

   k1 = (*f)(x0, *y);
   k2 = (*f)(x0+h12, *y + h12 * k1);
   k3 = (*f)(x0+h6, *y + h6 * k2 );
   k4 = (*f)(x0+a4*h, *y + h * ( b41*k1 + b43*k3) );
   k5 = (*f)(x0+a5*h,  *y + h * ( b51*k1 + b53*k3 + b54*k4) );
   k6 = (*f)(x0+a6*h, *y + h * ( b61*k1 + b63*k3 + b64*k4 + b65*k5) );
   k7 = (*f)(x0+a7*h, *y + h * ( b71*k1 + b73*k3 + b74*k4 + b75*k5 + b76*k6) );
   k8 = (*f)(x0+h, *y + h * ( b81*k1 + b83*k3 + b84*k4 + b85*k5 + b86*k6
                                                                   + b87*k7) );
   k9 = (*f)(x0+a9*h, *y + h * ( b91*k1 + b93*k3 + b94*k4 + b95*k5 + b96*k6
                                                                   + b97*k7) );
   k10 = (*f)(x0+h, *y + h * ( b10_1*k1 + b10_3*k3 + b10_4*k4 + b10_5*k5
                                          + b10_6*k6 + b10_7*k7 + b10_9*k9 ) );
   *(y+1) = *y +  h * (c1 * k1 + c4 * k4 + c5 * k5 + c7 * k7 + c8 * k8);
   return e1*k1 + e4*k4 + e5*k5 + e6*k6 + e7*k7 + e8*k8 + e9*k9 + e10*k10;
}
