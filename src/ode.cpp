# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "../include/ode.hpp"

//****************************************************************************80

int i4_sign ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    I4_SIGN returns the sign of an I4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 March 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the integer whose sign is desired.
//
//    Output, int I4_SIGN, the sign of I.
{
  int value;

  if ( i < 0 )
  {
    value = -1;
  }
  else
  {
    value = 1;
  }
  return value;
}
//****************************************************************************80

void intrp 
( 
  double x,
  double y[],
  double xout,
  double yout[],
  double ypout[], 
  int neqn,
  int kold,
  double phi[],
  double psi[] 
)

//****************************************************************************80
//
//  Purpose:
//
//    INTRP approximates the solution at XOUT by polynomial interpolation.
//
//  Discussion:
//
//    The methods in STEP approximate the solution near X by a polynomial.
//    This routine approximates the solution at XOUT by evaluating the
//    polynomial there.  Information defining this polynomial is passed
//    from STEP, so INTRP cannot be used alone.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 January 2012
//
//  Author:
//
//    Original FORTRAN77 version by Lawrence Shampine, Marilyn Gordon.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Lawrence Shampine, Marilyn Gordon,
//    Computer Solution of Ordinary Differential Equations:
//    The Initial Value Problem,
//    Freeman, 1975,
//    ISBN: 0716704617,
//    LC: QA372.S416.
//
//  Parameters:
//
//    Input, double X, the point where the solution has been computed.
//
//    Input, double Y[NEQN], the computed solution at X.
//
//    Input, double XOUT, the point at which the solution is desired.
//
//    Output, double YOUT[NEQN], the solution at XOUT.
//
//    Output, double YPOUT[NEQN], the derivative of the solution
//    at XOUT.
//
//    Input, int NEQN, the number of equations.
//
//    Input, int KOLD, the order used for the last
//    successful step.
//
//    Input, double PHI[NEQN*16], contains information about the
//    interpolating polynomial.
//
//    Input, double PSI[12], contains information about the
//    interpolating polynomial.
//
{
  double eta;
  double g[13];
  double gamma;
  double hi;
  int i;
  int j;
  int k;
  int ki;
  int l;
  double psijm1;
  double rho[13];
  double term;
  double w[13];

  hi = xout - x;
  ki = kold + 1;
//
//  Initialize W for computing G.
//
  for ( i = 1; i <= ki; i++ )
  {
    w[i-1] = 1.0 / ( double ) ( i );
  }
//
//  Compute G.
//
  g[0] = 1.0;
  rho[0] = 1.0;
  term = 0.0;

  for ( j = 2; j <= ki; j++ )
  {
    psijm1 = psi[j-2];
    gamma = ( hi + term ) / psijm1;
    eta = hi / psijm1;
    for ( i = 1; i <= ki + 1 - j; i++ )
    {
      w[i-1] = gamma * w[i-1] - eta * w[i];
    }
    g[j-1] = w[0];
    rho[j-1] = gamma * rho[j-2];
    term = psijm1;
  }
//
//  Interpolate.
//
  for ( k = 0; k < neqn; k++ )
  {
    ypout[k] = 0.0;
    yout[k] = 0.0;
  }

  for ( j = 1; j <= ki; j++ )
  {
    i = ki + 1 - j;
    for ( k = 0; k < neqn; k++ )
    {
      yout[k] = yout[k] + g[i-1] * phi[k+(i-1)*neqn];
      ypout[k] = ypout[k] + rho[i-1] * phi[k+(i-1)*neqn];
    }
  }

  for ( k = 0; k < neqn; k++ )
  {
    yout[k] = y[k] + hi * yout[k];
  }
  return;
}
//****************************************************************************80



double r8_epsilon ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_EPSILON returns the R8 roundoff unit.
//
//  Discussion:
//
//    The roundoff unit is a number R which is a power of 2 with the
//    property that, to the precision of the computer's arithmetic,
//      1 < 1 + R
//    but
//      1 = ( 1 + R / 2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 September 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_EPSILON, the R8 round-off unit.
//
{
  const double value = 2.220446049250313E-016;

  return value;
}
//****************************************************************************80

double r8_max ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MAX returns the maximum of two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MAX, the maximum of X and Y.
//
{
  double value;

  if ( y < x )
  {
    value = x;
  }
  else
  {
    value = y;
  }
  return value;
}
//****************************************************************************80

double r8_min ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MIN returns the minimum of two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MIN, the minimum of X and Y.
//
{
  double value;

  if ( y < x )
  {
    value = y;
  }
  else
  {
    value = x;
  }
  return value;
}
//****************************************************************************80

double r8_sign ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_SIGN returns the sign of an R8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the number whose sign is desired.
//
//    Output, double R8_SIGN, the sign of X.
//
{
  double value;

  if ( x < 0.0 )
  {
    value = -1.0;
  }
  else
  {
    value = 1.0;
  }
  return value;
}
//****************************************************************************80


void timestamp ( )

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    31 May 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    None
//
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct std::tm *tm_ptr;
  size_t len;
  std::time_t now;

  now = std::time ( NULL );
  tm_ptr = std::localtime ( &now );

  len = std::strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm_ptr );

  std::cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}