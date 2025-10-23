// $Id: sub_tauvalue.c,v 1.2 2012/06/14 12:19:50 jiayi Exp $
// return cluster profile  (no beam effect)
/* Based on beta model (cf. Melin (2006) A&A 459.341   y(x)=y_0(1+|x|^2/theta_c^2)^-((3beta-1)/2) */
/* Return the beta profile of the cluster */

#include <math.h>

#include "parameters.h"
#include "routines.h"

double tauvalue(int x, int y, double theta)
{
  double temp1,temp2;
  double r;
  temp1=Squ(theta);
  temp2=(1-3*BETA)/2.0;
  r = Squ(x)+Squ(y);
  //  if( r > 100*temp1)      // truncated model
    return (pow(1+r/temp1,temp2));
  //else
  //  return 0;
}
