//$Id: sub_calcbeta.c,v 1.1 2012/06/01 09:27:37 jiayi Exp $
// Calculate beta model and the fft

/* Return the tau profile of the cluster */
double tauvalue(int x, int y, double theta)
{
  double temp1,temp2;
  double r;
  temp1=Squ(theta);
  temp2=(1-3*beta)/2.0;
  r = Squ(x)+Squ(y);
  if( r < 10000*temp1)      // truncated model
    return (pow(1+r/temp1,temp2));
  else
    return 0;
}


// build beta model given map size and coresize.
// the beta map is pre-existed.
void buildbeta(float coresize,struct simmap* map, double *betamap)
{
  
}

void buildbetafft()
{
}
