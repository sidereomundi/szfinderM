// $Id: main_mknull.c,v 1.1 2012/06/04 14:03:20 jiayi Exp $
// generate sz-null y-map  (fill map with 0)

#include <stdio.h>
#include <stdlib.h>
#include "parameters.h"
#include "routines.h"

int main(int argc, char *argv[])
{
  float *img;
  long npix,i;
  char fname[200];
  npix=NDIM*NDIM;
  img = (float*)malloc(npix*sizeof(npix));
  for (i=0;i<npix;i++)
    img[i]=0;
  sprintf(fname,"sznull.fits");
  saveimagef(fname,img,NDIM,NDIM);
  return 0;
}
