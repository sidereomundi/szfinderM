// $Id: main_mkbeam.c,v 1.3 2012/07/10 13:08:50 jiayi Exp $
// changed from main_mkgauss.c


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "fftw3.h"
#include "routines.h"
#include "structures.h"
// build gaussian shape beam profile
// input:
//        fwhm [arcmin]
//        output filename
int main(int argc, char *argv[])
{
  float fwhm;
  float sigma,scale;
  long size,x,y,r,i;
  int xdim,ydim;
  struct simmap beammap;
  char outputimage[200]; //output image file
  
  fftw_complex *fdata;
  long ntot;
  
  float nor;
  if (argc < 4)  {
    printf("mkgaussian <fwhm[arcmin]> <output_name> [Normalization]\n");
    printf("using default normalization 1\n");
    nor = 1.0;
    if (argc <2 ) exit (0);
  }
  if (argc == 4)
  {
    sscanf(argv[3],"%f",&nor);
    printf("using normalization %f\n",nor); 
  }
  sscanf(argv[1],"%f",&fwhm);
  printf("%f\n",fwhm);
  initializemap(&beammap);
  fwhm/=beammap.pixsize;
  printf("# the beam size is %f pixels.\n",fwhm);
  sigma=fwhm/2.354;	// Based on Amy's code!
  
  xdim=beammap.nsidepix; ydim=beammap.nsidepix;
  size = beammap.npix;
  sigma=Squ(sigma);
  scale=1.0/(2.0*M_PI)/sigma*nor;// Based on Amy's code (Normalized)
  for (i=0;i<size;i++) {
    y=i/xdim; if (y>ydim/2) y=ydim-y;
    x=i%xdim; if (x>xdim/2) x=xdim-x;
    r=Squ(x)+Squ(y);
    beammap.img[i]=0.5*r/sigma;	// Based on Amy's code (r^2/2sigma^2)
    //    if (beammap.img[i]<40.0)
      beammap.img[i]=scale*exp(-beammap.img[i]);
    //    else beammap.img[i]=0.0;
  }
  
  /* ###### output image ###### */
  sprintf(outputimage,"%s.fits",argv[2]);
  savemap(&beammap, outputimage);
  
  /* ###### output fft data ###### */
  ntot=ydim*(xdim/2+1);
  fdata=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*ntot);
  fft(beammap.img,xdim,ydim,fdata);
  fftsave(fdata,ntot,argv[2]);
  fftw_free(fdata);
  cleanmap(&beammap);
  return 0;
}
