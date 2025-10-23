// $Id: sub_mkdirtmap.c,v 1.2 2012/11/27 13:50:16 jiayi Exp $
// add noise and beam effect

#include "fftw3.h"
#include "routines.h"
#include "structures.h"

void mkdirtmap(struct simmap* map, float *noise, fftw_complex *fbeam)
{
  int xdim,ydim,ntot;
  fftw_complex *fdata;
  double fre,fim;
  int i;
  
  /* ###### Convolve with the beam ###### */
  xdim=map->nsidepix;	ydim=map->nsidepix;
  ntot=ydim*(xdim/2+1);
  fdata=(fftw_complex *)fftw_malloc(sizeof(fftw_complex)*ntot);
  fft(map->img,xdim,ydim,fdata);
  for(i=0;i<ntot;i++) {
    fre=fdata[i][0];
    fim=fdata[i][1];
    fdata[i][0]=fre*fbeam[i][0]-fim*fbeam[i][1];
    fdata[i][1]=fre*fbeam[i][1]+fim*fbeam[i][0];
  }
  /* ###### Inverse FFT & Output the dirty map ###### */
  fft2imageaddnoise(fdata,xdim,ydim,noise,map->img);
  /* ###### Clear up ###### */
  fftw_free(fdata);
}

void removemean(struct simmap* map)
{
  long i,ntot;
  double mean;
  ntot = Squ(map->nsidepix);
  mean = 0;
  for(i=0;i<ntot;i++)
    mean += map->img[i];
  mean /=ntot;
  for(i=0;i<ntot;i++)
    map->img[i] -= mean;
}
