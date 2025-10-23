// $Id: sub_filter.c,v 1.2 2012/06/13 09:25:05 jiayi Exp $
// filter the map with filter
// rely on fft2image() in sub_fft2image.c

#include "fftw3.h"
#include "routines.h"

// changed to convolve with the conjugate of filter
void filter(fftw_complex *fdata, fftw_complex *ffilter,int xdim, int ydim, float *image)
{
  fftw_complex *temp;
  int i,ntot;
  
  ntot=ydim*(xdim/2+1);
  temp = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*ntot);
  
  /* ###### Convolve with the filter ###### */
  for(i=0;i<ntot;i++) {
//    temp[i][0]=fdata[i][0]*ffilter[i][0]-fdata[i][1]*ffilter[i][1];
//    temp[i][1]=fdata[i][0]*ffilter[i][1]+fdata[i][1]*ffilter[i][0];
    temp[i][0]=fdata[i][0]*ffilter[i][0]+fdata[i][1]*ffilter[i][1];
    temp[i][1]=-fdata[i][0]*ffilter[i][1]+fdata[i][1]*ffilter[i][0];
  }
  /* ###### Inverse FFT & Output the filtered map ###### */
  fft2image(temp,xdim,ydim,image);
  fftw_free(temp);
}

