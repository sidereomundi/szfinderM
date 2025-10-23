//$Id: sub_fft2image.c,v 1.1 2012/06/04 14:03:20 jiayi Exp $
// fft fourier space data back to image

#include <stdlib.h>
#include <string.h>
#include "fftw3.h"
#include "routines.h"

// Basic fft fourier space data back to image
// work for fft() result
void fft2image(fdata, xdim, ydim, image)
  fftw_complex *fdata;
  int xdim,ydim;
  float *image;
{
  long nsize;
  int i;
  fftw_plan plan;
  double *dimage;
  
  /* ###### Begin inverse FFT ###### */
  nsize = xdim * ydim;
  dimage=(double *) calloc(nsize,sizeof(double));
  plan = fftw_plan_dft_c2r_2d(ydim,xdim,fdata,dimage,FFTW_ESTIMATE);
  fftw_execute(plan);
  
  /* ###### Prepare space for image ###### */
  for (i=0;i<nsize;i++)
    image[i]=dimage[i]/nsize;
  
  /* ###### Clean up ###### */
  free(dimage);
  fftw_destroy_plan(plan);
}

// Basic fft fourier space data back to image
// work for fftc() result
void fftc2image(fftw_complex *fdata, int xdim, int ydim, char* file)
{
  long nsize;
  int i;
  float *image;
  char outputimage[200];
  fftw_plan plan;
  fftw_complex *dimage;
  /* ###### Begin inverse FFT ###### */
  printf("  Begin inverse FFT for %s. \n",file);
  nsize = xdim * ydim;
  dimage=(fftw_complex *) calloc(nsize,sizeof(fftw_complex));
  plan = fftw_plan_dft_2d(ydim,xdim,fdata,dimage,FFTW_BACKWARD,FFTW_ESTIMATE);
  fftw_execute(plan);
  
  /* ###### Prepare space for image ###### */
  image = (float *) calloc(nsize,sizeof(float));
  if (image == NULL){
    printf("ERROR: no memory!\n");
		exit(1);
  }
  for (i=0;i<nsize;i++)
    image[i]=dimage[i][0];
  /* ###### Write the result ###### */
  if (strncmp(&(file[strlen(file)-5]),".fits",5)){
    sprintf(outputimage,"%s.fits",file);
  }
  else {
    sprintf(outputimage,"%s",file);
  }
  printf("  Writing the image for %s. \n",file);
  saveimagef(outputimage,image,xdim,ydim);
  
  fftw_complex *ffdata;
  ffdata=(fftw_complex*)fftw_malloc(ydim*(xdim/2+1)*sizeof(fftw_complex));
  fft(image,xdim,ydim,ffdata);
  fftsave(ffdata,ydim*(xdim/2+1),file);
  
  /* ###### Clean up ###### */
  fftw_free(ffdata);
  fftw_free(dimage);
  free(image);
  fftw_destroy_plan(plan);
}

void fft2imageaddnoise(fftw_complex *fdata,int xdim,int ydim,float *noise,float *image)
{
  long nsize;
  int i;
  fftw_plan plan;
  double *dimage;
  /* ###### Begin inverse FFT ###### */
  nsize = ydim*xdim;
  dimage=(double *) calloc(nsize,sizeof(double));
  plan = fftw_plan_dft_c2r_2d(ydim,xdim,fdata,dimage,FFTW_ESTIMATE);
  fftw_execute(plan);
  /* ###### Adding the noise ###### */
  for (i=0;i<nsize;i++){
    image[i]=dimage[i]/nsize+noise[i];
  }
  
  /* ###### Clean up ###### */
  free(dimage);
  fftw_destroy_plan(plan);
}
