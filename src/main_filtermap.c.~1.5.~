// $Id: main_filtermap.c,v 1.5 2012/06/20 22:11:58 jiayi Exp $
// filtering the sz+noise map with pre-build filter profile.

// !! note !!
// the normalization is based on theoritical prediction
// which is lower for large theta_core case due to the smoothing scale

#include "fftw3.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "parameters.h"
#include "routines.h"
#include "structures.h"

int main(int argc, char *argv[])
{
  FILE *fsig;
  int Kcontrol;
  int ntot,xdim,ydim;
  float sigma0,theta;
  float *image;
  double mean,vari;
  char filtername[200],imagename[200];
  long pix,nsize;
  fftw_complex *ffilter,*fimage;
  struct simmap fmap;
  // general initialization
  if(argc<2){
    printf(" # Using default file szmock150.fits \n");
    sprintf(imagename,"%s","szmock150.fits");
  }
  else
    sscanf(argv[1],"%s",imagename);
  initializemap(&fmap);
  readimage(imagename,&image,&xdim,&ydim);
  if((xdim!=fmap.nsidepix)||(ydim!=fmap.nsidepix)){
    printf("size of image is wrong!\n");
    exit(1);
  }
  nsize=fmap.npix;
  mean =0;
  for(pix=0;pix<nsize;pix++)
    mean+=image[pix];
  mean /= nsize;
  for(pix=0;pix<nsize;pix++)
    image[pix]-=mean; // subtract 0 frequency mode
  
  printf("mean of CMB is %lf\n",mean);
  
  ntot=ydim*(xdim/2+1);
  ffilter=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*ntot);
  fimage=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*ntot);
  fft(image,xdim,ydim,fimage);
  fsig = fopen(FILTERSIGFILE,"r");
  if(fsig==NULL){
    printf("fail to open %s\n!",FILTERSIGFILE);
    exit(1);
  }
  for(Kcontrol=0;Kcontrol<NFILTER;Kcontrol++){
    // initializing case
    printf("# Initialize case %d\n",Kcontrol);
    fscanf(fsig,"%f %e", &theta, &sigma0);
    sprintf(filtername,"%s150_%d_fft.dat",FILTERFILE,Kcontrol);
    sprintf(imagename,"filtered%02d.fits",Kcontrol);
    fftread(ntot,ffilter,filtername);
    printf("  # Filtering \n");
    filter(fimage,ffilter,xdim,ydim,fmap.img);
    printf("  # Normalizeing \n");
    for(pix=0;pix<nsize;pix++)
      fmap.img[pix]/=-1;//sigma0;
    // image normalizing
    mean =0;
    vari =0;
    for (pix=0; pix< nsize; pix++) {
      mean += fmap.img[pix];
      vari += Squ(fmap.img[pix]);
    }
    mean = mean/(float)nsize;
    vari = sqrt(vari/(float)nsize - mean*mean);
    printf("case %d:: %e %e %e\n",Kcontrol,mean,vari,sigma0);
    for (pix=0; pix<nsize; pix++) {
 //     fmap.img[pix]/=sigma0;	
      fmap.img[pix]/=vari;	
    }
    savemap(&fmap,imagename);
  }
  fclose (fsig);
  return 0;
}

