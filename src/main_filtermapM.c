//$Id: main_filtermapM.c,v 1.3 2012/11/27 13:50:16 jiayi Exp $
// Multi-frequency filter

#include "fftw3.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "parameters.h"
#include "routines.h"
#include "structures.h"

int main(int argc, char* argv[])
{
  char imagename150[200],imagename90[200],imagename[200],filtername[200];
  struct simmap fmap;  //store filtered map
  float *image150,*image90; // store observed maps in frequency
  int xdim,ydim;
  long pix,nsize,ntot;
  double mean,vari; // check the map mean and sqrt(variance)
  fftw_complex *ffilter150,*ffilter90,*fimage90,*fimage150;
  int Kcontrol;
  float theta, sigma0;
  FILE *fsig;

  if(argc<2){
    printf(" # Using default file szmock[freq].fits \n");
    sprintf(imagename150,"%s","szmock150.fits");
    sprintf(imagename90,"%s","szmock9.fits");
  }
  else
  {
    sprintf(imagename150,"%s150.fits",argv[1]);
    sprintf(imagename90,"%s95.fits",argv[1]);
  }
  // preparing maps
  printf(" # preparing maps\n");
  initializemap(&fmap);
  readimage(imagename150,&image150,&xdim,&ydim);
  if((xdim!=fmap.nsidepix)||(ydim!=fmap.nsidepix)){
    printf("size of image is wrong!\n");
    exit(1);
  }
  readimage(imagename90,&image90,&xdim,&ydim);
  if((xdim!=fmap.nsidepix)||(ydim!=fmap.nsidepix)){
    printf("size of image is wrong!\n");
    exit(1);
  }
  nsize=fmap.npix;
  mean =0;
  for(pix=0;pix<nsize;pix++)
    mean+=image150[pix];
  mean /= nsize;
  for(pix=0;pix<nsize;pix++)
    image150[pix]-=mean; // subtract 0 frequency mode
  printf("mean of CMB is %lf (in 150 Ghz), ",mean);
  mean =0;
  for(pix=0;pix<nsize;pix++)
    mean+=image90[pix];
  mean /= nsize;
  for(pix=0;pix<nsize;pix++)
    image90[pix]-=mean; // subtract 0 frequency mode
  printf(" and %lf (in 90 Ghz)\n",mean);

  // fft the map
  ntot=ydim*(xdim/2+1);
  ffilter150=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*ntot);
  fimage150=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*ntot);
  fft(image150,xdim,ydim,fimage150);
  ffilter90=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*ntot);
  fimage90=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*ntot);
  fft(image90,xdim,ydim,fimage90);
  // load sigma file
  fsig = fopen(FILTERMSIGFILE,"r");
  if(fsig==NULL){
    printf("fail to open %s\n!",FILTERMSIGFILE);
    exit(1);
  }
  for(Kcontrol=0;Kcontrol<NFILTER;Kcontrol++)
  {
    // load filter
    printf("# Filtering case %d\n",Kcontrol);
    fscanf(fsig,"%f %e", &theta, &sigma0);
    sprintf(filtername,"%sM150_%d_fft.dat",FILTERFILE,Kcontrol);
    fftread(ntot,ffilter150,filtername);
    sprintf(filtername,"%sM90_%d_fft.dat",FILTERFILE,Kcontrol);
    fftread(ntot,ffilter90,filtername);
    // filtering
    filter(fimage150,ffilter150,xdim,ydim,image150);
    filter(fimage90,ffilter90,xdim,ydim,image90);
    // normalizing
    mean = 0;
    vari = 0;
    for (pix=0;pix<nsize;pix++)
    {
      fmap.img[pix]=-(image90[pix]+image150[pix]);
      mean += fmap.img[pix];
      vari += Squ(fmap.img[pix]);
    }
    mean = mean/(float)nsize;
    vari = sqrt(vari/(float)nsize - mean*mean);
    printf("case %d:: %e %e %e\n",Kcontrol,mean,vari,sigma0);
    for (pix=0; pix<nsize; pix++) {
      //      fmap.img[pix]/=sigma0;	
      fmap.img[pix]/=vari;
    }
    // save map
    sprintf(imagename,"filteredM%02d.fits",Kcontrol);
     savemap(&fmap,imagename);
  }

  // clean up
  fclose(fsig);
  free(image150);  free(image90);
  fftw_free(ffilter150);  fftw_free(ffilter90);
  fftw_free(fimage150);  fftw_free(fimage90);
  cleanmap(&fmap);

  return 0;
}
