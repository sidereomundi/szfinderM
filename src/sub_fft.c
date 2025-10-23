// $Id: sub_fft.c,v 1.2 2012/06/07 16:24:40 jiayi Exp $
// read in the frequency components written by filter.c mkdirmap.c

#include "fftw3.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Read in the complex fft component data written by fftsave() */
void fftread(long ntot, fftw_complex *fft, char *name)
{
  FILE *pfile;
  long i;
  long lSize;

  if (strncmp("none",name,4)) {
    pfile=fopen(name,"r");
    if(pfile==NULL){
      printf("reading %s failed\n",name);
      exit(1);
    }
    fseek(pfile,0,SEEK_END);
    lSize = ftell(pfile);
    rewind(pfile);
    if (lSize != ntot*sizeof(fftw_complex)){
      printf(" The size of %s is incorrect:%ld (%ld) \n",name,lSize,ntot*sizeof(fftw_complex));
      exit(1);
    }
    fread(fft,sizeof(fftw_complex),ntot,pfile);
    fclose(pfile);	
  }
  else{
    printf("No fft and fftspeq intput\n");
    for(i=0;i<ntot;i++){
      (fft)[i][0]=1;
      (fft)[i][1]=0;
    }
  }
}

// convert real image to complex 2D array (half complex plane)
void fft(float *idata, int xdim, int ydim, fftw_complex *fdata)
{
  double *dimage;
  long size,i;
  fftw_plan plan;
  size=xdim*ydim;
  dimage=(double *)calloc(sizeof(double),size);
  for (i=0;i<size;i++)
    dimage[i]=idata[i];
  /* set up Fourier array	*/
  plan=fftw_plan_dft_r2c_2d(ydim,xdim,dimage,fdata,FFTW_ESTIMATE);
  fftw_execute(plan);
  /* clean up */
  fftw_destroy_plan(plan);
  fftw_free(dimage);
}

// convert real image to complex 2D array (full complex plane)
void fftc(float *idata, int xdim, int ydim, fftw_complex *fdata)
{
	fftw_complex *dimage;
	int size,i;
	fftw_plan plan;
	size=xdim*ydim;
	dimage=(fftw_complex *)fftw_malloc(size*sizeof(fftw_complex));
	for (i=0;i<size;i++) {
		dimage[i][0]=idata[i];
		dimage[i][1]=0;
	}
/* set up Fourier array	*/
	plan=fftw_plan_dft_2d(ydim,xdim,dimage,fdata,FFTW_FORWARD,FFTW_ESTIMATE);
	fftw_execute(plan);
	for (i=0;i<size;i++){
		fdata[i][0]/=size;
		fdata[i][1]/=size;
	}
/* clean up */
	fftw_destroy_plan(plan);
	fftw_free(dimage);
}
