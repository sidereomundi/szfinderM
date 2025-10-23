// $Id: sub_mkcmb.c,v 1.1 2012/06/04 14:03:20 jiayi Exp $
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nr.h"
#include "fftw3.h"
#include "structures.h"
#include "parameters.h"
#include "routines.h"

// create differential CMB temperature map
void mkcmb(cllength,seed,imap,cmbmap,cl_l,cl_p,cl_py2,yp1,ypn)
     int cllength;
     long *seed;
     struct simmap  *imap;
     float *cmbmap,*cl_l,*cl_p,*cl_py2,yp1,ypn;
{
  float	l,lmax=0.0,lmin=9999999.9,pixsize,dT;	
  float	amplitude_r,amplitude_i;
  int	i,loc,length,u,v,u2,v2,loc2,halfimsize,vmin;
  
  fftw_complex *dimage,*fdata;
  fftw_plan plan;
  float  *amplitudes,*amplitudes_r,*amplitudes_i;
  double testtemp;
  int imsize,xdim,ydim;
 
  // initial setup
  imsize=imap->nsidepix;
  pixsize=imap->pixsize;
  length=imap->npix;

  pixsize*=M_PI/(180.0*60.0);
  /* create an array to hold the final image */
  //	printf(" Creating %d^2 single precision image template\n",imsize); 
  dimage=(fftw_complex*)fftw_malloc(length*sizeof(fftw_complex));
  
  fdata=(fftw_complex*)fftw_malloc(length*sizeof(fftw_complex));
  amplitudes = (float *)malloc(sizeof(float)*length);
  amplitudes_r = (float *)malloc(sizeof(float)*length);
  amplitudes_i = (float *)malloc(sizeof(float)*length);
  xdim=imsize; ydim=imsize;
  halfimsize=imsize/2;
  testtemp=imsize*pixsize;
  
  /* fill these arrays with Gaussian Random field */
  for (u=0;u<=halfimsize;u++) { 
    if (u!=0 && u!=halfimsize) vmin=-halfimsize+1;
    else vmin=0;
    for (v=vmin;v<=halfimsize;v++) { 
      /* calculate the location of mode at (u,v) */
      if (v>=0) loc=u+v*imsize;
      else loc=u+(v+imsize)*imsize;
      /* calculate distance from origin */
      l=sqrt(Squ((float)u)+Squ((float)v))/(testtemp);
      l*=2.*M_PI;
      if(l>lmax) lmax=l;
      /* we will set the monopole (u=0 and v=0 => l=0) to
       * zero before doing inverse Fourier transform. */
      if((u!=0&&v!=0)&&l<lmin) lmin=l;
      /* interpolate to get the correct dT */
      if (l<cl_l[cllength]) 
	splint(cl_l,cl_p,cl_py2,cllength,l,&dT);
      else dT=cl_p[cllength];
      /* get Gaussian random deviate using this dT */
      // normalize (imsize) and scaling (pixsize)
      amplitude_r=dT*gasdev(seed)/sqrt(2.)/((float)imsize)/pixsize;
      amplitude_i=dT*gasdev(seed)/sqrt(2.)/((float)imsize)/pixsize;
      amplitudes_r[loc]=amplitude_r;
      amplitudes_i[loc]=amplitude_i;
      /* set real component of Fourier mode */
      fdata[loc][0]=amplitude_r;
      /* set imaginary component of Fourier mode */
      fdata[loc][1]=amplitude_i;
      /* calculate the location of mode at (-u,-v) */
      u2=-u;v2=-v;
      if (u2<0) u2+=imsize;
      if (v2<0) v2+=imsize;
      loc2=u2+v2*imsize;
      /* set real component of Fourier mode */
      fdata[loc2][0]=fdata[loc][0];
      /* set imaginary component of Fourier mode */
      fdata[loc2][1]=-1.0*fdata[loc][1];
      amplitudes_r[loc2]=amplitude_r;
      amplitudes_i[loc2]=-amplitude_i;
    }
  }

  /* set u=v=0 component by hand */
  /*   data[0]=0 because we should have subtracted the monopole
      before calculating the spectrum (when u=v=0, l=0 -> monopole).
      We will add T_cmb to the final image */
  /* data[1]=0 because the final image is real (see Numerical 
      Recipes section 12-0 */
  fdata[0][0]=fdata[0][1]=0.0;
  
  /* set (u=0,v=N/2), (u=N/2,v=0) and (u=N/2,v=N/2) imaginary 
      components to zero */
  /* Note that loc=u+v*imsize, and the imaginary part is at 2*loc+1 */
  fdata[halfimsize*imsize][1]=0.0;
  fdata[halfimsize][1]=0.0;
  fdata[(halfimsize+halfimsize*imsize)][1]=0.0;
  
  /////////////////////////////////////////////////////
  /* Inverse transform to get CMB image */
  plan=fftw_plan_dft_2d(ydim,xdim,fdata,dimage,FFTW_BACKWARD,FFTW_ESTIMATE);
  fftw_execute(plan);
  
  /* copy real part of this image into image template */
  for (i=0;i<length;i++) {
    cmbmap[i]=dimage[i][0]+T_CMB;
  }
  
  fftw_destroy_plan(plan);
  fftw_free(fdata);
  fftw_free(dimage);
}
