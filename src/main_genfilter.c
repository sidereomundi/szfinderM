// $Id: main_genfilter.c,v 1.6 2013/06/14 09:40:38 jiayi Exp $
// generate the single frequency filter

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "fftw3.h"
#include "nr.h"
#include "structures.h"
#include "routines.h"
#include "parameters.h"

int main(int argc, char *argv[])
{
  char outputimage[200],beamimage[200];
  
  double theta;
  double pix2arcmin;
  double dT,Pnoise;
  
  fftw_complex *fbeam,*idata,*fdata;
  float *beamfile;
  double sigma;
  double Pbeam,Ptot;
  float Pcmb;
  
  int xdim,ydim;
  long ntot;
  int i,x,y,k;
  int halfdim,rdim;
  double size;
  
  int cllength;
  float *cl_l,*cl_p,*cl_py2,ypn,yp1,l;
  
  FILE *pfile;

  double temp1;

  struct simmap fmap;
  initializemap(&fmap);
  xdim=fmap.nsidepix; ydim=fmap.nsidepix;
  ntot=fmap.npix;
  size = ntot;
  halfdim = fmap.nsidepix/2;
  rdim = fmap.nsidepix;
  pix2arcmin = fmap.pixsize; 
  // work for 150 GHz
  /* Noise term */
  // work in power-density domain
  // devide the whole sky power by 4\pi
  dT=DT150*1.0e-6;	// K/arcmin
  dT/=pix2arcmin;		// K/arcmin  -> K/pixel
  Pnoise=Squ(dT)*Squ(pix2arcmin*M_PI/(180.*60.));  // power density

  // frequency base unit in l space
  temp1=fmap.nsidepix*pix2arcmin*M_PI/(180.0*60.0)/(2.0*M_PI);
  
  /* CMB term */
  /* read in from Cl.dat */
  /* Ref. mkcmb.c */
/* read in CMB power spectrum */
  printf(" # Reading CMB power spectrum from %s\n",clfile); fflush(stdout);
  clinput(&cl_l,&cl_p,&cl_py2,&yp1,&ypn,&cllength);

/* Beam term */
  sprintf(beamimage,"%s.fits",BEAM150FILE);
  readimage(beamimage,&beamfile,&xdim,&ydim);
  if ((xdim*ydim)!=ntot){
    printf("size of beam image is wrong!\n");
    exit(0);
  }
  fbeam=(fftw_complex*)fftw_malloc(ntot*sizeof(fftw_complex));
  fftc(beamfile,xdim,ydim,fbeam);
  
  for(i=0;i<ntot;i++) {
      fbeam[i][0]*=(float)ntot;
      fbeam[i][1]*=(float)ntot;
  } // normalize to power density. 
  
  pfile=fopen(FILTERSIGFILE,"w"); // write the filter variance for reference.
  
  // go through all filter size k 0~11
  for(k=0;k<NFILTER;k++){
    /* ####### Generate the tau profile ###### */
    /* Cluster radius */
    theta=thetasize(k)/pix2arcmin;
    fprintf(pfile,"%f ",theta);
    printf("# Calculating %f(arcmin)/%f(pixel) case\n",thetasize(k),theta);
    
    /* Based on beta model (cf. Melin (2006) A&A 459.341
       y(x)=y_0(1+|x|^2/theta_c^2)^-((3beta-1)/2) */
    for (i=0;i<size;i++) {
      y=i/xdim; if(y>halfdim) y = rdim-y;
      x=i%xdim; if(x>halfdim) x = rdim-x;
      fmap.img[i]=tauvalue(x,y,theta);
    }
    // output temporary beta profile
    // sprintf(outputimage,"beta_%d.fits",k);
    // savemap(&fmap,outputimage);

    idata=(fftw_complex*)fftw_malloc(ntot*sizeof(fftw_complex));
    fftc(fmap.img,xdim,ydim,idata);  // normalization is irrelevant  

    /* ####### Generate the filter profile ####### */
    sigma=0;
    fdata=(fftw_complex*)fftw_malloc(ntot*sizeof(fftw_complex));
    /* !! When the power of noise is zero, the filter will
       crush because the Pbeam goes to zero in high frequency.
       Analitically, the tauprofile goes to zero too.  !! */
    for (i=0;i<ntot;i++){
      y=i/xdim; if(y>halfdim) y = rdim-y;
      x=i%xdim; if(x>halfdim) x = rdim-x;
      l=sqrt(Squ((float)x)+Squ((float)y))/temp1;
      if (l==0)
	Pcmb=0;  // no cmb contribution at l=0.  remove sz map 0-mode component
      else if (l<cl_l[cllength])
	splint(cl_l,cl_p,cl_py2,cllength,l,&Pcmb);
      else Pcmb=cl_p[cllength];  // total power density at l mode
      Pcmb = Squ(Pcmb); 

      // Change for Melin 2005 paper, beamed tau
      Pbeam=(Squ(fbeam[i][0])+Squ(fbeam[i][1]));
      Ptot=((Pcmb)*Pbeam+Pnoise);
      //      sigma+=(Squ(idata[i][0]*fbeam[i][0])+Squ(idata[i][1]*fbeam[i][1]))/Ptot;
      fdata[i][0]=(idata[i][0]*fbeam[i][0]-idata[i][1]*fbeam[i][1])/Ptot;
      fdata[i][1]=(idata[i][0]*fbeam[i][1]+idata[i][1]*fbeam[i][0])/Ptot;
      sigma += (Squ(fdata[i][0])+Squ(fdata[i][1]))*Ptot;
      // imaginary part is tiny due to the symmetry.  verified by following line.
      //        if(l<10000)
      //       printf("%e %e\n",fdata[i][0],fdata[i][1]);
      //        printf("%e %e %e %e\n",l,Pcmb*Pbeam,Pnoise,(Squ(fdata[i][0])+Squ(fdata[i][1])));
    }
    // trying to normalize the filter.  // not check yet
    for (i=0;i<ntot;i++){
      fdata[i][0]/=(sigma*ntot);
      fdata[i][1]/=(sigma*ntot);
    }
    sigma=sqrt(1.0/sigma);
//    printf("filter sigma: %e\n",sigma*sqrt(2)*M_PI);
    fprintf(pfile,"%e\n",sigma*sqrt(2)*M_PI);
    
    /* ###### Output the filter fft data & image ###### */
    sprintf(outputimage,"%s150_%d",FILTERFILE,k);
    fftc2image(fdata,xdim,ydim,outputimage);  //output image.  use .dat only.
    // see sub_fft2image.c for detail.
    fftw_free(idata);
    fftw_free(fdata);
    }
  /* ###### Clear up ###### */
  fclose(pfile);
  fftw_free(fbeam);
  free(beamfile);
  free(cl_l);
  free(cl_p);
  free(cl_py2);
  return 0;
}


