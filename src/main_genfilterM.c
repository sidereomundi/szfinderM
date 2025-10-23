// $Id: main_genfilterM.c,v 1.5 2013/02/12 14:23:15 jiayi Exp $
// generate the multi-frequency filter

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "fftw3.h"
#include "nr.h"
#include "structures.h"
#include "routines.h"
#include "parameters.h"

double convertNoise(double, double);  // convert dT to fourier space // at the bottom

int main(int argc, char *argv[])
{
  int xdim,ydim,halfdim,rdim;
  long ntot;
  double pix2arcmin;
  struct simmap fmap;

  double Pnoise150, Pnoise90;
  double templ;  // previous temp1 in main_genfilter.c
  // cmb part
  int cllength;
  float *cl_l,*cl_p,*cl_py2,ypn,yp1,l;
  // beam part
  float *beam;
  char beamfile[200];
  fftw_complex *fbeam150,*fbeam90;
  // filter profile
  FILE *pfile;
  int k,ix,iy;
  long i;
  double theta;
  fftw_complex *idata,*fdata150,*fdata90;
  double sigma;
  float Pcmb;
  double Pbeam150,Pbeam90,b1b2real,b1b2image,P90,P150;
  double Pa,Pc_real,Pc_image,Pd,Ptot;
  double j150,j90,tempx,tempexpx;
  char outputfile[300];
  // map setup
  initializemap(&fmap);
  xdim = fmap.nsidepix; 
  ydim = fmap.nsidepix;
  ntot = fmap.npix;
  halfdim = fmap.nsidepix/2;
  rdim = fmap.nsidepix;
  pix2arcmin = fmap.pixsize;
  
  // initialize parameters
  /* noise term */
  if (argc < 3)
  {
    printf("!! using default noise %f@150, %f@95\n",DT150,DT90);
    Pnoise150 = convertNoise(DT150,pix2arcmin);
    Pnoise90 =  convertNoise(DT90,pix2arcmin);
  }
  else
  {
    sscanf(argv[1],"%lf",&templ);
    Pnoise150 = convertNoise(templ,pix2arcmin);
    sscanf(argv[2],"%lf",&templ);
    Pnoise90 = convertNoise(templ,pix2arcmin);
  }
  //frequency base unit in l space
  templ=xdim*pix2arcmin*M_PI/(180.0*60.0)/(2.0*M_PI);  // xdim=ydim

  /* CMB term
     read in from Cl.dat */
  printf(" # Reading CMB power spectrum from %s\n",clfile); fflush(stdout);
  clinput(&cl_l,&cl_p,&cl_py2,&yp1,&ypn,&cllength);

  /* Beam term */
  // 150Ghz
  printf(" # Importing Beam profile\n"); fflush(stdout);
  sprintf(beamfile,"%s.fits",BEAM150FILE);
  readimage(beamfile,&beam,&xdim,&ydim);
  if ((xdim*ydim)!=ntot){
    printf("size of beam image is wrong!\n");
    exit(0);
  }
  fbeam150=(fftw_complex*)fftw_malloc(ntot*sizeof(fftw_complex));
  fftc(beam,xdim,ydim,fbeam150);
  free(beam);
  // 90 Ghz
  sprintf(beamfile,"%s.fits",BEAM90FILE);
  readimage(beamfile,&beam,&xdim,&ydim);
  if ((xdim*ydim)!=ntot){
    printf("size of beam image is wrong!\n");
    exit(0);
  }
  fbeam90=(fftw_complex*)fftw_malloc(ntot*sizeof(fftw_complex));
  fftc(beam,xdim,ydim,fbeam90);
  free(beam);
  for(i=0;i<ntot;i++)
  {
    fbeam150[i][0]*=(float)ntot;
    fbeam150[i][1]*=(float)ntot;
    fbeam90[i][0]*=(float)ntot;
    fbeam90[i][1]*=(float)ntot;
  }
  printf(" # Generating filter profile\n"); fflush(stdout);
  // filter profile
  pfile=fopen(FILTERMSIGFILE,"w"); // write the filter variance for reference.
  // calculate relative strength of frequency (based on sub_mkbolo.cs)
  tempx=H*FREQUENCY150*1.0e9/KB/T_CMB;
  tempexpx = exp(tempx);
  j150 = tempx*(tempexpx+1.0)/(tempexpx-1.0)-4.0;
  j150 = -j150; // easy for visual check
  tempx=H*FREQUENCY90*1.0e9/KB/T_CMB;
  tempexpx = exp(tempx);
  j90 = tempx*(tempexpx+1.0)/(tempexpx-1.0)-4.0;
  j90 = -j90;

// go through all filter size k 0~11
  for(k=0;k<NFILTER;k++)
  {
    /* ####### Generate the tau profile ###### */
    /* Cluster radius */
    theta=thetasize(k)/pix2arcmin;
    fprintf(pfile,"%f ",theta);
    printf("# Calculating %f(arcmin)/%f(pixel) case\n",thetasize(k),theta);
    // build beta model in sub_tauvalue.c
    for(i=0;i<ntot;i++) {
      iy=i/xdim; if(iy>halfdim) iy = rdim-iy;
      ix=i%xdim; if(ix>halfdim) ix = rdim-ix;
      fmap.img[i]=tauvalue(ix,iy,theta);
    }

    idata=(fftw_complex*)fftw_malloc(ntot*sizeof(fftw_complex));
    fftc(fmap.img,xdim,ydim,idata);  // normalization is irrelevant  
    /* ####### Generate the filter profile ####### */
    fdata150=(fftw_complex*)fftw_malloc(ntot*sizeof(fftw_complex));
    fdata90=(fftw_complex*)fftw_malloc(ntot*sizeof(fftw_complex));
    for (i=0;i<ntot;i++)
    {
      iy=i/xdim; if(iy>halfdim) iy = rdim-iy;
      ix=i%xdim; if(ix>halfdim) ix = rdim-ix;
      l = sqrt(Squ((float)ix)+Squ((float)iy))/templ;
      if (l==0)
	Pcmb = 0;  // no cmb contribution at l=0. remove 0-mode in temperature map first in filtering.
      else if (l<cl_l[cllength])
	splint(cl_l,cl_p,cl_py2,cllength,l,&Pcmb);
      else Pcmb=cl_p[cllength];
      Pcmb = Squ(Pcmb);  // total power density at l mode
      /* build P(k) matrix */
      Pbeam150 = Squ(fbeam150[i][0])+Squ(fbeam150[i][1]);
      Pbeam90 = Squ(fbeam90[i][0])+Squ(fbeam90[i][1]);
      // 150-150
      Pa = Pnoise150 + Pbeam150*Pcmb;
      // 90-90
      Pd = Pnoise90 + Pbeam90*Pcmb;
      // 150-90
      // real part of Pc x B1B2*
      b1b2real = fbeam150[i][0]*fbeam90[i][0]+fbeam150[i][1]*fbeam90[i][1]; // real(B1B2*)
      Pc_real = b1b2real*Pcmb; 
      // image part of pc x B1B2*
      b1b2image = -fbeam150[i][0]*fbeam90[i][1]+fbeam150[i][1]*fbeam90[i][0];  // image(B1B2*)

      Pc_image = b1b2image*Pcmb;
      // |P|
      Ptot = Pnoise150*Pnoise90+Pbeam150*Pcmb*Pnoise90+Pbeam90*Pcmb*Pnoise150;
      // filter profile
      P90 = Pd*j150;
      P150 = Pa*j90;
      fdata150[i][0]=(idata[i][0]*(P90-b1b2real*Pcmb*j90)+b1b2image*Pcmb*j90*idata[i][1])/Ptot;
      fdata150[i][1]=(idata[i][1]*P90-Pcmb*j90*(b1b2real*idata[i][1]+b1b2image*idata[i][0]))/Ptot;
      fdata90[i][0]=(idata[i][0]*(P150-b1b2real*Pcmb*j150)-b1b2image*Pcmb*j150*idata[i][1])/Ptot;
      fdata90[i][1]=(idata[i][1]*P150-Pcmb*j150*(b1b2real*idata[i][1]-b1b2image*idata[i][0]))/Ptot;
      // sigma
      sigma += Squ(idata[i][0]+Squ(idata[i][1]))*(j150*P90+j90*P150-2*Pcmb*j150*j90*b1b2real)/Ptot;
    }

    for(i=0;i<ntot;i++)
    {
      fdata150[i][0]/=(sigma*ntot);
      fdata150[i][1]/=(sigma*ntot);
      fdata90[i][0]/=(sigma*ntot);
      fdata90[i][1]/=(sigma*ntot);
    }
    sigma = sqrt(1.0/sigma);  // still missing a normalization parameter
    fprintf(pfile,"%e\n",sigma);

    // save and clean
    sprintf(outputfile,"%sM150_%d",FILTERFILE,k);
    fftc2image(fdata150,xdim,ydim,outputfile);  //output image.  use .dat only.
    sprintf(outputfile,"%sM90_%d",FILTERFILE,k);
    fftc2image(fdata90,xdim,ydim,outputfile);  //output image.  use .dat only.    
    fftw_free(idata);
    fftw_free(fdata150);
    fftw_free(fdata90);
  }

  // clean up
  fclose(pfile);
  cleanmap(&fmap);
  fftw_free(fbeam150);
  fftw_free(fbeam90);
  free(cl_l);
  free(cl_p);
  free(cl_py2);
  
  return 0;
}

double convertNoise(double dT, double pix2arcmin)
{
  dT*=1.0e-6;	// K/arcmin
  dT/=pix2arcmin;		// K/arcmin  -> K/pixel
  return Squ(dT)*Squ(pix2arcmin*M_PI/(180.*60.));  // power density
}
