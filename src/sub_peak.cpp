//$Id: sub_peak.cpp,v 1.5 2013/10/16 10:32:02 jiayi Exp $
// subroutines for peak finder

#include <vector>
#include <stdlib.h>
#include "parameters.h"
#include "structures.h"
#include "routines.h"


using namespace std;

// check whether the pixel is a peak
bool checkpeak(long pos, long xdim, float *image)
{
  if((image[pos]>image[pos+1])&&(image[pos]>image[pos-1])
     &&(image[pos]>image[pos-xdim])&&(image[pos]>image[pos+xdim])
     &&(image[pos]>SIGMA_THRESHOD)){
    return true;
  }
  return false;
}

// build peak list
void peaklist(int coresize,vector<struct peak> *phead,float *img)
{
  int ix,iy,lowbound,highbound;
  long pix,nsidepix;
  struct peak temp;
  // set boundary
  if (EXTENDEDDIM == 0)
  {
    lowbound = 3;
    highbound = NDIM+EXTENDEDDIM-3;
  }
  else
  {
    lowbound = EXTENDEDDIM;
    highbound = NDIM+EXTENDEDDIM;
  }
  nsidepix = NDIM + 2*EXTENDEDDIM;
  for (ix = lowbound; ix<highbound; ix++)
    for (iy = lowbound; iy<highbound; iy++)
    {
      pix = ix+iy*nsidepix;
      if (checkpeak(pix,nsidepix,img))
      {
	temp.xpix=ix-EXTENDEDDIM;
	temp.ypix=iy-EXTENDEDDIM;
	temp.sn=img[pix];
	temp.coresize=thetasize(coresize);
	phead->push_back(temp);
      }
    }
}

 // remove duplicate peaks
void peakcombine(vector<struct peak> *phead)
{
  vector<struct peak>::iterator i1,i2,itemp;
  double pixelsize;  // [arcmin]
  double thetalimit;  // radius to search unique peak
  pixelsize = (double)FIELDSIZE/NDIM * 60;
  for(i1 = phead->begin(); i1<phead->end(); i1++)
    for( i2 = i1+1; i2<phead->end(); i2++)
    {
      thetalimit = i1->coresize>i2->coresize ? i1->coresize : i2->coresize;
      thetalimit = 1. > thetalimit ? 1. : thetalimit;   // select the max coresize of 1 arcmin
      thetalimit /= pixelsize;   // convert to pixel scale
      thetalimit *= thetalimit;  // direct comparable
      if ((Squ(i1->xpix-i2->xpix)+Squ(i1->ypix-i2->ypix)) < thetalimit)
      {
	if(i1->sn > i2->sn)
	{
	  itemp = i2-1;
	  phead->erase(i2);
	  i2 = itemp;
	}
	else
	{
	  itemp = i1-1;
	  phead->erase(i1);
	  i1 = itemp;
	}
      }
    }
}

void peakoutput(vector<struct peak> *phead,char*fname)
{
  vector <struct peak>::iterator ipeak;
  FILE *fp;
  fp = fopen(fname,"w");
  for (ipeak=phead->begin(); ipeak<phead->end(); ipeak++)
  {
    fprintf(fp,"%d %d %f %f\n",ipeak->xpix,ipeak->ypix,ipeak->sn,ipeak->coresize);
  }
  fclose(fp);
}

void peakinput(vector<struct peak> *phead,char*fname)
{
  struct peak temp;
  FILE *fp;
  char line[500];
  fp = fopen(fname,"r");
  if(fp == NULL)
  {
    printf(" failed to open %s\n",fname);
    exit(1);
  }
  while(fgets(line,500,fp)!=NULL)
  {
    sscanf(line,"%d %d %f %f",&(temp.xpix),&(temp.ypix),&(temp.sn),&(temp.coresize));
    phead->push_back(temp);
  }
  fclose(fp);
}
