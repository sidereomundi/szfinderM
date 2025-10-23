// $Id: main_szpeakfinder.cpp,v 1.5 2013/10/16 10:32:02 jiayi Exp $
// peakfinder code work for sz cluster finder
// generate the list of peaks from fits images

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include "structures.h"
#include "routines.h"
#include "routinesp.h"
#include "parameters.h"
using namespace std;

int main(int argc, char *argv[])
{
  char filename[200],fname[100],fname2[200];
  int id,xdim,ydim;
  float *img;
  vector<struct peak> phead;
  
  if(argc<2)
  {
    printf("szpeakfinder (maps name) (ID) \n using default (filtered), (szpeaks.dat)\n");
    sprintf(fname,"filtered");
  }
  else if (argc<3)
  {
    sprintf(fname,"%s",argv[1]);
    sprintf(fname2,"%s",SZPEAKLIST);
    printf("input %s ; output to default: %s\n",fname,filename);
  }
  else
  {
    sprintf(fname,"%s",argv[1]);
    sprintf(fname2,"szpeaks_%s.dat",argv[2]);
    printf("input %s ; output %s\n",fname,fname2);
  }

  for(id=0;id<NFILTER;id++)
  {
    sprintf(filename,"%s%02d.fits",fname,id);
    readimage(filename,&img,&xdim,&ydim);
    peaklist(id,&phead,img);
    free(img);
  }
  printf("total peaks %d\n",phead.size());
  // merge duplicate peaks
  peakcombine(&phead);
  printf("total peaks after merge %d\n",phead.size());

  peakoutput(&phead,fname2);
  
  return 0;
}
