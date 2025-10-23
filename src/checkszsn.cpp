//$Id: checkszsn.cpp,v 1.1 2013/10/16 10:32:02 jiayi Exp $
/*
create the SZ SN - mass relation file

input SZ cluster list from simulation
output SN 
*/

#include <vector>
#include <stdlib.h>
#include "structures.h"
#include "routines.h"
#include "routinesp.h"
#include "parameters.h"
using namespace std;

/* in used routines
 */

// predefined info
#if NFILTER != 12
#error "NFILTER should be 12"
#endif
const char*  clusterfile="clusterlist_";
const float Z_UP_LIMIT = 1.6;
const double M_LOW_LIMIT = 3e13;
const float searchsize[] = {1,1,1,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3.};

void extractlocalSN(float xpix,float ypix,int xdim,float *imgs[NFILTER], float *sn){
  // extract the SN at given position
  int i;
  long ipos;

  ipos = ((long)(xpix+0.5)) + (long)(ypix+0.5)*xdim; // pixel position
  for (i=0;i<NFILTER;i++){
    sn[i] = imgs[i][ipos];
  }
}


void extractMaxSN(float xpix, float ypix,int xdim, float *imgs[NFILTER], float *sn, float pos[2])
{
  // extract the max SN around the given position
  int i,ix,iy;
  float pixelsize = FIELDSIZE/xdim*60;  // pixel size in arcmin;
  int xlow,xhigh,ylow,yhigh;
  int maxid,mx,my; // store the max
  float maxsn;
  long ipos;

  maxsn = -100;
  for(i=0;i<NFILTER;i++){
    xlow = xpix-searchsize[i]/pixelsize;  xlow = xlow>0?xlow:0;
    xhigh = xpix+searchsize[i]/pixelsize;  xhigh = xhigh>=xdim?xhigh:xdim-1;
    ylow = ypix-searchsize[i]/pixelsize;  ylow = ylow>0?ylow:0;
    yhigh = ypix+searchsize[i]/pixelsize;  yhigh = yhigh>=xdim?yhigh:xdim-1; // assume xdim=ydim;
    for (ix=xlow;ix<xhigh;ix++)
      for(iy=ylow;iy<yhigh;iy++){
	ipos = ix+iy*xdim;
	if(imgs[i][ipos]>maxsn){
	  maxsn = imgs[i][ipos];
	  mx = ix; my = iy;
	  maxid = i;
	}
      }
  }
  ipos = mx+my*xdim;
  // assign the maximum to return
  for(i=0;i<NFILTER;i++){
    sn[i]=imgs[i][ipos];
  }
  pos[0]=mx;
  pos[1]=my;
}


int main(int argc, char * argv[])
{
  // control
  char filename[300];
  FILE *fp;
  int id,count,i;

  // read-in data
  vector<struct simcluster> chead;
  vector<struct simcluster>::iterator ic;
  float *imgs[NFILTER];
  int xdim,ydim;

  // output-data
  float sn[NFILTER];
  float pos[2];

  if(argc<4){
    printf("!! checkszsn <sim_id> <filtered_name_format> <outputname>!!\n");
    exit(1);
  }

  printf("# loading data %s\n",argv[1]);
  sprintf(filename,"%s%s.dat",clusterfile,argv[1]);
  inputcluster(filename,&chead);
  
  for(id=0;id<NFILTER;id++){
    sprintf(filename,argv[2],id);
    readimage(filename,&(imgs[id]),&xdim,&ydim);
  }
  
  printf("# looping clusters\n");
  
  fp = fopen(argv[3],"w");
  fprintf(fp,"# ihal, Mvir[Msun/h], redshift, xpos, ypos, xdect, ydect, SN_max, SN_local\n");
  count = 0;
  for(ic = chead.begin();ic<chead.end();ic++){
    if (ic->m < M_LOW_LIMIT || ic->redshift > Z_UP_LIMIT)
      continue;
    count ++;
    // load local SN
    extractMaxSN(ic->xpos,ic->ypos,xdim,imgs,sn,pos);
    fprintf(fp,"%d %e %f %f %f %f %f ",ic->ihal,ic->m,ic->redshift,ic->xpos,ic->ypos,pos[0],pos[1]);
    for(i=0;i<NFILTER;i++)
      fprintf(fp,"%f ",sn[i]);
    extractlocalSN(ic->xpos,ic->ypos,xdim,imgs,sn);
    for(i=0;i<NFILTER;i++)
      fprintf(fp,"%f ",sn[i]);
    fprintf(fp,"\n");
  }
  
  printf("# finishing analysis of %d clusters\n",count);
  fclose(fp);
  return 0;
}

