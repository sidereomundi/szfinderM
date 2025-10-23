//$Id: main_szpeakcmp.cpp,v 1.5 2013/10/16 10:32:02 jiayi Exp $

// compare the true table of clusters
// input <cluster list>
//       <peak list>
// output <peak mathc list>  :: peak position, cluster position, sn, coresize, mass, redshift, cluster id, match number
// currently match all the clusters we have
// future steps: 1. add mass limit 2. match to the most massive one

#include <vector>
#include <stdlib.h>
#include "parameters.h"
#include "structures.h"
#include "routinesp.h"
#include "routines.h"
using namespace std;

int main(int argc, char *argv[])
{
  vector<struct simcluster> chead;
  vector<struct simcluster>::iterator ic;
  vector<struct peak> phead;
  vector<struct peak>::iterator ip;
  FILE *fp;
  char filename[200];
  int nmatch;

  if (argc < 2 ){
    printf("input peak file ID\n");
    exit(1);
  }
  // input cluster list
  sprintf(filename,"clusterlist_%s.dat",argv[1]);
  inputcluster(filename,&chead);
  // input peak list
  sprintf(filename,"szpeaks_%s.dat",argv[1]);
  peakinput(&phead,filename);
  
  sprintf(filename,"matchpeak_%s.dat",argv[1]);
  fp = fopen(filename,"w");
  // match peaks to clusters

  double pixelsize = (double)FIELDSIZE/NDIM * 60;
  double thetalimit;
  for (ip=phead.begin(); ip<phead.end(); ip++)
  {
    nmatch = 0;
    thetalimit = (1>ip->coresize?1:ip->coresize)/pixelsize;
    thetalimit *= thetalimit;
    for (ic=chead.begin(); ic<chead.end(); ic++)
    {
      if((Squ(ip->xpix-ic->xpos)+Squ(ip->ypix-ic->ypos) <= thetalimit))
      {
	nmatch++;
	// peak position, cluster position, sn, coresize, mass, redshift, cluster id, match number
	fprintf(fp,"%d %d %f %f %f %f %lg %f %d %d\n",
		ip->xpix,ip->ypix,ic->xpos,ic->ypos,ip->sn,ip->coresize,ic->m,ic->redshift,ic->ihal,nmatch);
      }
    }
    if (nmatch == 0) 
	    fprintf(fp,"%d %d 0 0 %f %f 0 0 0 0\n",          
		         ip->xpix,ip->ypix,ip->sn,ip->coresize);
  }
  fclose(fp);
  return 0;
}
