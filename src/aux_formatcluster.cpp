//$Id: aux_formatcluster.cpp,v 1.8 2013/02/08 10:34:56 jiayi Exp $
// merge different cluster file into one.

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include "parameters.h"
#include "structures.h"
#include "routinesp.h"
using namespace std;

/*   work for old cluster.dat file 
int main(int argc, char*argv[])
{
  int i;
  vector<struct simcluster> clusterlist;
  vector<struct simcluster>::iterator ic;
  struct simcluster temp;
  struct geometry geo;
  FILE *filer,*filew;
  char fname[300],line[400];
  float tempz;
  if (argc != 2)
  {
    printf("reformat cluster catalogs <path> (no '/')\n");
    exit(1);
  }
  loadgeometry(argv[1],&geo);
  for (i=0;i<150;i++)
  {
    sprintf(fname,"%s/wmap.%03d.cluster.dat",argv[1],i);
    filer=fopen(fname,"r");
    if (filer != NULL)
    {
      tempz = findredshift(i,&geo);
      fgets(line,400,filer); // headline
      while(fgets(line,400,filer)!=NULL)
      {
	sscanf(line,"%d %f %f %lf %lf",&(temp.ihal),&(temp.xpos),&(temp.ypos),&(temp.rad),&(temp.m));
	temp.redshift = tempz;
	clusterlist.push_back(temp);
      }
      fclose(filer);
    }

  }
  filew = fopen(CLUSTERLIST,"w");
  for(ic=clusterlist.begin();ic<clusterlist.end();ic++)
    fprintf(filew,"%d %f %f %lg %f\n",ic->ihal,ic->xpos,ic->ypos,ic->m,ic->redshift);
  fclose(filew);
  
  return 0;
}

*/

 // work for new cluster.dat and cluster.sim.dat 
int main(int argc, char*argv[])
{
  int i;
  vector<struct simcluster> clusterlist;
  vector<struct simcluster>::iterator ic;
  struct simcluster temp;
  FILE *filer1,*filew,*filer2;
  char fname[300],line[400],line2[400];
  float tempz;
  double dtrash;
  float ftrash;
  int tempi;
  if (argc < 2)
  {
    printf("reformat cluster catalogs <path,(no '/' needed)> (optional output file_ID)\n");
    exit(1);
  }
  
  for (i=0;i<150;i++)
  {
    sprintf(fname,"%s/wmap.%03d.cluster.dat",argv[1],i);
    filer1=fopen(fname,"r");
    sprintf(fname,"%s/wmap.%03d.cluster.sim.dat",argv[1],i);
    filer2=fopen(fname,"r");   
    if ((filer1 != NULL) && (filer2 != NULL))
    {
      fgets(line,400,filer1); // headline
      fgets(line2,400,filer2);
      while(fgets(line,400,filer1)!=NULL)
      {
	if(fgets(line2,400,filer2) == NULL)
	{
	  printf(" cluster.dat and cluster.sim.dat are not matched to each other !!\n");
	  exit(3);
	}
	sscanf(line,"%d %f %f %f %lf %lf",&(temp.ihal),&(temp.xpos),&(temp.ypos),&ftrash,&(temp.rad),&(temp.m));
	sscanf(line2,"%d %f %f %f %f %f %f %lf %f %f %f %d",&tempi,&ftrash,&ftrash,&ftrash,&ftrash,&ftrash,&ftrash,
	       &dtrash,&ftrash,&(temp.redshift),&tempz,&(temp.flag));
	if(tempi != temp.ihal)
	{
	  printf(" order in cluster.dat and cluster.sim.dat are not matched!!\n");
	  printf(" file %d,  ihal: %d %d\n",i,temp.ihal,tempi);
	  exit(4);
	}
	clusterlist.push_back(temp);
      }
      fclose(filer1);
      fclose(filer2);
    }
  }
  if (argc ==2)
  {
    printf("using default filename: %s\n",CLUSTERLIST);
    filew = fopen(CLUSTERLIST,"w");
  }
  else
  {
    sprintf(fname,"clusterlist_%s.dat",argv[2]);
    printf("outputing to %s\n",fname);
    filew = fopen(fname,"w");
  }
  
  for(ic=clusterlist.begin();ic<clusterlist.end();ic++)
    fprintf(filew,"%d %f %f %lg %f\n",ic->ihal,ic->xpos,ic->ypos,ic->m,ic->redshift);
  fclose(filew);
  
  return 0;
}
