//$Id: sub_cluster.cpp,v 1.4 2013/10/16 10:32:02 jiayi Exp $

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include "parameters.h"
#include "structures.h"
using namespace std;

void loadgeometry(char *dir, struct geometry *geo)
{
  char trash[300],line[500];
  FILE *fp;
  int count=0;
  float trashf;
  sprintf(trash,"%s/wmap.geometry.dat",dir);
  fp = fopen(trash,"r");
  if (fp == NULL)
  {
    printf("error in open %s\n",trash);
    exit(1);
  }
  fgets(line,500,fp);
  fgets(line,500,fp);
  fgets(line,500,fp);  // omit first three lines
  while(fgets(line,500,fp)!=NULL)
  {
    if (count >= NZSLICE)
    {
      printf(" out of simulatin redshift max number \n");
      exit(2);
    }
    sscanf(line,"%d %f %lf %s",(geo->isnap+count),&trashf,(geo->zslice+count),trash);
    count ++;
  }
  fclose(fp);
}

float findredshift(int isnap, struct geometry *geo)
{
  int i;
  for (i=0;i<NZSLICE;i++)
    if(isnap == geo->isnap[i]) return geo->zslice[i];
  printf(" can't find the redshift in simulation %d\n",isnap);
  exit(3);
}

void inputcluster(char *filename,vector<struct simcluster> *chead)
{
  FILE *fp;
  char line[500];
  struct simcluster temp;
  fp = fopen(filename,"r");
  if (fp == NULL)
  {
    printf(" failed to open %s \n",filename);
    exit(1);
  }
  while(fgets(line,500,fp)!=NULL)
  {
    if(line[0]=='#') // omit comment line!
      continue;
    sscanf(line,"%d %f %f %lg %f",&(temp.ihal),&(temp.xpos),&(temp.ypos),&(temp.m),&(temp.redshift));
    temp.xpos *= NDIM;
    temp.ypos *= NDIM;
    chead->push_back(temp);
  }
  fclose(fp);
}
