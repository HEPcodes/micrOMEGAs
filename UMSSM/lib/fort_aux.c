#include<stdio.h>
#include"pmodel.h"
#include"pmodel_f.h"

#include<errno.h>
#include"../../include/micromegas_aux.h"

int readlesh_(char *fname, int *SMtb, int len)
{ char * cname=malloc(len+2);
  int err;
  fName2c(fname,cname,len);
  err= readLesH(cname, *SMtb);
  free(cname);
  return err;
}

int leshinput_(char *fname, int len)
{
  char * cname=malloc(len+2);
  int err;
  fName2c(fname,cname,len);
  err=lesHinput(cname);
  free(cname);
  return err;
}

void o1contents_(int *Nch)
{
  char fname[20];
  FILE*f;

  sprintf(fname,"%d.tmptxt",getpid());
  f=fopen(fname,"w");
  o1Contents(f);
  fclose(f);
  fortreread_(Nch,fname,strlen(fname));
  unlink(fname);
}

int  hbblocks_(char * fname,int len)
{
  char * cname=malloc(len+2);
  int err;
  fName2c(fname,cname,len);
  err=HBblocks(cname);
  free(cname);
  return err;           
}

int lilithf_(char*fname,int len)
{
   char * cname=malloc(len+2);
   int err;
   fName2c(fname,cname,len);
   err= LiLithF(cname);
   free(cname);
   return err; 
}

int umssmtools_(int *PDG_LSP) {return umssmtools(*PDG_LSP);}
