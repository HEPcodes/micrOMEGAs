#include<stdio.h>
#include<stdlib.h>
#include"micromegas.h"
#include"micromegas_aux.h"
#include"ic22.h"

double IC22nuAr_(double E)
{ if(E<=36) return 0;
  else if(E<72.7108)  return  4.403e-06 *1E-6;
  else if(E<105.737)  return  9.525e-05 *1E-6;
  else if(E<153.765)  return  0.0007292 *1E-6;
  else if(E<223.607)  return  0.002949  *1E-6;
  else if(E<325.172 ) return  0.009057  *1E-6;
  else if(E<472.871 ) return  0.02201   *1E-6;
  else if(E<687.656 ) return  0.04820   *1E-6;
  else /*if(E<1000 ) */return 0.1002    *1E-6;  
}

double IC22nuAr(double E)
{ 
  double lnE_[8]={  4.099, 4.474, 4.848, 5.223, 5.597, 5.972, 6.346, 6.720};
  double lnAr[8]={-12.333,-9.259,-7.224,-5.826,-4.704,-3.816,-3.032,-2.300};

  if(E<36) return 0;
  if(E>1E4)E=1.E4;
  return exp(polint3(log(E),8,lnE_,lnAr))*1E-6; 
}


double IC22sigma(double E)
{ double  lnE_[8]= {  4.099, 4.474, 4.848, 5.223, 5.597, 5.972, 6.346, 6.720};
   double  fiTab[8]={  3.69,  3.20,  2.96 ,  2.67, 2.47,  2.32,  2.24,  2.05}; 
//    double  fiTab[8]={  4.5,  3.9,  3.5 ,  3.2, 3.0, 2.8,  2.7,  2.5};
  return  (M_PI/180)*polint3(log(E),8,lnE_,fiTab);
}
  
double IC22nuBarAr_(double E)
{ if(E<50) return 0;
  else if(E<72.7108)  return  3.154E-06 *1E-6;
  else if(E<105.737)  return  8.150E-05 *1E-6;
  else if(E<153.765)  return  0.0005786 *1E-6;
  else if(E<223.607)  return  0.002301 *1E-6;
  else if(E<325.172 ) return  0.006607 *1E-6;
  else if(E<472.871 ) return  0.01649 *1E-6;
  else if(E<687.656 ) return  0.03537 *1E-6;
  else /*if(E<1000 ) */return 0.06991*1E-6; 
}

double IC22nuBarAr(double E)
{ 
   double lnE_[8]={  4.099, 4.474, 4.848, 5.223, 5.597, 5.972, 6.346, 6.720};
   double lnAr[8]={-12.667,-9.415,-7.455,-6.074,-5.020,-4.105,-3.342,-2.660};
   if(E<50) return 0;
   if(E>1E4)E=1.E4;
   return exp(polint3(log(E),8,lnE_,lnAr))*1E-6;      
}

double IC22BGdCos(double cs)
{
  double fi_bg[25]={2.034E+00,6.017E+00,9.895E+00,1.370E+01,1.725E+01,2.086E+01,2.431E+01,2.698E+01,2.937E+01,3.160E+01,
                    3.386E+01,3.595E+01,3.821E+01,4.022E+01,4.188E+01,4.368E+01,4.547E+01,4.727E+01,4.902E+01,5.080E+01,
                    5.227E+01,5.371E+01,5.399E+01,5.184E+01,4.927E+01};


// fi_bg angular distribution of backgraund simulation. From file  DarkSUSY/IC_data/BG_distributions_IC22.dat
                     
  int i;
  for(i=0;i<25;i++) fi_bg[i]/=cos(i*M_PI/180)-cos((i+1)*M_PI/180.);
  double cs_arr[25];
  for(i=0;i<25;i++) cs_arr[i]=cos((i+0.5)*M_PI/180.);
  return  polint2(cs,25,cs_arr,fi_bg);
}


IC22chanStr*IC22chan=NULL;

int IC22histRead(void)
{  FILE*F;
   char fname[300];
   int i;
   double E1,E2;

   if(IC22chan) return 0;
   
   sprintf(fname,"%s/sources/data_nu/ic22hist.dat",micrO);
        
   F=fopen(fname,"r");
   if(!F) return 1;
   

   IC22chan=malloc(21*sizeof(IC22chanStr));
         
   for(i=0;2==fscanf(F," %*s %lf %lf",&E1,&E2); i++)
   {  int j; 
      IC22chan[i].E1=pow(10,E1);  IC22chan[i].E2=pow(10,E2);
      for(j=0;j<17;j++) fscanf(F," %lf", IC22chan[i].prob+j);
   }   
   fclose(F);                                                    
}

static   double pp(double cs,  int n, double a)
   { double s= IC22BGdCos(cs)*IC22chan[0].prob[n];
     int i;
     for(i=1;i<=20;i++) if(IC22chan[i].n>0)
     { double s2=IC22chan[i].s2;
       s+= a*IC22chan[i].n*exp((cs-1)/s2)/s2*IC22chan[i].prob[n];
     } 
     return s;  
   }
                 
static   double L[16];
static    double P(double x)
   { double xx[16]={0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,1.1,1.2,1.3,1.4,1.5};
     return exp(polint3(x,16,xx,L));
   }
static   double LLL(double x)
       { double xx[16]={0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,1.1,1.2,1.3,1.4,1.5};
            return polint3(x,16,xx,L);
               }
   
   double pvalIC22(double * nu, double*NU, double *LL)
   {  int i;
      double cs,cs_;
      typedef struct{double cs; int n;} eventStr;
      static  eventStr*  events=NULL;
      static int Nev;
      double fiMax=8;   
      if(fiMax>10) { printf("Maximum angel reset 10 degrees\n"); fiMax=10;}
      if(!events)
      {  char fname[300];
         double dfi;
         int nch;
         events=malloc(180*sizeof(eventStr));
         sprintf(fname, "%s/sources/data_nu/IC22_events_25.dat",micrO); 
         FILE*F=fopen(fname,"r");
         cs_=cos(10./180.*M_PI);    
         for(Nev=0; fscanf(F," %lf %lf %d",&cs,&dfi, &nch)==3; ) if(cs>=cs_)
         { if(Nev==180) break;
           events[Nev].cs=cs;
           events[Nev].n=nch;
           Nev++; 
         }
         fclose(F);
         IC22histRead();  
      }
      cs_=cos(fiMax/180*M_PI);

      double  nu_[NZ], NU_[NZ];
      for(i=0;i<NZ;i++) {nu_[i]=nu[i];NU_[i]=NU[i];} 
      spectrMult(nu_,IC22nuAr);
      spectrMult(NU_,IC22nuBarAr); 
      addSpectrum(nu_,NU_);
      for(i=0;i<NZ;i++) NU_[i]=nu_[i];
      spectrMult(NU_,IC22sigma); spectrMult(NU_,IC22sigma);
      
      double M=nu_[0];
      for(i=1; i<=20;i++)
      {  
         double E1=IC22chan[i].E1,E2=IC22chan[i].E2;
         if(E1<M)
         {
           if(E2>M) E2=M;
           IC22chan[i].n=spectrInt(E1,E2,nu_);
           if(IC22chan[i].n<=0) IC22chan[i].s2=0; else
           IC22chan[i].s2=spectrInt(E1,E2,NU_)/IC22chan[i].n;
         } else {IC22chan[i].n=0; IC22chan[i].s2=0;} 
         IC22chan[i].n*=104/365./1.2; 
      }
      double Nbg=simpson(IC22BGdCos,cs_,1,1E-3);
      double Ns=0;
      for(i=1;i<20;i++) if(IC22chan[i].n>0) Ns+=IC22chan[i].n*(1-exp((cs_-1)/IC22chan[i].s2));
      double nData=0;
      for(i=0;i<Nev;i++) if(events[i].cs>cs_) nData++;
      int k;
      for(k=0;k<16;k++)
      { double a= 0.1*k;
        L[k]=nData*log(Nbg+a*Ns)-Nbg-a*Ns;
        for(i=0;i<Nev;i++)if(events[i].cs>cs_) L[k]+=log(pp(events[i].cs,events[i].n-10,a)/(Nbg+a*Ns)); 
      }
      
      for(k=1;k<16;k++) L[k]-=L[0];
      L[0]=0;
      
     if(L[15]>L[0]-0.1) return 1;
      double dI=exp(L[15])/(2*(L[10]-L[15]));
//displayFunc(P, 0,1.5,"P"); 
//for(int i=0;i<16;i++) printf("L[%d]=%E %E\n",i,L[i], LLL(0.1*i));   
//displayFunc(LLL,0,1.5,"LL");        
      double int1=simpson(P,0.,1 ,1E-3);
      double int2=simpson(P,1.,1.5,1E-2);
      if(LL) *LL=-log(P(1)/P(0));
      return (int2+dI)/(int1+int2+dI);
   }
#ifdef TEST_ANOMALY   
   double pvalIC22_randon(double * nu, double*NU, double fiMax, int in)
   {  int i;
      double cs,cs_;
      typedef struct{double cs; int n;} eventStr;
      static  eventStr*  events=NULL;
      static int Nev;
      
      if(!events)
      {  char fname[300];
         double dfi;
         int nch;
         events=malloc(1000*sizeof(eventStr));
         sprintf(fname, "%s/sources/data_nu/IC22_events_25.dat",micrO); 
         FILE*F=fopen(fname,"r");
         cs_=cos(25./180.*M_PI);    
         for(Nev=0; fscanf(F," %lf %lf %d",&cs,&dfi, &nch)==3; ) if(cs>=cs_)
         { if(Nev==1000) break;
           events[Nev].cs=cs;
           events[Nev].n=nch;
           Nev++; 
         }
         fclose(F);
         IC22histRead();  
      }
      cs_=cos(10./180.*M_PI);

if(in)
{
//printf("nData=%d\n",Nev);
for(;;)
{double  c=1;
  Nev=130+drand48()*100;
  
  for(i=182;i<=Nev;i++) c*=182./i;
  for(i=182;i>Nev;i--) c/=182./i;
  if(c>drand48()) break; 
}

printf("Nev=%d\n",Nev);
for(i=0;i<Nev;i++)
{ double cs;
  for(;;)
  { cs=cs_+drand48()*(1-cs_);
    if(IC22BGdCos(cs)/IC22BGdCos(1)>drand48()) break; 
  }
//printf("i=ok\n", i);  
  events[i].cs=cs;
}  
//printf("fill ok\n");
}
      double  nu_[NZ], NU_[NZ];
      for(i=0;i<NZ;i++) {nu_[i]=nu[i];NU_[i]=NU[i];} 
      spectrMult(nu_,IC22nuAr);
      spectrMult(NU_,IC22nuBarAr); 
      addSpectrum(nu_,NU_);
      for(i=0;i<NZ;i++) NU_[i]=nu_[i];
      spectrMult(NU_,IC22sigma); spectrMult(NU_,IC22sigma);
      
      double M=nu_[0];
      for(i=1; i<=20;i++)
      {  
         double E1=IC22chan[i].E1,E2=IC22chan[i].E2;
         if(E1<M)
         {
           if(E2>M) E2=M;
           IC22chan[i].n=spectrInt(E1,E2,nu_);
           if(IC22chan[i].n<=0) IC22chan[i].s2=0; else
           IC22chan[i].s2=spectrInt(E1,E2,NU_)/IC22chan[i].n;
         } else {IC22chan[i].n=0; IC22chan[i].s2=0;} 
         IC22chan[i].n*=104/365.*0.8; 
      }
      double Nbg=simpson(IC22BGdCos,cs_,1,1E-3);
      double Ns=0;
      for(i=1;i<20;i++) if(IC22chan[i].n>0) Ns+=IC22chan[i].n*(1-exp((cs_-1)/IC22chan[i].s2));
      double nData=0;
      for(i=0;i<Nev;i++) if(events[i].cs>cs_) nData++;
      int k;
      for(k=0;k<16;k++)
      { double a= 0.1*k;
        L[k]=nData*log(Nbg+a*Ns)-Nbg-a*Ns;
        for(i=0;i<Nev;i++)if(events[i].cs>cs_) L[k]+=log(pp(events[i].cs,events[i].n-10,a)/(Nbg+a*Ns)); 
      }
      for(k=1;k<15;k++) L[k]-=L[0];
      L[0]=0;
//displayFunc(P,0,1.5,"P");    
//displayFunc(LL,0,1.5,"LL");        
      double int1=simpson(P,0.,1. ,1E-3);
      double int2=simpson(P,1.,1.5,1E-2);
      return int2/(int1+int2);
   }
#endif

double  fluxFactorIC22(double pval,double *NU,double*NUbar)
{ 
  double nu[NZ],nu_[NZ];
  for(int i=0;i<NZ;i++) { nu[i]=NU[i]; nu_[i]=NUbar[i];}
  double f=1; 
  for(;;)
  {
     double p0= pvalIC22(nu,nu_, NULL);
     double x;
     if(p0==0) x=0.1; else  x=log(pval)/log(p0);
     f*=x;
     for(int i=1;i<NZ;i++) { nu[i]*=x; nu_[i]*=x;}
     if(fabs(x-1)<0.01) break;
  }
  return f;
}


double ic22nuar_(double* E)  { return IC22nuAr(*E);  }
double ic22nubarar_(double*E){ return IC22nuBarAr(*E);}

double ic22bgdcos_(double*cs){ return IC22BGdCos(*cs); }
double ic22sigma_(double*E)  { return IC22sigma(*E); }
                                     

double pvalic22_(double*nu,double*NU,double*L)
                             { return  pvalIC22(nu,NU,L); }

double fluxfactoric22_( double* pval,double *NU,double*NUbar)
                             { return fluxFactorIC22(*pval,NU,NUbar);} 
