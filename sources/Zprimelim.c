
#include<math.h>
#include"micromegas.h"
#include"micromegas_aux.h"

// Routines to derive ATLAS Z' limits with 3.2 fb-1 at sqrt(s) = 13 TeV

// Two parts are purely UMSSM-dependent, the rest of the code can be used to apply in any model with a Z'
// Nevertheless caution between the difference in the computation of cross section with hCollider+specific PDF set and the theoretical curve derived by the experimental collaboration, a modified version of the second UMSSM-dependent part might be needed


int Zprimelimits(void)
{
    
 if(!pdg2name(32)){printf("ERROR - there is no Zprime with PDG code = 32 in this model.\n");return 1;}
    
else
{    
double MZp=pMass(pdg2name(32));
char * Zpname=pdg2name(32);
printf("%s mass : %5.3f GeV\n",Zpname,MZp);
printf("==================================\n");
printf("==== Limits on the %s boson : ====\n",Zpname);
printf("==================================\n");
printf("Using Z' limits from dilepton final state with 3.2 fb-1 at sqrt(s) = 13 TeV\n");
printf("From ATLAS-CONF-2015-070 (http://cds.cern.ch/record/2114842)\n");
printf("==================================\n");
if(MZp<500.) {printf("ERROR - this routine cannot test %s mass below 500 GeV\n",Zpname);return 1;}

else
{

// First UMSSM-dependent step when MZ2 is large and tE6 not in [-0.6864;0.6595] to avoid useless tests of very heavy Zprime :
// Derived lower bound on the Z2 mass and the associated tE6 angle for which the ATLAS limits exclude a UMSSM point in the case Br(Z2 -> SM)=1.
double AngleE6[44]={+1.5708,+1.5259,+1.4810,+1.4362,+1.3913,+1.3464,+1.3015,+1.2566,+1.2118,+1.1669,+1.1220,+1.0771,+1.0322,+0.9874,+0.9425,+0.8976,+0.8527,+0.8078,+0.7630,+0.7181,+0.6595,0.,-0.6864,-0.7181,-0.7630,-0.8078,-0.8527,-0.8976,-0.9425,-0.9874,-1.0322,-1.0771,-1.1220,-1.1669,-1.2118,-1.2566,-1.2925,-1.3015,-1.3464,-1.3913,-1.4362,-1.4810,-1.5259,-1.5708};
double ZpE6[44]={2785.3,2793.8,2802.0,2809.7,2817.0,2823.8,2830.2,2836.3,2842.3,2848.2,2854.2,2860.4,2866.8,2873.5,2880.6,2888.1,2896.0,2904.5,2913.5,2922.8,2935.4,3080.,2934.7,2922.7,2904.3,2886.1,2867.8,2849.6,2831.7,2814.6,2798.5,2783.7,2771.1,2760.9,2753.7,2749.6,2748.5,2748.6,2750.4,2754.8,2761.0,2768.5,2776.8,2785.3};
double tE6,Mlim=0;
int E6model=findVal("tE6",&tE6);
if(E6model!=0)printf("Pure UMSSM-dependent test skipped\n");
else
{
if(tE6>=0.6595 || tE6<=-0.6864) {Mlim=polint3(tE6,44,AngleE6,ZpE6);printf("Interpolated MZ2 limit in the case Br(Z2 -> SM)=1, for tE6=%+5.4f rad : %5.3f GeV\n",tE6,Mlim);}
}


if(MZp > Mlim)
{
printf("-> below the %s mass chosen\n",Zpname);
printf("==================================\n");
}

else
{
if(Mlim==0.) printf("-> First derived lower bound not applicable here\n");
else printf("-> above the %s mass chosen\n",Zpname);
printf("Check whether this point is safe after including BSM branching fractions\n");

double BrZpll; 
double width,br;
txtList L;
width=pWidth(Zpname,&L);
printf("\n%s :   total width=%.2E[GeV]\n",Zpname,width);
printf("leptonic branching of %s : %4.2f %%\n",Zpname,100*(findBr(L,"e,E") + findBr(L,"m,M")));
BrZpll=(findBr(L,"e,E")+findBr(L,"m,M"))/2.;


double cs,Pcm=6500,pTmin=0,Qren=pTmin,Qfact=MZp;
int nf=3;
printf("Computation of pp -> %s +jet(pt>%.2E GeV)  at %.2E GeV :\n",Zpname,pTmin,Pcm); 

  char oldPDF[50];
  strcpy(oldPDF,pdfName);
  setPDT("cteq6l1");  
cs=hCollider(Pcm,1,nf,Qren, Qfact,Zpname,NULL,pTmin,0);
printf("cs(pp->%s)=%.2E[fb]\n",Zpname,cs*pow(10,3));
printf("cs(pp->%s)*Br(%s->ll)=%.3E[fb]\n",Zpname,Zpname,cs*pow(10,3)*BrZpll);
restorePDF(oldPDF);

// Second UMSSM-dependent step to be compatible with the sigma*Br(Z' -> l+l-)_th shown in ATLAS-CONF-2015-070 :
// Data from the theoretical curve sigma*Br(Z' -> l+l-) shown in ATLAS-CONF-2015-070 (see Fig 3c) for the model psi :
double ZpPsi[43]={.5,.55,.6,.65,.7,.75,.8,.85,.9,.95,1.,1.05,1.15,1.2,1.25,1.3,1.35,1.4,1.45,1.5,1.55,1.6,1.65,1.7,1.75,1.8,1.85,1.9,1.95,2.,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.,3.3,3.5,3.75};//Z_psi mass in TeV
double sigBPsiATLAS[43]={2333.,1645.,1172.,875.7,648.2,494.,391.4,301.2,243.3,194.7,155.7,128.3,86.18,72.37,60.19,51.04,42.87,36.7,30.53,26.14,22.38,18.98,16.57,14.19,12.39,10.71,9.351,8.085,7.128,6.163,4.743,3.651,2.837,2.205,1.730,1.358,1.076,.853,.669,.541,.272,.177,.104};//sigma*Br(Z' -> l+l-)_th
// Array of ratios between sigma*Br(Z' -> l+l-)_th and the sigma*Br(Z2 -> l+l-) obtained with hCollider using set 25 of NNPDF23_lo_as_0130.LHgrid for each Z' mass given in the array ZpPsi, all this in the UMSSM case Br(Z2 -> SM)=1 :
double RatiosigBPsi[43]={1.436,1.442,1.421,1.439,1.414,1.404,1.433,1.398,1.421,1.412,1.394,1.405,1.385,1.394,1.387,1.398,1.387,1.402,1.370,1.370,1.369,1.352,1.368,1.354,1.367,1.361,1.364,1.352,1.367,1.351,1.349,1.345,1.342,1.338,1.336,1.335,1.338,1.340,1.324,1.343,1.326,1.341,1.360};
// Get the rescaled sigma*Br(Z' -> l+l-) using a cubic interpolation of the array of ratios above : 
double sBth=cs*pow(10,3)*BrZpll*polint3(MZp*pow(10,-3),43,ZpPsi,RatiosigBPsi);
printf("Rescaled cs(pp->%s)*Br(%s->ll) with ATLAS data for (tE6,MZ2)=(%+5.4f rad, %5.3f GeV) = %.3E[fb]\n",Zpname,Zpname,tE6,MZp,sBth);

// Data from the observed limit curve sigma*Br(Z' -> l+l-) shown in ATLAS-CONF-2015-070 (see Fig 3c) :
double Zpmass[36]={.5,.55,.6,.65,.7,.75,.8,.85,.9,.95,1.,1.05,1.15,1.2,1.25,1.3,1.35,1.4,1.45,1.5,1.55,1.6,1.65,1.7,1.75,1.8,1.85,1.9,1.95,2.,2.2,2.7,2.8,3.,3.3,4.};//Z' mass in TeV
double sigBexp[36]={16.83,12.22,16.51,12.22,5.961,5.624,4.441,2.77,4.205,4.68,3.083,2.35,2.854,2.72,2.883,2.668,3.24,2.857,2.399,1.659,1.289,1.104,1.012,1.002,1.205,1.278,1.241,1.053,.946,0.91,
0.885,0.861,0.887,0.913,0.941,1.117};//sigma*Br(Z' -> l+l-)_exp
double sBexp=polint3(MZp*pow(10,-3),36,Zpmass,sigBexp);
printf("Observed limit cs(pp->Z')*Br(Z'->ll) for MZprime=%5.3f GeV : %.3E[fb]\n",MZp,sBexp);


if(sBexp >= sBth)
{
printf("==================================\n");
}
else
{
printf("%s too light, excluded by LHC limits on Z'\n",Zpname);
printf("==================================\n");
return 1;
}
}

return 0;
}
}
}

extern int zprimelimits_(void); //Fortran 
int  zprimelimits_(void) { return  Zprimelimits();}