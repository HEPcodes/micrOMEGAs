#ifndef  __IC22__
#define  __IC22__

#ifdef __cplusplus
extern "C" {
#endif 

extern double IC22nuAr_(double E);
// muon neutrino effective arrea in [km^2]. Step fuction realisation 
// based on DS/shared/DarkSusy/IC_data/nuEffArea_IC22.dat

extern double IC22nuBarAr_(double E);
// muon anti-neutrino effective arrea in [km^2]. Step fuction realisation 
// based on DS/shared/DarkSusy/IC_data/nuEffArea_IC22.dat 

extern double IC22nuAr(double E);    // smoothed  IC22nuAr_
extern double IC22nuBarAr(double E); // smoothed  IC22nuBarAr_

extern double IC22BGdCos(double cs);  // background  distribution for cos(Sun,neutrino)  
                                     // based on  DS/shared/DarkSusy/IC_data/BG_distributions_IC22.dat
extern double IC22sigma(double E);     // angle(degree) resolution as a function of energy 
                                     //  based on DS/shared/DarkSusy/IC_data/nuEffArea_IC22.dat
                                     
typedef   struct { double E1; double E2; double s2; double n; double prob[17];}  IC22chanStr;                                    

extern IC22chanStr*IC22chan;

extern int IC22histRead(void);  // reads file data_nu/ic22hist.dat to fill 21 elements of IC22nChan.
                                // zero^th one corresponds to background. Next 20 present  channel distribution for different 
                                // energy regions.      

extern double pvalIC22(double * nu, double*NU, double*L);

extern double  fluxFactorIC22(double pval,double *NU,double*NUbar);

#ifdef __cplusplus
}
#endif 


#endif
