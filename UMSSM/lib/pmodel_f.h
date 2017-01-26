#ifndef  __MSSM_F_
#define  __MSSM_F_


/*===============================================
      Les Houches accord interface
 =================================================*/
extern int readlesh_(char *fname , int*SM,  int len );
/* 
     integer function readLesH(fname,SM)
     character*(*) fname
     int SM
*/  

void o1contents_(int *N);
/*
  subroutine  o1Contents(int *N )
*/
    
#endif

