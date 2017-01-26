      SUBROUTINE CONSTSUSYPART_CPV(PAR,PROB)

*   Subroutine to check experimental constraints on SUSY particles
*
*   The required data files and numbers (inv. Z width, lower bounds on
*   sparticle masses) are transferred via COMMON/LEP/...

*      PROB(1) =/= 0   chargino too light
*      PROB(2) =/= 0   excluded by Z -> neutralinos
*      PROB(20) =/= 0  excluded by stop -> b l sneutrino
*      PROB(21) =/= 0  excluded by stop -> neutralino c
*      PROB(22) =/= 0  excluded by sbottom -> neutralino b
*      PROB(23) =/= 0  squark/gluino too light
*      PROB(24) =/= 0  selectron/smuon too light
*      PROB(25) =/= 0  stau too light


      IMPLICIT NONE

      INTEGER I,J,M
      INTEGER NhZind,NhZbb,NhZll,NhZinv,NhZjj,NhZgg
      INTEGER NhA4b,NhA4tau,NhA2b2tau,NhA2tau2b
      INTEGER NAAA6b,NAAA6tau,NAAZ4b,NAAZ4tau,NAAZ2b2tau
      INTEGER Ncccc02,Ncccc04,Ncccc05,Ncccc06,Ncccc08,Ncccc1
      INTEGER Nccgg02,Nccgg04,Nccgg05,Nccgg06,Nccgg08,Nccgg1
      INTEGER Ncctt02,Ncctt04,Ncctt05,Ncctt06,Ncctt08,Ncctt1
      INTEGER Ngggg02,Ngggg04,Ngggg05,Ngggg06,Ngggg08,Ngggg1
      INTEGER Nttgg02,Nttgg04,Nttgg05,Nttgg06,Nttgg08,Nttgg1
      INTEGER Ntttt02,Ntttt04,Ntttt05,Ntttt06,Ntttt08,Ntttt1
      INTEGER Nstblsn,Nstnc,Nsbnb,Nglsq

      DOUBLE PRECISION PAR(*),PROB(*),PI,MMIN
      DOUBLE PRECISION LAMBDA,X,Y,Z,O,Q,DZ,E1,E2,SIG,S,SQRS,CONV

      DOUBLE PRECISION GF,MZ,MW,g1,g2,alSMZ,S2TW,ALEMMZ
      DOUBLE PRECISION MGL
      DOUBLE PRECISION MCH2(2),U(2,2,2),V(2,2,2)
      DOUBLE PRECISION MNEU(5),NEU(5,5,2)
      DOUBLE PRECISION MST2(2),UT(2,2,2),MSB2(2),UB(2,2,2),MSL2(2),
     . UTAU(2,2,2),MSNT2
      DOUBLE PRECISION MSU2(2),MSD2(2),MSE2(2),MSNE2,MSMU2(2),
     . UMU(2,2,2)
      DOUBLE PRECISION MST2P(2),MSB2P(2),MSU2P(2),MSD2P(2)
      DOUBLE PRECISION GZMAX,GZINV,GZINVMAX,MCHAMIN,MCMIN,SIGNEU1
      DOUBLE PRECISION SIGNEU,MSLMIN,MSTMIN,MSQMIN,MGLMIN
      DOUBLE PRECISION hZind(1000,2),hZbb(1000,2),hZll(1000,2)
      DOUBLE PRECISION hZinv(1000,2),hZjj(1000,2),hZgg(1000,2)
      DOUBLE PRECISION hA4b(10000,3),hA4tau(10000,3)
      DOUBLE PRECISION hA2b2tau(10000,3),hA2tau2b(10000,3)
      DOUBLE PRECISION AAA6b(10000,3),AAA6tau(10000,3)
      DOUBLE PRECISION AAZ4b(10000,3),AAZ4tau(10000,3)
      DOUBLE PRECISION AAZ2b2tau(10000,3)
      DOUBLE PRECISION cccc02(100,2),cccc04(100,2),cccc05(100,2)
      DOUBLE PRECISION cccc06(100,2),cccc08(100,2),cccc1(100,2)
      DOUBLE PRECISION ccgg02(100,2),ccgg04(100,2),ccgg05(100,2)
      DOUBLE PRECISION ccgg06(100,2),ccgg08(100,2),ccgg1(100,2)
      DOUBLE PRECISION cctt02(100,2),cctt04(100,2),cctt05(100,2)
      DOUBLE PRECISION cctt06(100,2),cctt08(100,2),cctt1(100,2)
      DOUBLE PRECISION gggg02(100,2),gggg04(100,2),gggg05(100,2)
      DOUBLE PRECISION gggg06(100,2),gggg08(100,2),gggg1(100,2)
      DOUBLE PRECISION ttgg02(100,2),ttgg04(100,2),ttgg05(100,2)
      DOUBLE PRECISION ttgg06(100,2),ttgg08(100,2),ttgg1(100,2)
      DOUBLE PRECISION tttt02(100,2),tttt04(100,2),tttt05(100,2)
      DOUBLE PRECISION tttt06(100,2),tttt08(100,2),tttt1(100,2)
      DOUBLE PRECISION stblsn(100,2),stnc(100,2),sbnb(100,2)
      DOUBLE PRECISION glsq(100,2)

      COMMON/EWPAR/GF,MZ,MW,g1,g2,alSMZ,S2TW,ALEMMZ
      COMMON/GLUSPEC/MGL
      COMMON/CHASPEC/MCH2,U,V
      COMMON/NEUSPEC/MNEU,NEU
      COMMON/SFERM3SPEC/MST2,UT,MSB2,UB,MSL2,UTAU,MSNT2
      COMMON/SFERM1SPEC/MSU2,MSD2,MSE2,MSNE2,MSMU2,UMU
      COMMON/SFERMPSPEC/MST2P,MSB2P,MSU2P,MSD2P
      COMMON/LEP/GZMAX,GZINVMAX,MCHAMIN,MCMIN,SIGNEU1,SIGNEU,
     .      MSLMIN,MSTMIN,MSQMIN,MGLMIN,
     .      hZind,hZbb,hZll,hZinv,hZjj,hZgg,
     .      hA4b,hA4tau,hA2b2tau,hA2tau2b,
     .      AAA6b,AAA6tau,AAZ4b,AAZ4tau,AAZ2b2tau,
     .      cccc02,cccc04,cccc05,cccc06,cccc08,cccc1,
     .      ccgg02,ccgg04,ccgg05,ccgg06,ccgg08,ccgg1,
     .      cctt02,cctt04,cctt05,cctt06,cctt08,cctt1,
     .      gggg02,gggg04,gggg05,gggg06,gggg08,gggg1,
     .      ttgg02,ttgg04,ttgg05,ttgg06,ttgg08,ttgg1,
     .      tttt02,tttt04,tttt05,tttt06,tttt08,tttt1,
     .      stblsn,stnc,sbnb,glsq,
     .      NhZind,NhZbb,NhZll,NhZinv,NhZjj,NhZgg,
     .      NhA4b,NhA4tau,NhA2b2tau,NhA2tau2b,
     .      NAAA6b,NAAA6tau,NAAZ4b,NAAZ4tau,NAAZ2b2tau,
     .      Ncccc02,Ncccc04,Ncccc05,Ncccc06,Ncccc08,Ncccc1,
     .      Nccgg02,Nccgg04,Nccgg05,Nccgg06,Nccgg08,Nccgg1,
     .      Ncctt02,Ncctt04,Ncctt05,Ncctt06,Ncctt08,Ncctt1,
     .      Ngggg02,Ngggg04,Ngggg05,Ngggg06,Ngggg08,Ngggg1,
     .      Nttgg02,Nttgg04,Nttgg05,Nttgg06,Nttgg08,Nttgg1,
     .      Ntttt02,Ntttt04,Ntttt05,Ntttt06,Ntttt08,Ntttt1,
     .      Nstblsn,Nstnc,Nsbnb,Nglsq

c********************************************************************************

      LAMBDA(X,Y,Z)= DSQRT(X**2+Y**2+Z**2-2d0*X*Y-2d0*X*Z-2d0*Y*Z)

      PI=4d0*DATAN(1d0)

c      A: LEP searches

* Test on stop -> b l sneutrino

      I=1
      DOWHILE(stblsn(I,1).LE.dsqrt(MST2P(1)).AND. I.LT.Nstblsn)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. dsqrt(MST2P(1)).LT.stblsn(Nstblsn,1))THEN
       MMIN=stblsn(I-1,2)+(dsqrt(MST2P(1))-stblsn(I-1,1))
     .  /(stblsn(I,1)-stblsn(I-1,1))*(stblsn(I,2)-stblsn(I-1,2))
       PROB(20)=DDIM(1d0,dsqrt(MSNE2)/MMIN)
      ENDIF

* Test on stop -> neutralino c

      I=1
      DOWHILE(stnc(I,1).LE.dsqrt(MNEU(1)) .AND. I.LT.Nstnc)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. dsqrt(MNEU(1)).LT.stnc(Nstnc,1))THEN
       MMIN=stnc(I-1,2)+(dsqrt(MNEU(1))-stnc(I-1,1))
     .  /(stnc(I,1)-stnc(I-1,1))*(stnc(I,2)-stnc(I-1,2))
       PROB(21)=DDIM(1d0,dsqrt(MST2P(1))/MMIN)
      ENDIF

* Test on sbottom -> neutralino b

      I=1
      DOWHILE(sbnb(I,1).LE.dsqrt(MSB2P(1)) .AND. I.LT.Nsbnb)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. dsqrt(MSB2P(1)).LT.sbnb(Nsbnb,1))THEN
       MMIN=sbnb(I-1,2)+(dsqrt(MSB2P(1))-sbnb(I-1,1))
     .  /(sbnb(I,1)-sbnb(I-1,1))*(sbnb(I,2)-sbnb(I-1,2))
       PROB(22)=DDIM(1d0,dsqrt(MNEU(1))/MMIN)
      ENDIF

* Test on gluino/squark masses

      I=1
      DOWHILE(glsq(I,1).LE.DABS(MGL) .AND. I.LT.Nglsq)
       I=I+1
      ENDDO
      MMIN=MSQMIN
      IF(I.GT.1 .AND. DABS(MGL).LT.glsq(Nglsq,1))THEN
       MMIN=glsq(I-1,2)+(DABS(MGL)-glsq(I-1,1))
     .  /(glsq(I,1)-glsq(I-1,1))*(glsq(I,2)-glsq(I-1,2))
      ENDIF
      PROB(23)=DDIM(1d0,
     .          dsqrt(MIN(MSU2P(1),MSU2P(2),MSD2P(1),MSD2P(2)))/MMIN)
     .       +DDIM(1d0,DABS(MGL)/MGLMIN)

* Test on slepton masses

      PROB(24)=DDIM(1d0,dsqrt(MIN(MSE2(1),MSE2(2)))/MSLMIN)
      PROB(25)=DDIM(1d0,dsqrt(MSL2(1))/MSTMIN)

* Test on chargino mass

      PROB(1)=DDIM(1d0,dsqrt(MCH2(1))/MCHAMIN)

* Test on Z width into neutralinos

      IF(dsqrt(MNEU(1)).LT.MZ/2d0)THEN
       GZINV=MZ**3*GF/(12d0*DSQRT(2d0)*PI)
     .   *(NEU(1,3,1)**2-NEU(1,4,1)**2+NEU(1,3,2)**2-NEU(1,4,2)**2)**2
     .   *(1d0-4d0*MNEU(1)/MZ**2)**1.5
        PROB(2)=DDIM(GZINV/GZINVMAX,1d0)
      ENDIF

* Test on neutralinos

      SQRS=209d0
      S=SQRS**2
      CONV=.3894D9
      DZ = 1d0/(S-MZ**2)**2

      DO I=1,5
       DO J=I,5
        IF(dsqrt(MNEU(I))+dsqrt(MNEU(J)).LT.SQRS)THEN
         O= (NEU(I,3,1)*NEU(J,3,1)+NEU(I,3,2)*NEU(J,3,2)
     .      -NEU(I,4,1)*NEU(J,4,1)-NEU(I,4,2)*NEU(J,4,2))**2
     .     +(NEU(I,3,1)*NEU(J,3,2)-NEU(I,3,2)*NEU(J,3,1)
     .      -NEU(I,4,1)*NEU(J,4,2)+NEU(I,4,2)*NEU(J,4,1))**2
         E1= (S+MNEU(I)-MNEU(J))/(2d0*SQRS)
         E2= (S+MNEU(J)-MNEU(I))/(2d0*SQRS)
         Q= LAMBDA(S,MNEU(I),MNEU(J))/(2d0*SQRS)
         SIG= 1d0/(4d0*PI)*Q/SQRS*DZ
     .       *(E1*E2+Q**2/3d0-dsqrt(MNEU(I)*MNEU(J)))
     .       *(g1+g2)**2/4d0*O*(.25d0-S2TW+2d0*S2TW**2)
         SIG= CONV*SIG
         IF(I.EQ.J)SIG= SIG/2d0
         IF(I.EQ.1)THEN
          IF(J.GT.1)PROB(2)=PROB(2)+DDIM(SIG/SIGNEU1,1d0)
         ELSE
          PROB(2)=PROB(2)+DDIM(SIG/SIGNEU,1d0)
         ENDIF
        ENDIF
       ENDDO
      ENDDO


c      B: LHC limits -> not implemented

      END


      SUBROUTINE LEP_HIGGS_CPV(PAR,PROB)

***********************************************************************
*   Subroutine to check Higgs LEP constraints:
*      PROB(3) =/= 0   charged Higgs too light
*      PROB(4) =/= 0   excluded by ee -> hZ 
*      PROB(5) =/= 0   excluded by ee -> hZ, h -> bb
*      PROB(6) =/= 0   excluded by ee -> hZ, h -> tautau
*      PROB(7) =/= 0   excluded by ee -> hZ, h -> invisible 
*      PROB(8) =/= 0   excluded by ee -> hZ, h -> 2jets
*      PROB(9) =/= 0   excluded by ee -> hZ, h -> 2photons
*      PROB(10) =/= 0  excluded by ee -> hZ, h -> AA -> 4bs
*      PROB(11) =/= 0  excluded by ee -> hZ, h -> AA -> 4taus
*      PROB(12) =/= 0  excluded by ee -> hZ, h -> AA -> 2bs 2taus
*      PROB(13) =/= 0  excluded by Z -> hA (Z width)
*      PROB(14) =/= 0  excluded by ee -> hA -> 4bs
*      PROB(15) =/= 0  excluded by ee -> hA -> 4taus
*      PROB(16) =/= 0  excluded by ee -> hA -> 2bs 2taus
*      PROB(17) =/= 0  excluded by ee -> hA -> AAA -> 6bs
*      PROB(18) =/= 0  excluded by ee -> hA -> AAA -> 6taus
*      PROB(19) =/= 0  excluded by ee -> Zh -> ZAA -> Z + light pairs
*      PROB(41) =/= 0  excluded by ee -> hZ, h -> AA -> 4taus (ALEPH analysis)
***********************************************************************

      IMPLICIT NONE

      INTEGER I,J,K
      INTEGER NhZind,NhZbb,NhZll,NhZinv,NhZjj,NhZgg
      INTEGER NhA4b,NhA4tau,NhA2b2tau,NhA2tau2b
      INTEGER NAAA6b,NAAA6tau,NAAZ4b,NAAZ4tau,NAAZ2b2tau
      INTEGER Ncccc02,Ncccc04,Ncccc05,Ncccc06,Ncccc08,Ncccc1
      INTEGER Nccgg02,Nccgg04,Nccgg05,Nccgg06,Nccgg08,Nccgg1
      INTEGER Ncctt02,Ncctt04,Ncctt05,Ncctt06,Ncctt08,Ncctt1
      INTEGER Ngggg02,Ngggg04,Ngggg05,Ngggg06,Ngggg08,Ngggg1
      INTEGER Nttgg02,Nttgg04,Nttgg05,Nttgg06,Nttgg08,Nttgg1
      INTEGER Ntttt02,Ntttt04,Ntttt05,Ntttt06,Ntttt08,Ntttt1
      INTEGER Nstblsn,Nstnc,Nsbnb,Nglsq

      DOUBLE PRECISION PAR(*),PROB(*),PI,S,SQRS,SQR2,CONV
      DOUBLE PRECISION LAMBDA,X,Y,Z,MA,MH,BRINV,DMA,Rmax
      DOUBLE PRECISION MAtest,ceff,HMIN,HMAX,M1,M2,R

      DOUBLE PRECISION GF,MZ,MW,g1,g2,alSMZ,S2TW,ALEMMZ
      DOUBLE PRECISION tanb,cosb,sinb,vu,vd,CMASS
      DOUBLE PRECISION MHC2,XC(2,2),MH0(5),XH(5,5),MA2
      DOUBLE PRECISION WIDTH(5),HCWIDTH
      DOUBLE PRECISION BRJJ(5),BREE(5),BRMM(5),BRLL(5),BRSS(5),
     . BRCC(5),BRBB(5),BRTT(5),BRWW(5),BRZZ(5),BRGG(5),BRZG(5)
      DOUBLE PRECISION BRHHH(5,10),BRHCHC(5),BRHAZ(5,4),BRHCW(5),
     . BRHIGGS(5)
      DOUBLE PRECISION BRNEU(5,5,5),BRCHA(5,3),BRHSQ(5,10),BRHSL(5,7),
     . BRSUSY(5)
      DOUBLE PRECISION GZ,GZMAX,GZINV,GZINVMAX,MCHAMIN,MCMIN,SIGNEU1
      DOUBLE PRECISION SIGNEU,MSLMIN,MSTMIN,MSQMIN,MGLMIN
      DOUBLE PRECISION hZind(1000,2),hZbb(1000,2),hZll(1000,2)
      DOUBLE PRECISION hZinv(1000,2),hZjj(1000,2),hZgg(1000,2)
      DOUBLE PRECISION hA4b(10000,3),hA4tau(10000,3)
      DOUBLE PRECISION hA2b2tau(10000,3),hA2tau2b(10000,3)
      DOUBLE PRECISION AAA6b(10000,3),AAA6tau(10000,3)
      DOUBLE PRECISION AAZ4b(10000,3),AAZ4tau(10000,3)
      DOUBLE PRECISION AAZ2b2tau(10000,3)
      DOUBLE PRECISION cccc02(100,2),cccc04(100,2),cccc05(100,2)
      DOUBLE PRECISION cccc06(100,2),cccc08(100,2),cccc1(100,2)
      DOUBLE PRECISION ccgg02(100,2),ccgg04(100,2),ccgg05(100,2)
      DOUBLE PRECISION ccgg06(100,2),ccgg08(100,2),ccgg1(100,2)
      DOUBLE PRECISION cctt02(100,2),cctt04(100,2),cctt05(100,2)
      DOUBLE PRECISION cctt06(100,2),cctt08(100,2),cctt1(100,2)
      DOUBLE PRECISION gggg02(100,2),gggg04(100,2),gggg05(100,2)
      DOUBLE PRECISION gggg06(100,2),gggg08(100,2),gggg1(100,2)
      DOUBLE PRECISION ttgg02(100,2),ttgg04(100,2),ttgg05(100,2)
      DOUBLE PRECISION ttgg06(100,2),ttgg08(100,2),ttgg1(100,2)
      DOUBLE PRECISION tttt02(100,2),tttt04(100,2),tttt05(100,2)
      DOUBLE PRECISION tttt06(100,2),tttt08(100,2),tttt1(100,2)
      DOUBLE PRECISION stblsn(100,2),stnc(100,2),sbnb(100,2)
      DOUBLE PRECISION glsq(100,2)
      DOUBLE PRECISION XIN(6),HM1(6),HM2(6),HM3(6),HM4(6),HM5(6)

      DATA XIN/.2d0,.25d0,.4d0,.5d0,.7d0,1d0/
      DATA HM1/70d0,94d0,99.7d0,102.3d0,105.2d0,108.4d0/
      DATA HM2/70d0,93d0,98.6d0,101.6d0,105d0,108d0/
      DATA HM3/70d0,92.3d0,98.3d0,100.3d0,104.6d0,107.6d0/
      DATA HM4/70d0,90.8d0,98d0,100.1d0,104d0,107d0/
      DATA HM5/70d0,86.7d0,97d0,99d0,103.5d0,106d0/

      COMMON/EWPAR/GF,MZ,MW,g1,g2,alSMZ,S2TW,ALEMMZ
      COMMON/TBPAR/tanb,cosb,sinb,vu,vd
      COMMON/HISPEC/MHC2,XC,MH0,XH,MA2
      COMMON/HIWIDTH/WIDTH,HCWIDTH
      COMMON/HNSMBR/BRJJ,BREE,BRMM,BRLL,BRSS,BRCC,BRBB,BRTT,BRWW,BRZZ,
     . BRGG,BRZG
      COMMON/HNHIBR/BRHHH,BRHCHC,BRHAZ,BRHCW,BRHIGGS
      COMMON/HNSUSYBR/BRNEU,BRCHA,BRHSQ,BRHSL,BRSUSY
      COMMON/LEP/GZMAX,GZINVMAX,MCHAMIN,MCMIN,SIGNEU1,SIGNEU,
     .      MSLMIN,MSTMIN,MSQMIN,MGLMIN,
     .      hZind,hZbb,hZll,hZinv,hZjj,hZgg,
     .      hA4b,hA4tau,hA2b2tau,hA2tau2b,
     .      AAA6b,AAA6tau,AAZ4b,AAZ4tau,AAZ2b2tau,
     .      cccc02,cccc04,cccc05,cccc06,cccc08,cccc1,
     .      ccgg02,ccgg04,ccgg05,ccgg06,ccgg08,ccgg1,
     .      cctt02,cctt04,cctt05,cctt06,cctt08,cctt1,
     .      gggg02,gggg04,gggg05,gggg06,gggg08,gggg1,
     .      ttgg02,ttgg04,ttgg05,ttgg06,ttgg08,ttgg1,
     .      tttt02,tttt04,tttt05,tttt06,tttt08,tttt1,
     .      stblsn,stnc,sbnb,glsq,
     .      NhZind,NhZbb,NhZll,NhZinv,NhZjj,NhZgg,
     .      NhA4b,NhA4tau,NhA2b2tau,NhA2tau2b,
     .      NAAA6b,NAAA6tau,NAAZ4b,NAAZ4tau,NAAZ2b2tau,
     .      Ncccc02,Ncccc04,Ncccc05,Ncccc06,Ncccc08,Ncccc1,
     .      Nccgg02,Nccgg04,Nccgg05,Nccgg06,Nccgg08,Nccgg1,
     .      Ncctt02,Ncctt04,Ncctt05,Ncctt06,Ncctt08,Ncctt1,
     .      Ngggg02,Ngggg04,Ngggg05,Ngggg06,Ngggg08,Ngggg1,
     .      Nttgg02,Nttgg04,Nttgg05,Nttgg06,Nttgg08,Nttgg1,
     .      Ntttt02,Ntttt04,Ntttt05,Ntttt06,Ntttt08,Ntttt1,
     .      Nstblsn,Nstnc,Nsbnb,Nglsq

      LAMBDA(X,Y,Z)= DSQRT(X**2+Y**2+Z**2-2d0*X*Y-2d0*X*Z-2d0*Y*Z)

      PI=4d0*DATAN(1d0)
      SQR2=DSQRT(2d0)
      SQRS=209d0
      S=SQRS**2
      CONV=.3894D9

*   Test on charged Higgs

      CMASS=dsqrt(MHC2)
      PROB(3)=DDIM(1d0,CMASS/MCMIN)

*   Tests on neutral Higgses

c      I- Associated Z / CP-even Higgs production (Higgsstrahlung)
c          ee -> Z* -> hZ

      DO I=1,5

      IF(I.GT.1)THEN
       IF(DABS(dsqrt(MH0(I))-dsqrt(MH0(I-1))).LT.3d0)THEN
        MH=dsqrt(MH0(I-1))
        R=R+(XH(I,1)*vu+XH(I,2)*vd)**2/(vu**2+vd**2)
       ENDIF
      ELSE
       MH=dsqrt(MH0(I))
       R=(XH(I,1)*vu+XH(I,2)*vd)**2/(vu**2+vd**2)
      ENDIF

      IF(MH+MZ.LT.SQRS)THEN

*  ee -> hZ flavor independent

       Rmax=1d0
       J=1
       DOWHILE(hZind(J,1).LE.MH .AND. J.LT.NhZind)
        J=J+1
       ENDDO
       IF(J.GT.1 .AND. MH.LT.hZind(NhZind,1))
     .    Rmax=hZind(J-1,2)+(MH-hZind(J-1,1))/(hZind(J,1)
     .       -hZind(J-1,1))*(hZind(J,2)-hZind(J-1,2))
       PROB(4)=PROB(4)+DDIM(R/Rmax,1d0)

*  ee -> hZ, h -> bb

       Rmax=1d0
       J=1
       DOWHILE(hZbb(J,1).LE.MH .AND. J.LT.NhZbb)
        J=J+1
       ENDDO
       IF(J.GT.1 .AND. MH.LT.hZbb(NhZbb,1))
     .    Rmax=hZbb(J-1,2)+(MH-hZbb(J-1,1))/(hZbb(J,1)-hZbb(J-1,1))
     .       *(hZbb(J,2)-hZbb(J-1,2))
       PROB(5)=PROB(5)+DDIM(R*BRBB(I)/Rmax,1d0)

*  ee -> hZ, h -> tautau

       Rmax=1d0
       J=1
       DOWHILE(hZll(J,1).LE.MH .AND. J.LT.NhZll)
        J=J+1
       ENDDO
       IF(J.GT.1 .AND. MH.LT.hZll(NhZll,1))
     .    Rmax=hZll(J-1,2)+(MH-hZll(J-1,1))/(hZll(J,1)-hZll(J-1,1))
     .       *(hZll(J,2)-hZll(J-1,2))
       PROB(6)=PROB(6)+DDIM(R*BRLL(I)/Rmax,1d0)

*  ee -> hZ, h -> invisible

       Rmax=1d0
       J=1
       DOWHILE(hZinv(J,1).LE.MH .AND. J.LT.NhZinv)
        J=J+1
       ENDDO
       IF(J.GT.1 .AND. MH.LT.hZinv(NhZinv,1))
     .    Rmax=hZinv(J-1,2)+(MH-hZinv(J-1,1))/(hZinv(J,1)-hZinv(J-1,1))
     .       *(hZinv(J,2)-hZinv(J-1,2))
       BRINV=BRNEU(I,1,1)
       IF(I.GT.1)
     .  BRINV=BRINV+BRHHH(I,1)*BRNEU(1,1,1)**2
       IF(I.GT.2)
     .  BRINV=BRINV+BRHHH(I,2)*BRNEU(1,1,1)*BRNEU(2,1,1)
     .             +BRHHH(I,3)*BRNEU(2,1,1)**2
       IF(I.GT.3)
     .  BRINV=BRINV+BRHHH(I,4)*BRNEU(1,1,1)*BRNEU(3,1,1)
     .             +BRHHH(I,5)*BRNEU(2,1,1)*BRNEU(3,1,1)
     .             +BRHHH(I,6)*BRNEU(3,1,1)**2
       IF(I.GT.4)
     .  BRINV=BRINV+BRHHH(I,7)*BRNEU(1,1,1)*BRNEU(4,1,1)
     .             +BRHHH(I,8)*BRNEU(2,1,1)*BRNEU(4,1,1)
     .             +BRHHH(I,9)*BRNEU(3,1,1)*BRNEU(4,1,1)
     .             +BRHHH(I,10)*BRNEU(4,1,1)**2
       PROB(7)=PROB(7)+DDIM(R*BRINV/Rmax,1d0)

*  ee -> hZ, h -> 2jets

       Rmax=1d0
       J=1
       DOWHILE(hZjj(J,1).LE.MH .AND. J.LT.NhZjj)
        J=J+1
       ENDDO
       IF(J.GT.1 .AND. MH.LT.hZjj(NhZjj,1))
     .    Rmax=hZjj(J-1,2)+(MH-hZjj(J-1,1))/(hZjj(J,1)-hZjj(J-1,1))
     .       *(hZjj(J,2)-hZjj(J-1,2))
       PROB(8)=PROB(8)
     .       +DDIM(R*(BRJJ(I)+BRSS(I)+BRCC(I)+BRBB(I))/Rmax,1d0)

*  ee -> hZ, h -> 2photons

       Rmax=1d0
       J=1
       DOWHILE(hZgg(J,1).LE.MH .AND. J.LT.NhZgg)
        J=J+1
       ENDDO
       IF(J.GT.1 .AND. MH.LT.hZgg(NhZgg,1))
     .    Rmax=hZgg(J-1,2)+(MH-hZgg(J-1,1))/(hZgg(J,1)-hZgg(J-1,1))
     .       *(hZgg(J,2)-hZgg(J-1,2))
       PROB(9)=PROB(9)+DDIM(R*BRGG(I)/Rmax,1d0)

*  ee -> hZ, h -> AA -> 4bs

       MA=dsqrt(MH0(1))
       IF(MH.GT.2d0*MA)THEN
        Rmax=1d0
        DO K=1,NAAZ4b
         IF(AINT(MH).EQ.AAZ4b(K,1) .AND. AINT(MA).EQ.AAZ4b(K,2))THEN
          Rmax=AAZ4b(K,3)
          GOTO 1
         ENDIF
        ENDDO
 1      PROB(10)=PROB(10)+DDIM(R*BRHHH(I,1)*BRBB(1)**2/Rmax,1d0)
       ENDIF

       MA=dsqrt(MH0(2))
       IF(MH.GT.2d0*MA)THEN
        Rmax=1d0
        DO K=1,NAAZ4b
         IF(AINT(MH).EQ.AAZ4b(K,1) .AND. AINT(MA).EQ.AAZ4b(K,2))THEN
          Rmax=AAZ4b(K,3)
          GOTO 2
         ENDIF
        ENDDO
 2      PROB(10)=PROB(10)+DDIM(R*BRHHH(I,3)*BRBB(2)**2/Rmax,1d0)
       ENDIF

*  ee -> hZ, h -> AA -> 4taus

       MA=dsqrt(MH0(1))
       IF(MH.GT.2d0*MA)THEN
        Rmax=1d0
        DO K=1,NAAZ4tau
         IF(AINT(MH).EQ.AAZ4tau(K,1) .AND.
     .     AINT(MA).EQ.AAZ4tau(K,2))THEN
          Rmax=AAZ4tau(K,3)
          GOTO 11
         ENDIF
        ENDDO
* Apply only for MA < 9.4GeV where A <-> eta_b mixing is absent:
 11     IF(MA.LE.9.4d0) 
     . PROB(11)=PROB(11)+DDIM(R*BRHHH(I,1)*BRLL(1)**2/Rmax,1d0)
       ENDIF

       MA=dsqrt(MH0(2))
       IF(MH.GT.2d0*MA)THEN
        Rmax=1d0
        DO K=1,NAAZ4tau
         IF(AINT(MH).EQ.AAZ4tau(K,1) .AND.
     .     AINT(MA).EQ.AAZ4tau(K,2))THEN
          Rmax=AAZ4tau(K,3)
          GOTO 12
         ENDIF
        ENDDO
 12     PROB(11)=PROB(11)+DDIM(R*BRHHH(I,3)*BRLL(2)**2/Rmax,1d0)
       ENDIF

*  ee -> hZ, h -> AA -> 4taus (new ALEPH analysis)

       MA=dsqrt(MH0(1))
       RMAX=1d0
       IF(MH.GE.70d0)THEN
        IF(MA.GE.4d0.AND.MA.LT.6d0)THEN
             DMA=(6d0-MA)/2d0
             K=0
 991     K=K+1
             HMIN=HM2(K)+(HM1(K)-HM2(K))*DMA
             HMAX=HM2(K+1)+(HM1(K+1)-HM2(K+1))*DMA
             IF(MH.LE.HMAX)THEN
              RMAX=XIN(K)+(XIN(K+1)-XIN(K))*(MH-HMIN)/(HMAX-HMIN)
             ELSEIF(K.LT.5)THEN
          GOTO 991
         ENDIF
        ELSEIF(MA.GE.6d0.AND.MA.LT.8d0)THEN
         DMA=(8d0-MA)/2d0
              K=0
 992     K=K+1
         HMIN=HM3(K)+(HM2(K)-HM3(K))*DMA
             HMAX=HM3(K+1)+(HM2(K+1)-HM3(K+1))*DMA
             IF(MH.LE.HMAX)THEN
              RMAX=XIN(K)+(XIN(K+1)-XIN(K))*(MH-HMIN)/(HMAX-HMIN)
             ELSEIF(K.LT.5)THEN
          GOTO 992
         ENDIF
        ELSEIF(MA.GE.8d0.AND.MA.LT.10d0)THEN
              DMA=(10d0-MA)/2d0
              K=0
 993     K=K+1
             HMIN=HM4(K)+(HM3(K)-HM4(K))*DMA
             HMAX=HM4(K+1)+(HM3(K+1)-HM4(K+1))*DMA
             IF(MH.LE.HMAX)THEN
              RMAX=XIN(K)+(XIN(K+1)-XIN(K))*(MH-HMIN)/(HMAX-HMIN)
             ELSEIF(K.LT.5)THEN
          GOTO 993
         ENDIF
        ELSEIF(MA.GE.10d0.AND.MA.LE.12d0)THEN
              DMA=(12d0-MA)/2d0
         K=0
 994     K=K+1
             HMIN=HM5(K)+(HM4(K)-HM5(K))*DMA
             HMAX=HM5(K+1)+(HM4(K+1)-HM5(K+1))*DMA
             IF(MH.LE.HMAX)THEN
              RMAX=XIN(K)+(XIN(K+1)-XIN(K))*(MH-HMIN)/(HMAX-HMIN)
             ELSEIF(K.LT.5)THEN
          GOTO 994
         ENDIF
        ENDIF
       ENDIF
* Apply only for MA < 9.4GeV where A <-> eta_b mixing is absent:
       IF(MA.LE.9.4d0) 
     .   PROB(41)=PROB(41)+DDIM(R*BRHHH(I,1)*BRLL(1)**2/Rmax,1d0)


       MA=dsqrt(MH0(2))
       RMAX=1d0
       IF(MH.GE.70d0)THEN
        IF(MA.GE.4d0.AND.MA.LT.6d0)THEN
             DMA=(6d0-MA)/2d0
             K=0
 995     K=K+1
             HMIN=HM2(K)+(HM1(K)-HM2(K))*DMA
             HMAX=HM2(K+1)+(HM1(K+1)-HM2(K+1))*DMA
             IF(MH.LE.HMAX)THEN
              RMAX=XIN(K)+(XIN(K+1)-XIN(K))*(MH-HMIN)/(HMAX-HMIN)
             ELSEIF(K.LT.5)THEN
          GOTO 995
         ENDIF
        ELSEIF(MA.GE.6d0.AND.MA.LT.8d0)THEN
             DMA=(8d0-MA)/2d0
             K=0
 996     K=K+1
             HMIN=HM3(K)+(HM2(K)-HM3(K))*DMA
             HMAX=HM3(K+1)+(HM2(K+1)-HM3(K+1))*DMA
             IF(MH.LE.HMAX)THEN
              RMAX=XIN(K)+(XIN(K+1)-XIN(K))*(MH-HMIN)/(HMAX-HMIN)
             ELSEIF(K.LT.5)THEN
          GOTO 996
         ENDIF
        ELSEIF(MA.GE.8d0.AND.MA.LT.10d0)THEN
             DMA=(10d0-MA)/2d0
             K=0
 997     K=K+1
             HMIN=HM4(K)+(HM3(K)-HM4(K))*DMA
             HMAX=HM4(K+1)+(HM3(K+1)-HM4(K+1))*DMA
             IF(MH.LE.HMAX)THEN
              RMAX=XIN(K)+(XIN(K+1)-XIN(K))*(MH-HMIN)/(HMAX-HMIN)
             ELSEIF(K.LT.5)THEN
          GOTO 997
         ENDIF
        ELSEIF(MA.GE.10d0.AND.MA.LE.12d0)THEN
             DMA=(12d0-MA)/2d0
             K=0
 998     K=K+1
             HMIN=HM5(K)+(HM4(K)-HM5(K))*DMA
             HMAX=HM5(K+1)+(HM4(K+1)-HM5(K+1))*DMA
             IF(MH.LE.HMAX)THEN
              RMAX=XIN(K)+(XIN(K+1)-XIN(K))*(MH-HMIN)/(HMAX-HMIN)
             ELSEIF(K.LT.5)THEN
          GOTO 998
         ENDIF
        ENDIF
       ENDIF
       PROB(41)=PROB(41)+DDIM(R*BRHHH(I,3)*BRLL(2)**2/Rmax,1d0)

*  ee -> hZ, h -> AA -> 2bs 2taus

       MA=dsqrt(MH0(1))
       IF(MH.GT.2d0*MA)THEN
        Rmax=1d0
        DO K=1,NAAZ2b2tau
         IF(AINT(MH).EQ.AAZ2b2tau(K,1) .AND.
     .     AINT(MA).EQ.AAZ2b2tau(K,2))THEN
          Rmax=AAZ2b2tau(K,3)
          GOTO 21
         ENDIF
        ENDDO
* Apply only for MA < 9.4GeV or MA > 10.5GeV without A <-> eta_b mixing:
 21    IF(MA.LE.9.4d0.OR.MA.GE.10.5d0) 
     .   PROB(12)=PROB(12)
     .   +DDIM(R*BRHHH(I,1)*2d0*BRLL(1)*BRBB(1)/Rmax,1d0)
       ENDIF

       MA=dsqrt(MH0(2))
       IF(MH.GT.2d0*MA)THEN
        Rmax=1d0
        DO K=1,NAAZ2b2tau
         IF(AINT(MH).EQ.AAZ2b2tau(K,1) .AND.
     .     AINT(MA).EQ.AAZ2b2tau(K,2))THEN
          Rmax=AAZ2b2tau(K,3)
          GOTO 22
         ENDIF
        ENDDO
 22     PROB(12)=PROB(12)
     .   +DDIM(R*BRHHH(I,3)*2d0*BRBB(2)*BRLL(2)/Rmax,1d0)
       ENDIF


*  ee -> hZ -> AAZ -> Z + light pairs

       MA=dsqrt(MH0(1))
       IF(MH.LT.2d0*MA .OR. MH.LT.40d0 .OR. MH.GT.90d0
     .   .OR. MA.LT.2d0 .OR. MA.GT.12d0)GOTO 752

*      AA -> cccc

       ceff=R*BRHHH(I,1)*BRCC(1)**2
       IF(ceff.GE.0.0.AND.ceff.LT.0.2)GOTO 102
       IF(ceff.GE.0.2.AND.ceff.LT.0.4)GOTO 202
       IF(ceff.GE.0.4.AND.ceff.LT.0.5)GOTO 302
       IF(ceff.GE.0.5.AND.ceff.LT.0.6)GOTO 402
       IF(ceff.GE.0.6.AND.ceff.LT.0.8)GOTO 502
       IF(ceff.GE.0.8.AND.ceff.LT.1.0)GOTO 602
       GOTO 702

 102   Rmax=.2d0
       J=1
       DOWHILE(cccc02(J,1).LE.MH .AND. J.LT.Ncccc02)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LT.cccc02(Ncccc02,1))
     .  MAtest=cccc02(J-1,2)+(MH-cccc02(J-1,1))/
     .   (cccc02(J,1)-cccc02(J-1,1))
     .   *(cccc02(J,2)-cccc02(J-1,2))
       PROB(19)=PROB(19)+5d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
       GOTO 702

 202   Rmax=.4d0
       J=1
       DOWHILE(cccc04(J,1).LE.MH .AND. J.LT.Ncccc04)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LT.cccc04(Ncccc04,1))
     .  MAtest=cccc04(J-1,2)+(MH-cccc04(J-1,1))/
     .   (cccc04(J,1)-cccc04(J-1,1))
     .   *(cccc04(J,2)-cccc04(J-1,2))
       PROB(19)=PROB(19)+5d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
       GOTO 702

 302   Rmax=.5d0
       J=1
       DOWHILE(cccc05(J,1).LE.MH .AND. J.LT.Ncccc05)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LT.cccc05(Ncccc05,1))
     .  MAtest=cccc05(J-1,2)+(MH-cccc05(J-1,1))/
     .   (cccc05(J,1)-cccc05(J-1,1))
     .   *(cccc05(J,2)-cccc05(J-1,2))
       PROB(19)=PROB(19)+10d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
       GOTO 702

 402   Rmax=.6d0
       J=1
       DOWHILE(cccc06(J,1).LE.MH .AND. J.LT.Ncccc06)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LT.cccc06(Ncccc06,1))
     .  MAtest=cccc06(J-1,2)+(MH-cccc06(J-1,1))/
     .   (cccc06(J,1)-cccc06(J-1,1))
     .   *(cccc06(J,2)-cccc06(J-1,2))
       PROB(19)=PROB(19)+10d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
       GOTO 702

 502   Rmax=.8d0
       J=1
       DOWHILE(cccc08(J,1).LE.MH .AND. J.LT.Ncccc08)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LT.cccc08(Ncccc08,1))
     .  MAtest=cccc08(J-1,2)+(MH-cccc08(J-1,1))/
     .   (cccc08(J,1)-cccc08(J-1,1))
     .   *(cccc08(J,2)-cccc08(J-1,2))
       PROB(19)=PROB(19)+5d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
       GOTO 702

 602   Rmax=1d0
       J=1
       DOWHILE(cccc1(J,1).LE.MH .AND. J.LT.Ncccc1)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LT.cccc1(Ncccc1,1))
     .  MAtest=cccc1(J-1,2)+(MH-cccc1(J-1,1))/
     .   (cccc1(J,1)-cccc1(J-1,1))
     .   *(cccc1(J,2)-cccc1(J-1,2))
       PROB(19)=PROB(19)+5d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)

 702   CONTINUE

*      AA -> ccjj

       ceff=R*BRHHH(I,1)*BRCC(1)*(BRJJ(1)+BRSS(1))*2d0
       IF(ceff.GE.0.0.AND.ceff.LT.0.2)GOTO 112
       IF(ceff.GE.0.2.AND.ceff.LT.0.4)GOTO 212
       IF(ceff.GE.0.4.AND.ceff.LT.0.5)GOTO 312
       IF(ceff.GE.0.5.AND.ceff.LT.0.6)GOTO 412
       IF(ceff.GE.0.6.AND.ceff.LT.0.8)GOTO 512
       IF(ceff.GE.0.8.AND.ceff.LT.1.0)GOTO 612
       GOTO 712

 112   Rmax=.2d0
       J=1
       DOWHILE(ccgg02(J,1).LE.MH .AND. J.LT.Nccgg02)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LT.ccgg02(Nccgg02,1))
     .  MAtest=ccgg02(J-1,2)+(MH-ccgg02(J-1,1))/
     .   (ccgg02(J,1)-ccgg02(J-1,1))
     .   *(ccgg02(J,2)-ccgg02(J-1,2))
       PROB(19)=PROB(19)+5d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
       GOTO 712

 212   Rmax=.4d0
       J=1
       DOWHILE(ccgg04(J,1).LE.MH .AND. J.LT.Nccgg04)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LT.ccgg04(Nccgg04,1))
     .  MAtest=ccgg04(J-1,2)+(MH-ccgg04(J-1,1))/
     .   (ccgg04(J,1)-ccgg04(J-1,1))
     .   *(ccgg04(J,2)-ccgg04(J-1,2))
       PROB(19)=PROB(19)+5d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
       GOTO 712

 312   Rmax=.5d0
       J=1
       DOWHILE(ccgg05(J,1).LE.MH .AND. J.LT.Nccgg05)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LT.ccgg05(Nccgg05,1))
     .  MAtest=ccgg05(J-1,2)+(MH-ccgg05(J-1,1))/
     .   (ccgg05(J,1)-ccgg05(J-1,1))
     .   *(ccgg05(J,2)-ccgg05(J-1,2))
       PROB(19)=PROB(19)+10d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
       GOTO 712

 412   Rmax=.6d0
       J=1
       DOWHILE(ccgg06(J,1).LE.MH .AND. J.LT.Nccgg06)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LT.ccgg06(Nccgg06,1))
     .    MAtest=ccgg06(J-1,2)+(MH-ccgg06(J-1,1))/
     .   (ccgg06(J,1)-ccgg06(J-1,1))
     .   *(ccgg06(J,2)-ccgg06(J-1,2))
       PROB(19)=PROB(19)+10d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
       GOTO 712

 512   Rmax=.8d0
       J=1
       DOWHILE(ccgg08(J,1).LE.MH .AND. J.LT.Nccgg08)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LT.ccgg08(Nccgg08,1))
     .  MAtest=ccgg08(J-1,2)+(MH-ccgg08(J-1,1))/
     .   (ccgg08(J,1)-ccgg08(J-1,1))
     .   *(ccgg08(J,2)-ccgg08(J-1,2))
       PROB(19)=PROB(19)+5d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
       GOTO 712

 612   Rmax=1d0
       J=1
       DOWHILE(ccgg1(J,1).LE.MH .AND. J.LT.Nccgg1)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LT.ccgg1(Nccgg1,1))
     .  MAtest=ccgg1(J-1,2)+(MH-ccgg1(J-1,1))/
     .   (ccgg1(J,1)-ccgg1(J-1,1))
     .   *(ccgg1(J,2)-ccgg1(J-1,2))
       PROB(19)=PROB(19)+5d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
      
 712   CONTINUE

*      AA -> cctautau

       ceff=R*BRHHH(I,1)*BRCC(1)*BRLL(1)*2d0
       IF(ceff.GE.0.0.AND.ceff.LT.0.2)GOTO 122
       IF(ceff.GE.0.2.AND.ceff.LT.0.4)GOTO 222
       IF(ceff.GE.0.4.AND.ceff.LT.0.5)GOTO 322
       IF(ceff.GE.0.5.AND.ceff.LT.0.6)GOTO 422
       IF(ceff.GE.0.6.AND.ceff.LT.0.8)GOTO 522
       IF(ceff.GE.0.8.AND.ceff.LT.1.0)GOTO 622
       GOTO 722

 122   Rmax=.2d0
       J=1
       DOWHILE(cctt02(J,1).LE.MH .AND. J.LT.Ncctt02)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LT.cctt02(Ncctt02,1))
     .  MAtest=cctt02(J-1,2)+(MH-cctt02(J-1,1))/
     .   (cctt02(J,1)-cctt02(J-1,1))
     .   *(cctt02(J,2)-cctt02(J-1,2))
       PROB(19)=PROB(19)+5d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
       GOTO 722

 222   Rmax=.4d0
       J=1
       DOWHILE(cctt04(J,1).LE.MH .AND. J.LT.Ncctt04)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LT.cctt04(Ncctt04,1))
     .  MAtest=cctt04(J-1,2)+(MH-cctt04(J-1,1))/
     .   (cctt04(J,1)-cctt04(J-1,1))
     .   *(cctt04(J,2)-cctt04(J-1,2))
       PROB(19)=PROB(19)+5d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
       GOTO 722

 322   Rmax=.5d0
       J=1
       DOWHILE(cctt05(J,1).LE.MH .AND. J.LT.Ncctt05)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LT.cctt05(Ncctt05,1))
     .  MAtest=cctt05(J-1,2)+(MH-cctt05(J-1,1))/
     .   (cctt05(J,1)-cctt05(J-1,1))
     .   *(cctt05(J,2)-cctt05(J-1,2))
       PROB(19)=PROB(19)+10d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
       GOTO 722

 422   Rmax=.6d0
       J=1
       DOWHILE(cctt06(J,1).LE.MH .AND. J.LT.Ncctt06)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LT.cctt06(Ncctt06,1))
     .  MAtest=cctt06(J-1,2)+(MH-cctt06(J-1,1))/
     .   (cctt06(J,1)-cctt06(J-1,1))
     .   *(cctt06(J,2)-cctt06(J-1,2))
       PROB(19)=PROB(19)+10d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
       GOTO 722

 522   Rmax=.8d0
       J=1
       DOWHILE(cctt08(J,1).LE.MH .AND. J.LT.Ncctt08)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LT.cctt08(Ncctt08,1))
     .  MAtest=cctt08(J-1,2)+(MH-cctt08(J-1,1))/
     .   (cctt08(J,1)-cctt08(J-1,1))
     .   *(cctt08(J,2)-cctt08(J-1,2))
       PROB(19)=PROB(19)+5d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
       GOTO 722

 622   Rmax=1d0
       J=1
       DOWHILE(cctt1(J,1).LE.MH .AND. J.LT.Ncctt1)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LT.cctt1(Ncctt1,1))
     .  MAtest=cctt1(J-1,2)+(MH-cctt1(J-1,1))/
     .   (cctt1(J,1)-cctt1(J-1,1))
     .   *(cctt1(J,2)-cctt1(J-1,2))
       PROB(19)=PROB(19)+5d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
      
 722   CONTINUE

*      AA -> jjjj

       ceff=R*BRHHH(I,1)*(BRJJ(1)+BRSS(1))**2
       IF(ceff.GE.0.0.AND.ceff.LT.0.2)GOTO 132
       IF(ceff.GE.0.2.AND.ceff.LT.0.4)GOTO 232
       IF(ceff.GE.0.4.AND.ceff.LT.0.5)GOTO 332
       IF(ceff.GE.0.5.AND.ceff.LT.0.6)GOTO 432
       IF(ceff.GE.0.6.AND.ceff.LT.0.8)GOTO 532
       IF(ceff.GE.0.8.AND.ceff.LT.1.0)GOTO 632
       GOTO 732

 132   Rmax=.2d0
       J=1
       DOWHILE(gggg02(J,1).LE.MH .AND. J.LT.Ngggg02)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LT.gggg02(Ngggg02,1))
     .  MAtest=gggg02(J-1,2)+(MH-gggg02(J-1,1))/
     .   (gggg02(J,1)-gggg02(J-1,1))
     .   *(gggg02(J,2)-gggg02(J-1,2))
       PROB(19)=PROB(19)+5d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
       GOTO 732

 232   Rmax=.4d0
       J=1
       DOWHILE(gggg04(J,1).LE.MH .AND. J.LT.Ngggg04)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LT.gggg04(Ngggg04,1))
     .  MAtest=gggg04(J-1,2)+(MH-gggg04(J-1,1))/
     .   (gggg04(J,1)-gggg04(J-1,1))
     .   *(gggg04(J,2)-gggg04(J-1,2))
       PROB(19)=PROB(19)+5d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
       GOTO 732

 332   Rmax=.5d0
       J=1
       DOWHILE(gggg05(J,1).LE.MH .AND. J.LT.Ngggg05)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LT.gggg05(Ngggg05,1))
     .  MAtest=gggg05(J-1,2)+(MH-gggg05(J-1,1))/
     .   (gggg05(J,1)-gggg05(J-1,1))
     .   *(gggg05(J,2)-gggg05(J-1,2))
       PROB(19)=PROB(19)+10d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
       GOTO 732

 432   Rmax=.6d0
       J=1
       DOWHILE(gggg06(J,1).LE.MH .AND. J.LT.Ngggg06)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LT.gggg06(Ngggg06,1))
     .  MAtest=gggg06(J-1,2)+(MH-gggg06(J-1,1))/
     .   (gggg06(J,1)-gggg06(J-1,1))
     .   *(gggg06(J,2)-gggg06(J-1,2))
       PROB(19)=PROB(19)+10d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
       GOTO 732

 532   Rmax=.8d0
       J=1
       DOWHILE(gggg08(J,1).LE.MH .AND. J.LT.Ngggg08)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LT.gggg08(Ngggg08,1))
     .  MAtest=gggg08(J-1,2)+(MH-gggg08(J-1,1))/
     .   (gggg08(J,1)-gggg08(J-1,1))
     .   *(gggg08(J,2)-gggg08(J-1,2))
       PROB(19)=PROB(19)+5d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
       GOTO 732

 632   Rmax=1d0
       J=1
       DOWHILE(gggg1(J,1).LE.MH .AND. J.LT.Ngggg1)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LT.gggg1(Ngggg1,1))
     .  MAtest=gggg1(J-1,2)+(MH-gggg1(J-1,1))/
     .   (gggg1(J,1)-gggg1(J-1,1))
     .   *(gggg1(J,2)-gggg1(J-1,2))
       PROB(19)=PROB(19)+5d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
      
 732   CONTINUE

*      AA -> tautaujj

       ceff=R*BRHHH(I,1)*BRLL(1)*(BRJJ(1)+BRSS(1))*2d0
       IF(ceff.GE.0.0.AND.ceff.LT.0.2)GOTO 142
       IF(ceff.GE.0.2.AND.ceff.LT.0.4)GOTO 242
       IF(ceff.GE.0.4.AND.ceff.LT.0.5)GOTO 342
       IF(ceff.GE.0.5.AND.ceff.LT.0.6)GOTO 442
       IF(ceff.GE.0.6.AND.ceff.LT.0.8)GOTO 542
       IF(ceff.GE.0.8.AND.ceff.LT.1.0)GOTO 642
       GOTO 742

 142   Rmax=.2d0
       J=1
       DOWHILE(ttgg02(J,1).LE.MH .AND. J.LT.Nttgg02)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LT.ttgg02(Nttgg02,1))
     .  MAtest=ttgg02(J-1,2)+(MH-ttgg02(J-1,1))/
     .   (ttgg02(J,1)-ttgg02(J-1,1))
     .   *(ttgg02(J,2)-ttgg02(J-1,2))
       PROB(19)=PROB(19)+5d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
       GOTO 742

 242   Rmax=.4d0
       J=1
       DOWHILE(ttgg04(J,1).LE.MH .AND. J.LT.Nttgg04)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LT.ttgg04(Nttgg04,1))
     .  MAtest=ttgg04(J-1,2)+(MH-ttgg04(J-1,1))/
     .   (ttgg04(J,1)-ttgg04(J-1,1))
     .   *(ttgg04(J,2)-ttgg04(J-1,2))
       PROB(19)=PROB(19)+5d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
       GOTO 742

 342   Rmax=.5d0
       J=1
       DOWHILE(ttgg05(J,1).LE.MH .AND. J.LT.Nttgg05)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LT.ttgg05(Nttgg05,1))
     .  MAtest=ttgg05(J-1,2)+(MH-ttgg05(J-1,1))/
     .   (ttgg05(J,1)-ttgg05(J-1,1))
     .   *(ttgg05(J,2)-ttgg05(J-1,2))
       PROB(19)=PROB(19)+10d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
       GOTO 742

 442   Rmax=.6d0
       J=1
       DOWHILE(ttgg06(J,1).LE.MH .AND. J.LT.Nttgg06)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LT.ttgg06(Nttgg06,1))
     .  MAtest=ttgg06(J-1,2)+(MH-ttgg06(J-1,1))/
     .   (ttgg06(J,1)-ttgg06(J-1,1))
     .   *(ttgg06(J,2)-ttgg06(J-1,2))
       PROB(19)=PROB(19)+10d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
       GOTO 742

 542   Rmax=.8d0
       J=1
       DOWHILE(ttgg08(J,1).LE.MH .AND. J.LT.Nttgg08)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LT.ttgg08(Nttgg08,1))
     .  MAtest=ttgg08(J-1,2)+(MH-ttgg08(J-1,1))/
     .   (ttgg08(J,1)-ttgg08(J-1,1))
     .   *(ttgg08(J,2)-ttgg08(J-1,2))
       PROB(19)=PROB(19)+5d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
       GOTO 742

 642   Rmax=1d0
       J=1
       DOWHILE(ttgg1(J,1).LE.MH .AND. J.LT.Nttgg1)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LT.ttgg1(Nttgg1,1))
     .  MAtest=ttgg1(J-1,2)+(MH-ttgg1(J-1,1))/
     .   (ttgg1(J,1)-ttgg1(J-1,1))
     .   *(ttgg1(J,2)-ttgg1(J-1,2))
       PROB(19)=PROB(19)+5d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
      
 742   CONTINUE

*      AA -> tautautautau

       ceff=R*BRHHH(I,1)*BRLL(1)**2
       IF(ceff.GE.0.0.AND.ceff.LT.0.2)GOTO 152
       IF(ceff.GE.0.2.AND.ceff.LT.0.4)GOTO 252
       IF(ceff.GE.0.4.AND.ceff.LT.0.5)GOTO 352
       IF(ceff.GE.0.5.AND.ceff.LT.0.6)GOTO 452
       IF(ceff.GE.0.6.AND.ceff.LT.0.8)GOTO 552
       IF(ceff.GE.0.8.AND.ceff.LT.1.0)GOTO 652
       GOTO 752

 152   Rmax=.2d0
       J=1
       DOWHILE(tttt02(J,1).LE.MH .AND. J.LT.Ntttt02)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LT.tttt02(Ntttt02,1))
     .  MAtest=tttt02(J-1,2)+(MH-tttt02(J-1,1))/
     .   (tttt02(J,1)-tttt02(J-1,1))
     .   *(tttt02(J,2)-tttt02(J-1,2))
       PROB(19)=PROB(19)+5d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
       GOTO 752

 252   Rmax=.4d0
       J=1
       DOWHILE(tttt04(J,1).LE.MH .AND. J.LT.Ntttt04)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LT.tttt04(Ntttt04,1))
     .  MAtest=tttt04(J-1,2)+(MH-tttt04(J-1,1))/
     .   (tttt04(J,1)-tttt04(J-1,1))
     .   *(tttt04(J,2)-tttt04(J-1,2))
       PROB(19)=PROB(19)+5d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
       GOTO 752

 352   Rmax=.5d0
       J=1
       DOWHILE(tttt05(J,1).LE.MH .AND. J.LT.Ntttt05)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LT.tttt05(Ntttt05,1))
     .  MAtest=tttt05(J-1,2)+(MH-tttt05(J-1,1))/
     .   (tttt05(J,1)-tttt05(J-1,1))
     .   *(tttt05(J,2)-tttt05(J-1,2))
       PROB(19)=PROB(19)+10d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
       GOTO 752

 452   Rmax=.6d0
       J=1
       DOWHILE(tttt06(J,1).LE.MH .AND. J.LT.Ntttt06)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LT.tttt06(Ntttt06,1))
     .  MAtest=tttt06(J-1,2)+(MH-tttt06(J-1,1))/
     .   (tttt06(J,1)-tttt06(J-1,1))
     .   *(tttt06(J,2)-tttt06(J-1,2))
       PROB(19)=PROB(19)+10d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
       GOTO 752

 552   Rmax=.8d0
       J=1
       DOWHILE(tttt08(J,1).LE.MH .AND. J.LT.Ntttt08)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LT.tttt08(Ntttt08,1))
     .  MAtest=tttt08(J-1,2)+(MH-tttt08(J-1,1))/
     .   (tttt08(J,1)-tttt08(J-1,1))
     .   *(tttt08(J,2)-tttt08(J-1,2))
       PROB(19)=PROB(19)+5d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
       GOTO 752

 652   Rmax=1d0
       J=1
       DOWHILE(tttt1(J,1).LE.MH .AND. J.LT.Ntttt1)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LT.tttt1(Ntttt1,1))
     .  MAtest=tttt1(J-1,2)+(MH-tttt1(J-1,1))/
     .   (tttt1(J,1)-tttt1(J-1,1))
     .   *(tttt1(J,2)-tttt1(J-1,2))
       PROB(19)=PROB(19)+5d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
      
 752   CONTINUE
      ENDIF

      ENDDO


c      II- Pair production: ee -> Z* -> h_i+h_j

      DO I=1,5
      DO J=1,I

      MH=dsqrt(MH0(I))
      MA=dsqrt(MH0(J))
      R=(XH(I,4)*(XH(J,1)*vd-XH(J,2)*vu)
     .  -XH(J,4)*(XH(I,1)*vd-XH(I,2)*vu))**2/(vu**2+vd**2)
      IF(I.eq.J)R=R/2d0

*  Z width

      IF(MH+MA.LT.MZ)THEN
       GZ=SQR2*GF*R*LAMBDA(MZ**2,MA**2,MH**2)**3/(48d0*PI*MZ**3)
        PROB(13)=PROB(13)+DDIM(GZ/GZMAX,1d0)
      ENDIF

*  ee -> hA -> 4bs

      IF(MH+MA.LT.SQRS)THEN

       M1=MIN(MH,MA)
       M2=MAX(MH,MA)

       Rmax=1d0
       DO K=1,NhA4b
        IF(AINT(M2).EQ.hA4b(K,1) .AND. AINT(M1).EQ.hA4b(K,2))THEN
         Rmax=hA4b(K,3)
         GOTO 3
        ENDIF
       ENDDO
 3     PROB(14)=PROB(14)+DDIM(R*BRBB(I)*BRBB(J)/Rmax,1d0)

*  ee -> hA -> 4taus

       Rmax=1d0
       DO K=1,NhA4tau
        IF(AINT(M2).EQ.hA4tau(K,1) .AND. AINT(M1).EQ.hA4tau(K,2))THEN
         Rmax=hA4tau(K,3)
         GOTO 4
        ENDIF
       ENDDO
 4     PROB(15)=PROB(15)+DDIM(R*BRLL(I)*BRLL(J)/Rmax,1d0)

*  ee -> hA -> 2b 2taus

       Rmax=1d0
       IF(MA.GT.MH)THEN
        DO K=1,NhA2b2tau
         IF(AINT(M2).EQ.hA2b2tau(K,1) .AND.
     .     AINT(M1).EQ.hA2b2tau(K,2))THEN
          Rmax=hA2b2tau(K,3)
          GOTO 5
         ENDIF
        ENDDO
 5      PROB(16)=PROB(16)+DDIM(R*BRLL(I)*BRBB(J)/Rmax,1d0)
       ENDIF

       Rmax=1d0
       IF(MH.GT.MA)THEN
        DO K=1,NhA2b2tau
         IF(AINT(M2).EQ.hA2b2tau(K,1) .AND.
     .     AINT(M1).EQ.hA2b2tau(K,2))THEN
          Rmax=hA2b2tau(K,3)
          GOTO 6
         ENDIF
        ENDDO
 6      PROB(16)=PROB(16)+DDIM(R*BRBB(I)*BRLL(J)/Rmax,1d0)
       ENDIF

       Rmax=1d0
       IF(MA.GT.MH)THEN
        DO K=1,NhA2tau2b
         IF(AINT(M2).EQ.hA2tau2b(K,1) .AND.
     .     AINT(M1).EQ.hA2tau2b(K,2))THEN
          Rmax=hA2tau2b(K,3)
          GOTO 7
         ENDIF
        ENDDO
 7      PROB(16)=PROB(16)+DDIM(R*BRBB(I)*BRLL(J)/Rmax,1d0)
       ENDIF

       Rmax=1d0
       IF(MH.GT.MA)THEN
        DO K=1,NhA2tau2b
         IF(AINT(M2).EQ.hA2tau2b(K,1) .AND.
     .     AINT(M1).EQ.hA2tau2b(K,2))THEN
          Rmax=hA2tau2b(K,3)
          GOTO 8
         ENDIF
        ENDDO
 8      PROB(16)=PROB(16)+DDIM(R*BRLL(I)*BRBB(J)/Rmax,1d0)
       ENDIF

*  ee -> hA -> AAA -> 6bs

       IF(MH.GT.2d0*MA)THEN
        Rmax=1d0
        DO K=1,NAAA6b
         IF(AINT(MH).EQ.AAA6b(K,1) .AND. AINT(MA).EQ.AAA6b(K,2))THEN
          Rmax=AAA6b(K,3)
          GOTO 9
         ENDIF
        ENDDO
 9      PROB(17)=PROB(17)+DDIM(R*BRHHH(I,J*(J+1)/2)*BRBB(J)**3/Rmax,
     .   1d0)

*  ee -> hA -> AAA -> 6taus

        Rmax=1d0
        DO K=1,NAAA6tau
         IF(AINT(MH).EQ.AAA6tau(K,1) .AND.
     .     AINT(MA).EQ.AAA6tau(K,2))THEN
          Rmax=AAA6tau(K,3)
          GOTO 10
         ENDIF
        ENDDO
 10     PROB(18)=PROB(18)+DDIM(R*BRHHH(I,J*(J+1)/2)*BRLL(J)**3/Rmax,
     .   1d0)
      ENDIF

      ENDIF

      ENDDO
      ENDDO

      END


      SUBROUTINE TEVATRON_CHIGGS_CPV(PAR,PROB)

***********************************************************************
*   Subroutine to check TeVatron constraints on a charged Higgs:
*      PROB(42) =/= 0  excluded by top -> b H+, H+ -> c s (CDF, D0)
*      PROB(43) =/= 0  excluded by top -> b H+, H+ -> tau nu_tau (D0)
*      PROB(44) =/= 0  excluded by top -> b H+, H+ -> W+ A1, A1 -> 2taus (CDF)
***********************************************************************

      IMPLICIT NONE

      INTEGER I,J

      DOUBLE PRECISION PAR(*),PROB(*),MA,CMASS,brtoplim
      DOUBLE PRECISION minf(9),msup(9),binf(9),bsup(9)
      DOUBLE PRECISION mh4m(4),mh4p(4),br4m(4),br4p(4),brtau(8)
      DOUBLE PRECISION mh7(5),br7(5),br9(5)
      DOUBLE PRECISION brtopcs,brtoptau,brtopa1

      DOUBLE PRECISION MHC2,XC(2,2),MH0(5),XH(5,5),MA2
      DOUBLE PRECISION brtopbw,brtopbh,brtopneutrstop(5,2)
      DOUBLE PRECISION HCBRM,HCBRL,HCBRSU,HCBRBU,HCBRSC,HCBRBC,HCBRBT
      DOUBLE PRECISION HCBRWH(5),HCBRWHT
      DOUBLE PRECISION BRJJ(5),BREE(5),BRMM(5),BRLL(5),BRSS(5),
     . BRCC(5),BRBB(5),BRTT(5),BRWW(5),BRZZ(5),BRGG(5),BRZG(5)

      DATA minf/60d0,70d0,80d0,100d0,110d0,120d0,130d0,140d0,150d0/
      DATA msup/70d0,80d0,100d0,110d0,120d0,130d0,140d0,150d0,155d0/
      DATA binf/.09d0,0d0,.21d0,.21d0,.15d0,.12d0,.08d0,.1d0,.20d0/
      DATA bsup/.12d0,0d0,.21d0,.15d0,.12d0,.08d0,.1d0,.13d0,.19d0/
      DATA mh4m/85d0,90d0,100d0,120d0/
      DATA mh4p/90d0,100d0,120d0,140d0/
      DATA br4m/.5d0,.33d0,.27d0,.35d0/
      DATA br4p/.33d0,.27d0,.35d0,.52d0/
      DATA mh7/90d0,100d0,120d0,140d0,160d0/
      DATA br7/.13d0,.1d0,.13d0,.2d0,.3d0/
      DATA br9/.11d0,.08d0,.09d0,.14d0,.21d0/
      DATA brtau/.16d0,.15d0,.16,.17d0,.175,.18d0,.19d0,.18d0/

      COMMON/HISPEC/MHC2,XC,MH0,XH,MA2
      COMMON/BR_top2body/brtopbw,brtopbh,brtopneutrstop
      COMMON/HCSMBR/HCBRM,HCBRL,HCBRSU,HCBRBU,HCBRSC,HCBRBC,HCBRBT
      COMMON/HCHIBR/HCBRWH,HCBRWHT
      COMMON/HNSMBR/BRJJ,BREE,BRMM,BRLL,BRSS,BRCC,BRBB,BRTT,BRWW,BRZZ,
     . BRGG,BRZG

***********************************************************************

      CMASS=dsqrt(MHC2)

* top -> H+ b, H+ -> c s (from CDF, 0907.1269, and D0, 0908.1811)
       
       brtopcs=brtopbh*HCBRSC
       PROB(42)=0d0
       
       DO I=1,9
        IF((minf(i).LE.CMASS).AND.(CMASS.LE.msup(i))) THEN
        brtoplim=binf(i)+
     .      (CMASS-minf(i))*(bsup(i)-binf(i))/(msup(i)-minf(i))
         PROB(42)=PROB(42)+DDIM(brtopcs/brtoplim,1d0)
        ENDIF
       ENDDO

* top -> H+ b, H+ -> tau nu_tau (from D0, 0908.1811)

       brtoptau=brtopbh*HCBRL
       PROB(43)=0d0
    
       DO I=1,7
        IF((minf(i+2).LE.CMASS).AND.(CMASS.LE.msup(i+2))) THEN
        brtoplim=brtau(i)+
     .    (CMASS-minf(i+2))*(brtau(i+1)-brtau(i))/(msup(i+2)-minf(i+2))
         PROB(43)=PROB(43)+DDIM(brtoptau/brtoplim,1d0)
        ENDIF
       ENDDO

* top -> H+ b, H+ -> W+ A_1, A_1 -> 2taus (from CDF Note 10104)

       PROB(44)=0d0

      DO J=1,4
       MA=dsqrt(MH0(J))
       brtopa1=brtopbh*HCBRWH(J)*BRLL(J)

      IF(MA.LE.4d0) then     
       DO I=1,4
        IF((mh4m(i).LE.CMASS).AND.(CMASS.LE.mh4p(i))) THEN
          brtoplim=br4m(i)+
     .      (CMASS-mh4m(i))*(br4p(i)-br4m(i))/(mh4p(i)-mh4m(i))
         PROB(44)=PROB(44)+DDIM(brtopa1/brtoplim,1d0)
        ENDIF
       ENDDO
      ENDIF

      IF((MA.gt.4d0).and.(MA.LE.7d0)) then     
       DO I=1,3
        IF((mh7(i).LE.CMASS).AND.(CMASS.LE.mh7(i+1))) THEN
          brtoplim=br4p(i)+(MA-4d0)*(br7(i)-br4p(i))/3d0
     .     +(br4p(i+1)-br4p(i)
     .      +(br7(i+1)-br7(i)-br4p(i+1)+br4p(i))*(MA-4d0)/3d0)
     .     *(CMASS-mh7(i))/(mh7(i+1)-mh7(i))
         PROB(44)=PROB(44)+DDIM(brtopa1/brtoplim,1d0)
        ENDIF
       ENDDO
      ENDIF

      IF((MA.gt.7d0).and.(MA.LE.9d0)) then     
       DO I=1,4
        IF((mh7(i).LE.CMASS).AND.(CMASS.LE.mh7(i+1))) THEN
          brtoplim=br7(i)+(MA-7d0)*(br9(i)-br7(i))/2d0
     .     +(br7(i+1)-br7(i)
     .      +(br9(i+1)-br9(i)-br7(i+1)+br7(i))*(MA-7d0)/2d0)
     .     *(CMASS-mh7(i))/(mh7(i+1)-mh7(i))
         PROB(44)=PROB(44)+DDIM(brtopa1/brtoplim,1d0)
        ENDIF
       ENDDO
      ENDIF

      ENDDO

      END


      SUBROUTINE LHC_HIGGS_CPV(PAR,PROB)

***********************************************************************
*   Subroutine to check Higgs LHC constraints:
*      PROB(45) =/= 0: Bound on Br(t->bH+)*BR(H+->tau nu)
*      PROB(46) =/= 0:  No Higgs in the MHmin-MHmax GeV range
*      PROB(47) =/= 0:  chi2gam > chi2max
*      PROB(48) =/= 0:  chi2bb > chi2max
*      PROB(49) =/= 0:  chi2zz > chi2max
*      PROB(51) =/= 0: excluded by ggF/bb->H/A->tautau (CMS-PAS-HIG-13-021)
*      PROB(52) =/= 0: Excluded H_125->AA->4mu (CMS)
***********************************************************************

      IMPLICIT NONE

      CHARACTER*256 FILENAME,EXPCON_PATH,catpath

      INTEGER I,I1,J,J1,NX,NY,JBAR,JBARbb,K,K1,K2,K3
      PARAMETER(NX=18)
      PARAMETER(NY=6)

      DOUBLE PRECISION PAR(*),PROB(*)
      DOUBLE PRECISION HMAS(NX),XSM(NX),XSMbb(NX),LCMS(NX),LCMSbb(NX)
      DOUBLE PRECISION LATLASgg(NX),LATLASbb(NX)
      DOUBLE PRECISION MH(5),XSMH(5),LCMSH(5),SIG1(5),LATLASH(5)
      DOUBLE PRECISION XSMHbb(5),LCMSHbb(5),SIGbb(5),LATLASHbb(5)
      DOUBLE PRECISION DEL,SIGTOT,MBAR,LCMSMB,LATLASMB
      DOUBLE PRECISION SIGTOTbb,MBARbb,LCMSMBbb,LATLASMBbb
      DOUBLE PRECISION X1(NY),LSIGBR(NY),LIMIT,LHC_TBH
      DOUBLE PRECISION D1,D2,CJ2,BRHTOAA,MHcen,masstest
      DOUBLE PRECISION agg,bgg,cgg,mugcengg,muvcengg,SSIG(6)
      DOUBLE PRECISION abb,bbb,cbb,mugcenbb,muvcenbb
      DOUBLE PRECISION azz,bzz,czz,mugcenzz,muvcenzz
      DOUBLE PRECISION BRJJSM,BREESM,BRMMSM,BRLLSM,BRSSSM,BRCCSM
      DOUBLE PRECISION BRBBSM,BRTTSM,BRWWSM,BRZZSM,BRGGSM,BRZGSM
      DOUBLE PRECISION ggHgg(100,2),dummy(100,4),SMXS(100,2)

      DOUBLE PRECISION MHC,XC(2,2),MH0(5),XH(5,5),MA2
      DOUBLE PRECISION CU(5),CUP(5),CD(5),CDP(5),CB(5),CBP(5),CJ(5),
     . CJP(5),CI(5),CG(5),CGP(5),CV(5),CZG(5),CZGP(5)
      DOUBLE PRECISION BRJJ(5),BREE(5),BRMM(5),BRLL(5),BRSS(5),
     . BRCC(5),BRBB(5),BRTT(5),BRWW(5),BRZZ(5),BRGG(5),BRZG(5)
      DOUBLE PRECISION BRHHH(5,10),BRHCHC(5),BRHAZ(5,4),BRHCW(5),
     . BRHIGGS(5)
      DOUBLE PRECISION BRNEU(5,5,5),BRCHA(5,3),BRHSQ(5,10),BRHSL(5,7),
     . BRSUSY(5)
      DOUBLE PRECISION brtopbw,brtopbh,brtopneutrstop(5,2)
      DOUBLE PRECISION HCBRM,HCBRL,HCBRSU,HCBRBU,HCBRSC,HCBRBC,HCBRBT
      DOUBLE PRECISION MHmin,MHmax,chi2max,chi2gam,chi2bb,chi2zz
      DOUBLE PRECISION SIG(5,10),muH2

      DATA HMAS/90d0,100d0,120d0,140d0,160d0,180d0,200d0,250d0,300d0,
     . 350d0,400d0,450d0,500d0,600d0,700d0,800d0,900d0,1000d0/

* SM Higgs ggF prod. cross sect. at 8 TeV from 
* https://twiki.cern.ch/twiki/bin/view/LHCPhysics/
* CERNYellowReportPageAt8TeV#gluon_gluon_Fusion_Process
      DATA XSM/36.32d0,29.68d0,20.8d0,15.42d0,11.96d0,8.98d0,7.081d0,
     . 4.783d0,3.594d0,3.401d0,2.921d0,2.002d0,1.283d0,.523d0,.229d0,
     . .1097d0,.0571d0,.032d0/

* SM Higgs bbH prod. cross sect. at 8 TeV from 
* https://twiki.cern.ch/twiki/bin/view/LHCPhysics/
* /CrossSectionsFigures#MSSM_WG_plots (estimated)
      DATA XSMbb/0.56d0,0.42d0,0.25d0,0.15d0,8.7d-2,5.1d-2,3.2d-2,
     . 1.2d-2,4.5d-3,2.8d-3,1.4d-3,8.2d-4,4.9d-4,2.6d-4,1.2d-4,
     . 5.2d-5,2.4d-5,1.2d-5/

* Upper limit on ggF->H->tautau (8 TeV) from CMS-PAS-HIG-13-021, Table 7
      DATA LCMS/50.2d0,31.3d0,7.38d0,2.27d0,.845d0,.549d0,.517d0,.315d0,
     . .15d0,.112d0,.103d0,.607d-1,.385d-1,.193d-1,.143d-1,.115d-1,
     . .923d-2,.865d-2/

* Upper limit on bbH->tautau (8 TeV) from CMS-PAS-HIG-13-021, Table 8
      DATA LCMSbb/6.03d0,4.14d0,1.76d0,1.25d0,.814d0,.659d0,.553d0,
     . .217d0,
     . .975d-1,.638d-1,.613d-1,.431d-1,.320d-1,.203d-1,.173d-1,.166d-1,
     . .146d-1,.133d-1/

* Upper limit on ggF->H->tautau (8 TeV) from ATLAS-CONF-2014-049, Fig. 7
      DATA LATLASgg/29.1d0,24.0d0,5.25d0,2.02d0,1.39d0,1.00d0,.794d0,
     .  .281d0,.127d0,.112d0,.773d-1,.400d-1,.240d-1,.177d-1,.127d-1,
     .  .993d-2,.840d-2,.735d-2/

* Upper limit on bbH->tautau (8 TeV) from ATLAS-CONF-2014-049, Fig. 7
      DATA LATLASbb/6.32d0,6.32d0,2.73d0,1.27d0,.966d0,.606d0,.393d0,
     .  .305d0,.116d0,.101d0,.656d-1,.363d-1,.238d-1,.159d-1,.117d-1,
     .  .943d-2,.785d-2,.716d-2/

      DATA X1/.25d0,.5d0,.75d0,1d0,2d0,3d0/
      DATA LSIGBR/3.8d-5,3.6d-5,4.1d-5,4.2d-5,4.4d-5,4.6d-5/

      COMMON/HISPEC/MHC,XC,MH0,XH,MA2
      COMMON/HNSMCOUP/CU,CUP,CD,CDP,CB,CBP,CJ,CJP,CI,CG,CGP,CV,CZG,CZGP
      COMMON/HNSMBR/BRJJ,BREE,BRMM,BRLL,BRSS,BRCC,BRBB,BRTT,BRWW,BRZZ,
     . BRGG,BRZG
      COMMON/HNHIBR/BRHHH,BRHCHC,BRHAZ,BRHCW,BRHIGGS
      COMMON/HNSUSYBR/BRNEU,BRCHA,BRHSQ,BRHSL,BRSUSY
      COMMON/BR_top2body/brtopbw,brtopbh,brtopneutrstop
      COMMON/HCSMBR/HCBRM,HCBRL,HCBRSU,HCBRBU,HCBRSC,HCBRBC,HCBRBT
      COMMON/HIGGSFIT/MHmin,MHmax,chi2max,chi2gam,chi2bb,chi2zz
      COMMON/LHCSIGCPV/SIG
      COMMON/HIGGSMS/muH2

***********************************************************************

c       I- Constraints from ggF/bb->H/A->tautau from CMS-PAS-HIG-13-021

* Loop over 5 Higgses
      DO I=1,5
        MH(I)=dsqrt(MH0(I))
        J=1
        DOWHILE(HMAS(J).LE.MH(I) .AND. J.LT.NX)
          J=J+1
        ENDDO
        IF(J.GE.2 .AND. MH(I).LT.HMAS(NX)) THEN
        XSMH(I)=0d0
      LCMSH(I)=1d10
        XSMHbb(I)=0d0
      LCMSHbb(I)=1d10

* SM Higgs ggF prod. cross sect.:
          XSMH(I)=XSM(J-1)+(MH(I)-HMAS(J-1))/
     .      (HMAS(J)-HMAS(J-1))*(XSM(J)-XSM(J-1))

* ggF Signal cross section*BR:
          SIG1(I)=(CJ(I)**2+CJP(I)**2)*BRLL(I)*XSMH(I)

* CMS ggF limit:
          LCMSH(I)=LCMS(J-1)+(MH(I)-HMAS(J-1))/
     .      (HMAS(J)-HMAS(J-1))*(LCMS(J)-LCMS(J-1))
          PROB(51)=PROB(51)+DDIM(1d0,LCMSH(I)/SIG1(I))

* ATLAS ggF limit:
          LATLASH(I)=LATLASgg(J-1)+(MH(I)-HMAS(J-1))/
     .      (HMAS(J)-HMAS(J-1))*(LATLASgg(J)-LATLASgg(J-1))

* Correct for jump in Fig.7 at MA=200 GeV: J=8, 
* modif. LATLASgg(J-1)=LATLASgg(7)=.96D0 and not .794d0:
        IF(J.EQ.8) THEN
            LATLASH(I)=.96D0+(MH(I)-HMAS(J-1))/
     .        (HMAS(J)-HMAS(J-1))*(LATLASgg(J)-.96D0)
        ENDIF
          PROB(51)=PROB(51)+DDIM(1D0,LATLASH(I)/SIG1(I))
***

* SM Higgs bbH prod. cross sect.:
          XSMHbb(I)=XSMbb(J-1)+(MH(I)-HMAS(J-1))/
     .      (HMAS(J)-HMAS(J-1))*(XSMbb(J)-XSMbb(J-1))
* bbH Signal cross section*BR:
          SIGbb(I)=(CB(I)**2+CBP(I)**2)*XSMHbb(I)*BRLL(I)
* CMS Hbb limit:
          LCMSHbb(I)=LCMSbb(J-1)+(MH(I)-HMAS(J-1))/
     .      (HMAS(J)-HMAS(J-1))*(LCMSbb(J)-LCMSbb(J-1))
          PROB(51)=PROB(51)+DDIM(1d0,LCMSHbb(I)/SIGbb(I))
* ATLAS Hbb limit:
          LATLASHbb(I)=LATLASbb(J-1)+(MH(I)-HMAS(J-1))/
     .      (HMAS(J)-HMAS(J-1))*(LATLASbb(J)-LATLASbb(J-1))
* Correct for jump in Fig.7 at MA=200 GeV: J=8, 
* modif. LATLASbb(J-1)=LATLASbb(7)=.858D0 and not .393d0:
        IF(J.EQ.8) THEN
            LATLASHbb(I)=.858d0+(MH(I)-HMAS(J-1))/
     .        (HMAS(J)-HMAS(J-1))*(LATLASbb(J)-.858d0)
        ENDIF
          PROB(51)=PROB(51)+DDIM(1d0,LATLASHbb(I)/SIGbb(I))

*********************************************************
* Combine signal rates of 2 Higgses
          DO I1=1,I-1
            J1=1
            DOWHILE(HMAS(J1).LE.MH(I1) .AND. J1.LT.NX)
              J1=J1+1
            ENDDO
          IF(MH(I1).GE.HMAS(1) .AND. MH(I1).LT.HMAS(NX)) THEN

* Average masses weighted by the signal rates (MBAR for ggF, MBARbb for bbH):
             MBAR=(SIG1(I)*MH(I)+SIG1(I1)*MH(I1))/(SIG1(I)+SIG1(I1))
             JBAR=1
             DOWHILE(HMAS(JBAR).LE.MBAR.AND.JBAR.LT.NX)
               JBAR=JBAR+1
             ENDDO
             MBARbb=(SIGbb(I)*MH(I)
     .         +SIGbb(I1)*MH(I1))/(SIGbb(I)+SIGbb(I1))
             JBARbb=1
             DOWHILE(HMAS(JBARbb).LE.MBARbb.AND.JBARbb.LT.NX)
               JBARbb=JBARbb+1
             ENDDO
* DEL=mass difference divided by a (small) resolution squared:
* [DEL < 1 only if |MH(I)-MH(I1)| < (MH(I)+MH(I1))/15d0;
*  otherwise the combined signal rate is small]
             DEL=((MH(I)-MH(I1))/(MH(I)+MH(I1))*15d0)**2
****

* Estimate of the combined ggF signal rates:
             SIGTOT=SIG1(I)+SIG1(I1)
     .         -SIG1(I)*SIG1(I1)*DEL/(SIG1(I)+SIG1(I1))
* Continue only if SIGTOT > 0 and 90<MBAR<1000
*      and |MH(I)-MH(I1)|/MBAR<0.20:
             IF(SIGTOT.GT.0D0.AND.MBAR.GE.HMAS(1).AND.
     .         MBAR.LT.HMAS(NX).AND.
     .           dabs(MH(I)-MH(I1))/MBAR.LE.0.20d0) THEN
* CMS ggF limit at MBAR:
               LCMSMB=LCMS(JBAR-1)+(MBAR-HMAS(JBAR-1))/
     .          (HMAS(JBAR)-HMAS(JBAR-1))*(LCMS(JBAR)-LCMS(JBAR-1))
               PROB(51)=PROB(51)+DDIM(1d0,LCMSMB/SIGTOT)
* ATLAS ggF limit at MBAR:
               LATLASMB=LATLASgg(JBAR-1)+(MBAR-HMAS(JBAR-1))/
     .       (HMAS(JBAR)-HMAS(JBAR-1))*(LATLASgg(JBAR)-LATLASgg(JBAR-1))
* Correct for jump in Fig.7 at MA=200 GeV: JBAR=8, 
* modif. LATLASgg(JBAR-1)=LATLASgg(7)=.96D0 and not .794d0:
               IF(J.EQ.8) THEN
                   LATLASMB=.96D0+(MBAR-HMAS(JBAR-1))/
     .             (HMAS(JBAR)-HMAS(JBAR-1))*(LATLASgg(JBAR)-.96D0)
               ENDIF
               PROB(51)=PROB(51)+DDIM(1D0,LATLASMB/SIGTOT)
             ENDIF
****

* Estimate of the combined bbH signal rates:
             SIGTOTbb=SIGbb(I)+SIGbb(I1)
     .         -SIGbb(I)*SIGbb(I1)*DEL/(SIGbb(I)+SIGbb(I1))
* Continue only if SIGTOTbb > 0 and 90<MBARbb<1000
*      and |MH(I)-MH(I1)|/MBARbb<0.20:
             IF(SIGTOTbb.GT.0D0.AND.MBARbb.GE.HMAS(1).AND.
     .           MBARbb.LT.HMAS(NX).AND.
     .             dabs(MH(I)-MH(I1))/MBARbb.LE.0.20d0) THEN
* CMS bbH limit at MBARbb:
               LCMSMBbb=LCMSbb(JBARbb-1)+(MBARbb-HMAS(JBARbb-1))/
     .           (HMAS(JBARbb)-HMAS(JBARbb-1))*
     .             (LCMSbb(JBARbb)-LCMSbb(JBARbb-1))
               PROB(51)=PROB(51)+DDIM(1d0,LCMSMBbb/SIGTOTbb)
* ATLAS bbH limit at MBARbb:
               LATLASMBbb=LATLASbb(JBARbb-1)+(MBARbb-HMAS(JBARbb-1))/
     .           (HMAS(JBARbb)-HMAS(JBARbb-1))*
     .             (LATLASbb(JBARbb)-LATLASbb(JBARbb-1))
* Correct for jump in Fig.7 at MA=200 GeV: JBARbb=8, 
* modif. LATLASbb(JBARbb-1)=LATLASbb(7)=.858D0 and not .393d0:
              IF(J.EQ.8) THEN
                  LATLASMBbb=.858D0+(MBARbb-HMAS(JBARbb-1))/
     .              (HMAS(JBARbb)-HMAS(JBARbb-1))*
     .                (LATLASbb(JBARbb)-.858D0)
              ENDIF

               PROB(51)=PROB(51)+DDIM(1D0,LATLASMBbb/SIGTOTbb)
             ENDIF
           ENDIF
           ENDDO
         ENDIF
!        write(*,*) "Prob(51):",prob(51)
      ENDDO
!      write(*,*) "Prob(51):",prob(51)


c        II- Bound on Br(t->bH+)*BR(H+->tau nu)

      PROB(45)=DDIM(brtopbh*HCBRL/LHC_TBH(dsqrt(MHC)),1d0)


c        III- Constraints from ggF->H/A->gamgam from ATLAS-CONF-2014-031, M_H/A < 122

*   The EXPCON_PATH variable is set:
      CALL getenv('EXPCON_PATH',EXPCON_PATH)
      if(EXPCON_PATH.eq.' ')  EXPCON_PATH='../EXPCON'

* Read ATLAS upper limit
* ggHgg(I,1): Higgs mass
* ggHgg(I,2): upper limit in fb
      FILENAME=catpath(EXPCON_PATH,'ggHgg.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1801 READ(11,*,END=1802,ERR=2)(dummy(I,J),J=1,4)
      ggHgg(I,1)=dummy(I,1)
      ggHgg(I,2)=dummy(I,4)
      I=I+1
      GOTO 1801
 1802 CLOSE(11)

* Read SM Higgs ggF production cross section (60-122 GeV)
* SMXS(I,1): Higgs mass
* SMXS.65-122.dat:  ggF production cross section in pb
* SMXS(I,2): ggF production cross section in fb
      FILENAME=catpath(EXPCON_PATH,'SMXS.65-122.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1811 READ(11,*,END=1812,ERR=2)(SMXS(I,J),J=1,2)
      SMXS(I,2)=1d3*SMXS(I,2)
      I=I+1
      GOTO 1811
 1812 CLOSE(11)

* Loop over 5 Higgses
      DO I=1,5
        J=1
        DOWHILE(SMXS(J,1).LE.dsqrt(MH0(I)) .AND. J.LT.70)
          J=J+1
        ENDDO

       IF(J.GE.2 .AND. dsqrt(MH0(I)).LT.122d0) THEN
        XSMH(I)=0d0
        LATLASH(I)=1d10
* SM Higgs ggF prod. cross sect.:
          XSMH(I)=SMXS(J-1,2)+(dsqrt(MH0(I))-SMXS(J-1,1))/
     .      (SMXS(J,1)-SMXS(J-1,1))*(SMXS(J,2)-SMXS(J-1,2))
* ggF Signal cross section*BR(H->gamgam):
          SIG(I,8)=(CJ(I)**2+CJP(I)**2)*BRGG(I)*XSMH(I)
* ATLAS limit:
          LATLASH(I)=ggHgg(J-1,2)+(dsqrt(MH0(I))-SMXS(J-1,1))/
     .      (SMXS(J,1)-SMXS(J-1,1))*(ggHgg(J,2)-ggHgg(J-1,2))
          PROB(53)=PROB(53)+DDIM(1d0,LATLASH(I)/SIG(I,8))

        ENDIF
        ENDDO


c       IV- Higgs test at ~125 GeV

c      1) Data

* adding linearly 1 GeV exp. + 2 GeV theor. errors
      MHmin=DSQRT(muH2)-3d0
      MHmax=DSQRT(muH2)+3d0
      chi2max=6.18d0

* From J. Bernon with run-II data (sep. 2016):

c Chi^2 from gammagamma:
      agg=17.47d0
      bgg=3.17d0
      cgg=6.26d0
      mugcengg=1.18d0
      muvcengg=1.07d0

c Chi^2 from bb/tautau:
      abb=4.81d0
      bbb=2.69d0
      cbb=21.57d0
      mugcenbb=1.27d0
      muvcenbb=0.86d0

c Chi^2 from ZZ/WW:
      azz=36.40d0
      bzz=4.56d0
      czz=8.16d0
      mugcenzz=1.11d0
      muvcenzz=1.37d0

c      2) Higgs signals / SM

      DO I=1,5

       DO J=1,10
        SIG(I,J)=0d0
       ENDDO

       CALL HDECAY(dsqrt(MH0(I)),BRJJSM,BREESM,BRMMSM,BRLLSM,BRSSSM,
     .      BRCCSM,BRBBSM,BRTTSM,BRWWSM,BRZZSM,BRGGSM,BRZGSM)

*   H -> tautau
* VBF/VH
       IF(BRLLSM.NE.0d0)SIG(I,1)=CV(I)**2*BRLL(I)/BRLLSM
* ggF
       IF(BRLLSM.NE.0d0)SIG(I,2)=(CJ(I)**2+CJP(I)**2)*BRLL(I)/BRLLSM
       
*   H -> bb
* VBF/VH
       IF(BRBBSM.NE.0d0)SIG(I,3)=CV(I)**2*BRBB(I)/BRBBSM
* ttH
       IF(BRBBSM.NE.0d0)SIG(I,4)=(CU(I)**2+CUP(I)**2)*BRBB(I)/BRBBSM

*   H -> ZZ/WW
* VBF/VH
       IF(BRZZSM.NE.0d0)SIG(I,5)=CV(I)**2*BRZZ(I)/BRZZSM
* ggF
       IF(BRZZSM.NE.0d0)SIG(I,6)=(CJ(I)**2+CJP(I)**2)*BRZZ(I)/BRZZSM
       
*   H -> gammagamma
* VBF/VH
       IF(BRGGSM.NE.0d0)SIG(I,7)=CV(I)**2*BRGG(I)/BRGGSM
* ggF
       IF(BRGGSM.NE.0d0)SIG(I,8)=(CJ(I)**2+CJP(I)**2)*BRGG(I)/BRGGSM

*   H -> invisible
* VBF/VH
       SIG(I,9)=CV(I)**2*BRNEU(I,1,1)
* ggF
       SIG(I,10)=(CJ(I)**2+CJP(I)**2)*BRNEU(I,1,1)

      ENDDO

c      3) Chi^2 test of the couplings

      MHcen=(MHmin+MHmax)/2d0
      D1=(MHmax-MHmin)/2d0
      masstest=2d0
      K=0
      K1=0
      K2=0
      K3=0
      J=1

      DO I=1,6
       SSIG(I)=0d0
      ENDDO

      DO I=1,5
      masstest=min(masstest,dabs(dsqrt(MH0(I))-MHcen)/D1)
       IF(dabs(dsqrt(MH0(I))-MHcen).le.D1)THEN
        if(K.eq.2)then
         K3=I
         K=3
        endif
        if(K.eq.1)then
         K2=I
         K=2
        endif
        if(K.eq.0)then
         K1=I
         K=1
        endif
       ELSE
        if(K.eq.0.and.I.gt.1)then
         if(dabs(dsqrt(MH0(I))-MHcen)
     .  .lt.dabs(dsqrt(MH0(I-1))-MHcen))J=I
        endif
       ENDIF
      ENDDO

      IF(K.eq.0)THEN
       chi2gam=agg*(SIG(J,1)-mugcengg)**2+cgg*(SIG(J,2)-muvcengg)**2
     .     +2d0*bgg*(SIG(J,1)-mugcengg)*(SIG(J,2)-muvcengg)
       chi2bb=abb*(SIG(J,3)-mugcenbb)**2+cbb*(SIG(J,4)-muvcenbb)**2
     .    +2d0*bbb*(SIG(J,3)-mugcenbb)*(SIG(J,4)-muvcenbb)
       chi2zz=azz*(SIG(J,5)-mugcenzz)**2+czz*(SIG(J,6)-muvcenzz)**2
     .    +2d0*bzz*(SIG(J,5)-mugcenzz)*(SIG(J,6)-muvcenzz)
      ENDIF

      IF(K.eq.1)THEN
       chi2gam=agg*(SIG(K1,1)-mugcengg)**2+cgg*(SIG(K1,2)-muvcengg)**2
     .     +2d0*bgg*(SIG(K1,1)-mugcengg)*(SIG(K1,2)-muvcengg)
       chi2bb=abb*(SIG(K1,3)-mugcenbb)**2+cbb*(SIG(K1,4)-muvcenbb)**2
     .    +2d0*bbb*(SIG(K1,3)-mugcenbb)*(SIG(K1,4)-muvcenbb)
       chi2zz=azz*(SIG(K1,5)-mugcenzz)**2+czz*(SIG(K1,6)-muvcenzz)**2
     .    +2d0*bzz*(SIG(K1,5)-mugcenzz)*(SIG(K1,6)-muvcenzz)
      ENDIF

      IF(K.eq.2)THEN
       chi2gam=min(agg*(SIG(K1,1)-mugcengg)**2
     .                +cgg*(SIG(K1,2)-muvcengg)**2
     .     +2d0*bgg*(SIG(K1,1)-mugcengg)*(SIG(K1,2)-muvcengg),
     .     agg*(SIG(K2,1)-mugcengg)**2+cgg*(SIG(K2,2)-muvcengg)**2
     .     +2d0*bgg*(SIG(K2,1)-mugcengg)*(SIG(K2,2)-muvcengg))
       chi2zz=min(azz*(SIG(K1,5)-mugcenzz)**2
     .                +czz*(SIG(K1,6)-muvcenzz)**2
     .    +2d0*bzz*(SIG(K1,5)-mugcenzz)*(SIG(K1,6)-muvcenzz),
     .      azz*(SIG(K2,5)-mugcenzz)**2+czz*(SIG(K2,6)-muvcenzz)**2
     .    +2d0*bzz*(SIG(K2,5)-mugcenzz)*(SIG(K2,6)-muvcenzz))

        SSIG(1)=SIG(K1,8)+SIG(K2,8)
        SSIG(2)=SIG(K1,7)+SIG(K2,7)
        SSIG(3)=SIG(K1,2)+SIG(K2,2)
        SSIG(4)=SIG(K1,3)+SIG(K2,3)
        SSIG(5)=SIG(K1,6)+SIG(K2,6)
        SSIG(6)=SIG(K1,5)+SIG(K2,5)

        if(dabs(dsqrt(MH0(K1))-dsqrt(MH0(K2))).le.3d0)then
       chi2gam=agg*(SSIG(1)-mugcengg)**2+cgg*(SSIG(2)-muvcengg)**2
     .     +2d0*bgg*(SSIG(1)-mugcengg)*(SSIG(2)-muvcengg)
       chi2zz=azz*(SSIG(5)-mugcenzz)**2+czz*(SSIG(6)-muvcenzz)**2
     .    +2d0*bzz*(SSIG(5)-mugcenzz)*(SSIG(6)-muvcenzz)
        endif

       chi2bb=abb*(SSIG(3)-mugcenbb)**2+cbb*(SSIG(4)-muvcenbb)**2
     .    +2d0*bbb*(SSIG(3)-mugcenbb)*(SSIG(4)-muvcenbb)
      ENDIF

      IF(K.eq.3)THEN
       chi2gam=min(agg*(SIG(K1,1)-mugcengg)**2
     .                 +cgg*(SIG(K1,2)-muvcengg)**2
     .     +2d0*bgg*(SIG(K1,1)-mugcengg)*(SIG(K1,2)-muvcengg),
     .     agg*(SIG(K2,1)-mugcengg)**2+cgg*(SIG(K2,2)-muvcengg)**2
     .     +2d0*bgg*(SIG(K2,1)-mugcengg)*(SIG(K2,2)-muvcengg),
     .     agg*(SIG(K3,1)-mugcengg)**2+cgg*(SIG(K3,2)-muvcengg)**2
     .     +2d0*bgg*(SIG(K3,1)-mugcengg)*(SIG(K3,2)-muvcengg))
       chi2zz=min(azz*(SIG(K1,5)-mugcenzz)**2
     .                 +czz*(SIG(K1,6)-muvcenzz)**2
     .    +2d0*bzz*(SIG(K1,5)-mugcenzz)*(SIG(K1,6)-muvcenzz),
     .      azz*(SIG(K2,5)-mugcenzz)**2+czz*(SIG(K2,6)-muvcenzz)**2
     .    +2d0*bzz*(SIG(K2,5)-mugcenzz)*(SIG(K2,6)-muvcenzz),
     .      azz*(SIG(K3,5)-mugcenzz)**2+czz*(SIG(K3,6)-muvcenzz)**2
     .    +2d0*bzz*(SIG(K3,5)-mugcenzz)*(SIG(K3,6)-muvcenzz))

        SSIG(1)=SIG(K1,8)+SIG(K2,8)+SIG(K3,8)
        SSIG(2)=SIG(K1,7)+SIG(K2,7)+SIG(K3,7)
        SSIG(3)=SIG(K1,2)+SIG(K2,2)+SIG(K3,2)
        SSIG(4)=SIG(K1,3)+SIG(K2,3)+SIG(K3,3)
        SSIG(5)=SIG(K1,6)+SIG(K2,6)+SIG(K3,6)
        SSIG(6)=SIG(K1,5)+SIG(K2,5)+SIG(K3,5)

        if(dabs(dsqrt(MH0(K1))-dsqrt(MH0(K3))).le.3d0)then
       chi2gam=agg*(SSIG(1)-mugcengg)**2+cgg*(SSIG(2)-muvcengg)**2
     .     +2d0*bgg*(SSIG(1)-mugcengg)*(SSIG(2)-muvcengg)
       chi2zz=azz*(SSIG(5)-mugcenzz)**2+czz*(SSIG(6)-muvcenzz)**2
     .    +2d0*bzz*(SSIG(5)-mugcenzz)*(SSIG(6)-muvcenzz)
        endif

       chi2bb=abb*(SSIG(3)-mugcenbb)**2+cbb*(SSIG(4)-muvcenbb)**2
     .    +2d0*bbb*(SSIG(3)-mugcenbb)*(SSIG(4)-muvcenbb)

        if(dabs(dsqrt(MH0(K1))-dsqrt(MH0(K3))).gt.3d0)then
         if(dabs(dsqrt(MH0(K1))-dsqrt(MH0(K2))).lt.3d0.and.
     .       dabs(dsqrt(MH0(K2))-dsqrt(MH0(K3))).gt.3d0)then
         SSIG(1)=SIG(K1,8)+SIG(K2,8)
         SSIG(2)=SIG(K1,7)+SIG(K2,7)
         SSIG(3)=SIG(K1,2)+SIG(K2,2)
         SSIG(4)=SIG(K1,3)+SIG(K2,3)
         SSIG(5)=SIG(K1,6)+SIG(K2,6)
         SSIG(6)=SIG(K1,5)+SIG(K2,5)
         chi2gam=min(agg*(SSIG(1)-mugcengg)**2
     .                       +cgg*(SSIG(2)-muvcengg)**2
     .     +2d0*bgg*(SSIG(1)-mugcengg)*(SSIG(2)-muvcengg),
     .     agg*(SIG(K3,1)-mugcengg)**2+cgg*(SIG(K3,2)-muvcengg)**2
     .     +2d0*bgg*(SIG(K3,1)-mugcengg)*(SIG(K3,2)-muvcengg))
         chi2zz=min(azz*(SSIG(5)-mugcenzz)**2
     .                       +czz*(SSIG(6)-muvcenzz)**2
     .    +2d0*bzz*(SSIG(5)-mugcenzz)*(SSIG(6)-muvcenzz),
     .      azz*(SIG(K3,5)-mugcenzz)**2+czz*(SIG(K3,6)-muvcenzz)**2
     .    +2d0*bzz*(SIG(K3,5)-mugcenzz)*(SIG(K3,6)-muvcenzz))
         endif
         if(dabs(dsqrt(MH0(K1))-dsqrt(MH0(K2))).gt.3d0.and.
     .       dabs(dsqrt(MH0(K2))-dsqrt(MH0(K3))).lt.3d0)then
         SSIG(1)=SIG(K3,8)+SIG(K2,8)
         SSIG(2)=SIG(K3,7)+SIG(K2,7)
         SSIG(3)=SIG(K3,2)+SIG(K2,2)
         SSIG(4)=SIG(K3,3)+SIG(K2,3)
         SSIG(5)=SIG(K3,6)+SIG(K2,6)
         SSIG(6)=SIG(K3,5)+SIG(K2,5)
         chi2gam=min(agg*(SSIG(1)-mugcengg)**2
     .                              +cgg*(SSIG(2)-muvcengg)**2
     .     +2d0*bgg*(SSIG(1)-mugcengg)*(SSIG(2)-muvcengg),
     .     agg*(SIG(K1,1)-mugcengg)**2+cgg*(SIG(K1,2)-muvcengg)**2
     .     +2d0*bgg*(SIG(K1,1)-mugcengg)*(SIG(K1,2)-muvcengg))
         chi2zz=min(azz*(SSIG(5)-mugcenzz)**2
     .                              +czz*(SSIG(6)-muvcenzz)**2
     .    +2d0*bzz*(SSIG(5)-mugcenzz)*(SSIG(6)-muvcenzz),
     .      azz*(SIG(K1,5)-mugcenzz)**2+czz*(SIG(K1,6)-muvcenzz)**2
     .    +2d0*bzz*(SIG(K1,5)-mugcenzz)*(SIG(K1,6)-muvcenzz))
         endif
         if(dabs(dsqrt(MH0(K1))-dsqrt(MH0(K2))).lt.3d0.and.
     .       dabs(dsqrt(MH0(K2))-dsqrt(MH0(K3))).lt.3d0)then
         SSIG(1)=SIG(K1,8)+SIG(K2,8)
         SSIG(2)=SIG(K1,7)+SIG(K2,7)
         SSIG(3)=SIG(K1,2)+SIG(K2,2)
         SSIG(4)=SIG(K1,3)+SIG(K2,3)
         SSIG(5)=SIG(K1,6)+SIG(K2,6)
         SSIG(6)=SIG(K1,5)+SIG(K2,5)
         chi2gam=agg*(SSIG(1)-mugcengg)**2+cgg*(SSIG(2)-muvcengg)**2
     .     +2d0*bgg*(SSIG(1)-mugcengg)*(SSIG(2)-muvcengg)
         chi2zz=azz*(SSIG(5)-mugcenzz)**2+czz*(SSIG(6)-muvcenzz)**2
     .    +2d0*bzz*(SSIG(5)-mugcenzz)*(SSIG(6)-muvcenzz)
         SSIG(1)=SIG(K3,8)+SIG(K2,8)
         SSIG(2)=SIG(K3,7)+SIG(K2,7)
         SSIG(3)=SIG(K3,2)+SIG(K2,2)
         SSIG(4)=SIG(K3,3)+SIG(K2,3)
         SSIG(5)=SIG(K3,6)+SIG(K2,6)
         SSIG(6)=SIG(K3,5)+SIG(K2,5)
         chi2gam=min(agg*(SSIG(1)-mugcengg)**2
     .                  +cgg*(SSIG(2)-muvcengg)**2
     .     +2d0*bgg*(SSIG(1)-mugcengg)*(SSIG(2)-muvcengg),chi2gam)
         chi2zz=min(azz*(SSIG(5)-mugcenzz)**2
     .                  +czz*(SSIG(6)-muvcenzz)**2
     .    +2d0*bzz*(SSIG(5)-mugcenzz)*(SSIG(6)-muvcenzz),chi2zz)
         endif
        endif

      ENDIF

      IF(masstest.gt.1d0)PROB(46)=masstest

      PROB(47)=DDIM(chi2gam/chi2MAX,1d0)
      PROB(48)=DDIM(chi2bb/chi2MAX,1d0)
      PROB(49)=DDIM(chi2zz/chi2MAX,1d0)


c        V- CMS constraints on h[125]->2A->4mu

      PROB(52)=0d0

      DO I1=2,5

      CJ2=0d0
      BRHTOAA=0d0

      IF(dabs(dsqrt(MH0(I1))-MHcen).le.D1)THEN
       CJ2=CJ(I1)**2+CJP(I1)**2
       BRHTOAA=BRHHH(I1,1)
      ENDIF

      I=1
      DOWHILE(dsqrt(MH0(1)).GE.X1(1).AND.X1(I).LE.dsqrt(MH0(1))
     .                                 .AND.I.LT.NY)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND.dsqrt(MH0(1)).LT.X1(NY))THEN
       LIMIT=LSIGBR(I-1)
     . +(dsqrt(MH0(1))-X1(I-1))/(X1(I)-X1(I-1))*(LSIGBR(I)-LSIGBR(I-1))
       PROB(52)=PROB(52)+DDIM(CJ2*BRHTOAA*BRMM(1)**2/LIMIT,1d0)
      ENDIF

      ENDDO

      RETURN

*   Error catch

 1    WRITE(*,*)"Cannot find the file ",FILENAME
      STOP

 2    WRITE(*,*)"Read error in the file ",FILENAME
      STOP

      END


      SUBROUTINE LHC_HSMAA_LEPTONS_CPV(PROB)

*     PROB(52) =/= 0: excluded 

      IMPLICIT NONE

      CHARACTER*256 FILENAME,EXPCON_PATH,catpath

      INTEGER NHAATAUS1,NHAATAUS2,NHAABMU
      INTEGER NHAAMUS1,NHAAMUS2,NHAAMUS3
      INTEGER NHAAMUS4,I,J,K,NM
      PARAMETER(NM=14)

      DOUBLE PRECISION MHC,XC(2,2),MH0(5),XH(5,5),MA2
      DOUBLE PRECISION PROB(*)
      DOUBLE PRECISION D1,D2,CJ2BRHTOAA,MH,MHcen,LIMIT
      DOUBLE PRECISION MHSM(NM),SMXS_8TeV(NM),SMXS_125_8TeV,SMXS_J_8TeV
      DOUBLE PRECISION LOWBOUND,HIGHBOUND

      DOUBLE PRECISION BRJJ(5),BREE(5),BRMM(5),BRLL(5),BRSS(5),
     . BRCC(5),BRBB(5),BRTT(5),BRWW(5),BRZZ(5),BRGG(5),BRZG(5)
      DOUBLE PRECISION BRHHH(5,10),BRHCHC(5),BRHAZ(5,4),BRHCW(5),
     . BRHIGGS(5)
      DOUBLE PRECISION CU(5),CUP(5),CD(5),CDP(5),CB(5),CBP(5),CJ(5),
     . CJP(5),CI(5),CG(5),CGP(5),CV(5),CZG(5),CZGP(5)
      DOUBLE PRECISION MHmin,MHmax,chi2max,chi2gam,chi2bb,chi2zz
      DOUBLE PRECISION HAATAUS1(92,2),HAATAUS2(92,2),HAABMU(76,2)
      DOUBLE PRECISION HAAMUS1(20,2),HAAMUS2(12,2),HAAMUS3(7,2)
      DOUBLE PRECISION HAAMUS4(67,4)

      COMMON/HISPEC/MHC,XC,MH0,XH,MA2
      COMMON/HNSMBR/BRJJ,BREE,BRMM,BRLL,BRSS,BRCC,BRBB,BRTT,BRWW,BRZZ,
     . BRGG,BRZG
      COMMON/HNHIBR/BRHHH,BRHCHC,BRHAZ,BRHCW,BRHIGGS
      COMMON/HNSMCOUP/CU,CUP,CD,CDP,CB,CBP,CJ,CJP,CI,CG,CGP,CV,CZG,CZGP
      COMMON/HIGGSFIT/MHmin,MHmax,chi2max,chi2gam,chi2bb,chi2zz
      COMMON/LHCHAA/HAATAUS1,HAATAUS2,HAABMU,
     .      HAAMUS1,HAAMUS2,HAAMUS3,HAAMUS4,
     .      NHAATAUS1,NHAATAUS2,NHAABMU,
     .      NHAAMUS1,NHAAMUS2,NHAAMUS3,NHAAMUS4

      DATA MHSM/85d0,90d0,95d0,100d0,105d0,110d0,115d0,120d0,125d0,
     . 130d0,135d0,140d0,145d0,150d0/
      DATA SMXS_8TeV/3.940d1,3.526d1,3.175d1,2.873d1,2.611d1,2.383d1,
     . 2.184d1,2.008d1,1.851d1,1.712d1,1.587d1,1.475d1,1.375d1,1.284d1/           ! in pb

* Determining the SM-like Higgs and its couplings / BR / XS

      CJ2BRHTOAA=0d0
      MHcen=(MHmin+MHmax)/2d0
      MH=0d0
      D2=0d0

      DO I=2,5
       D1=DDIM(dsqrt(MH0(I))/MHMAX,1d0)-DDIM(1d0,dsqrt(MH0(I))/MHMIN)
       IF(D1.EQ.0d0)THEN
        CJ2BRHTOAA=CJ2BRHTOAA+(CJ(I)**2+CJP(I)**2)*BRHHH(I,1)
        MH=MH+(CJ(I)**2+CJP(I)**2)*dsqrt(MH0(I))
        D2=D2+CJ(I)**2+CJP(I)**2
       ENDIF
      ENDDO
      MH=MH/Max(D2,1d-10)
      D1=DDIM(MH/MHMAX,1d0)-DDIM(1d0,MH/MHMIN)
      IF(D1.ne.0d0)THEN
       MH=MHcen
      ENDIF

      SMXS_125_8TeV=0d0
      I=1
      DOWHILE(MH.GE.MHSM(1) .AND. MHSM(I).LE.MH .AND. I.LT.NM)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. MH.LT.MHSM(NM))THEN
       SMXS_125_8TeV=SMXS_8TeV(I-1)
     .    +(MH-MHSM(I-1))/(MHSM(I)-MHSM(I-1))
     .         *(SMXS_8TeV(I)-SMXS_8TeV(I-1))
      ENDIF

* Constraints from ggF->HSM->AA->2mu2tau from 1505.01609 (ATLAS), M_A < 50GeV
* BR(A->2mu) converted in BR(A->2tau), normalized to SM XS

      I=1
      DOWHILE(dsqrt(MH0(1)).GE.HAATAUS1(1,1)
     .        .AND. HAATAUS1(I,1).LE.dsqrt(MH0(1))
     .        .AND. I.LT.NHAATAUS1)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. dsqrt(MH0(1)).LT.HAATAUS1(NHAATAUS1,1))THEN
       LIMIT=HAATAUS1(I-1,2)
     .    +(dsqrt(MH0(1))-HAATAUS1(I-1,1))
     .         /(HAATAUS1(I,1)-HAATAUS1(I-1,1))
     .         *(HAATAUS1(I,2)-HAATAUS1(I-1,2))
       PROB(52)=PROB(52)+DDIM(CJ2BRHTOAA*BRLL(1)**2/LIMIT,1d0)
      ENDIF


* Constraints from ggF->HSM->AA->4tau from 1510.06534 (CMS), 4GeV < M_A < 8GeV

      I=1
      DOWHILE(dsqrt(MH0(1)).GE.HAATAUS2(1,1)
     .        .AND. HAATAUS2(I,1).LE.dsqrt(MH0(1))
     .        .AND. I.LT.NHAATAUS2)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. dsqrt(MH0(1)).LT.HAATAUS2(NHAATAUS2,1))THEN
       LIMIT=HAATAUS2(I-1,2)
     .    +(dsqrt(MH0(1))-HAATAUS2(I-1,1))
     .         /(HAATAUS2(I,1)-HAATAUS2(I-1,1))
     .         *(HAATAUS2(I,2)-HAATAUS2(I-1,2))
       PROB(52)=PROB(52)
     .          +DDIM(SMXS_125_8TeV*CJ2BRHTOAA*BRLL(1)**2/LIMIT,1d0)
      ENDIF


* Constraints from ggF->HSM->AA->2mu2b from CMS-PAS-HIG-14-041, 25GeV < M_A < 63GeV
* normalized to SM XS

      I=1
      DOWHILE(dsqrt(MH0(1)).GE.HAABMU(1,1)
     .        .AND. HAABMU(I,1).LE.dsqrt(MH0(1))
     .        .AND. I.LT.NHAABMU)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. dsqrt(MH0(1)).LT.HAABMU(NHAABMU,1))THEN
       LIMIT=HAABMU(I-1,2)
     .    +(dsqrt(MH0(1))-HAABMU(I-1,1))
     .         /(HAABMU(I,1)-HAABMU(I-1,1))
     .         *(HAABMU(I,2)-HAABMU(I-1,2))
       PROB(52)=PROB(52)
     .          +DDIM(CJ2BRHTOAA*BRMM(1)*BRBB(1)/LIMIT,1d0)
      ENDIF


* Constraints from ggF->HSM->AA->2mu2tau from CMS-PAS-HIG-15-011, 19GeV < M_A < 57GeV
* BR(A->2tau) converted in BR(A->2mu), normalized to SM XS

      I=1
      DOWHILE(dsqrt(MH0(1)).GE.HAAMUS1(1,1)
     .        .AND. HAAMUS1(I,1).LE.dsqrt(MH0(1))
     .        .AND. I.LT.NHAAMUS1)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. dsqrt(MH0(1)).LT.HAAMUS1(NHAAMUS1,1))THEN
       LIMIT=HAAMUS1(I-1,2)
     .    +(dsqrt(MH0(1))-HAAMUS1(I-1,1))
     .         /(HAAMUS1(I,1)-HAAMUS1(I-1,1))
     .         *(HAAMUS1(I,2)-HAAMUS1(I-1,2))
       PROB(52)=PROB(52)
     .          +DDIM(CJ2BRHTOAA*BRMM(1)**2/LIMIT,1d0)
      ENDIF


* Constraints from ggF->HSM->AA->4tau from CMS-PAS-HIG-14-022, 5GeV < M_A < 14GeV
* BR(A->2tau) converted in BR(A->2mu), as shown in Fig7 of CMS-PAS-HIG-15-011,
* normalized to SM XS

      I=1
      DOWHILE(dsqrt(MH0(1)).GE.HAAMUS2(1,1)
     .        .AND. HAAMUS2(I,1).LE.dsqrt(MH0(1))
     .        .AND. I.LT.NHAAMUS2)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. dsqrt(MH0(1)).LT.HAAMUS2(NHAAMUS2,1))THEN
       LIMIT=HAAMUS2(I-1,2)
     .    +(dsqrt(MH0(1))-HAAMUS2(I-1,1))
     .         /(HAAMUS2(I,1)-HAAMUS2(I-1,1))
     .         *(HAAMUS2(I,2)-HAAMUS2(I-1,2))
       PROB(52)=PROB(52)
     .          +DDIM(CJ2BRHTOAA*BRMM(1)**2/LIMIT,1d0)
      ENDIF


* Constraints from ggF->HSM->AA->4tau from CMS-PAS-HIG-14-019, 3GeV < M_A < 8GeV
* BR(A->2tau) converted in BR(A->2mu), as shown in Fig7 of CMS-PAS-HIG-15-011,
* normalized to SM XS

      I=1
      DOWHILE(dsqrt(MH0(1)).GE.HAAMUS3(1,1)
     .        .AND. HAAMUS3(I,1).LE.dsqrt(MH0(1))
     .        .AND. I.LT.NHAAMUS3)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. dsqrt(MH0(1)).LT.HAAMUS3(NHAAMUS3,1))THEN
       LIMIT=HAAMUS3(I-1,2)
     .    +(dsqrt(MH0(1))-HAAMUS3(I-1,1))
     .         /(HAAMUS3(I,1)-HAAMUS3(I-1,1))
     .         *(HAAMUS3(I,2)-HAAMUS3(I-1,2))
       PROB(52)=PROB(52)
     .          +DDIM(CJ2BRHTOAA*BRMM(1)**2/LIMIT,1d0)
      ENDIF


* Constraints from ggF->H->AA->4mu from 1506.00424, 0.25GeV < M_A < 3.55GeV, 85GeV < m_H < 150GeV

      CJ2BRHTOAA=0d0

        DO J=2,5

      IF(J.ge.3 .and. dabs(dsqrt(MH0(J))-dsqrt(MH0(J-1))).lt.3d0)THEN
       IF(J.ge.4 .and. dabs(dsqrt(MH0(J))-dsqrt(MH0(J-2))).lt.3d0)THEN
        MH=((CJ(J-2)**2+CJP(J-2)**2)*dsqrt(MH0(J-2))
     .    +(CJ(J-1)**2+CJP(J-1)**2)*dsqrt(MH0(J-1))
     .    +(CJ(J)**2+CJP(J)**2)*dsqrt(MH0(J)))
     .  /Max(CJ(J)**2+CJ(J-1)**2+CJ(J-2)**2+CJP(J)**2+CJP(J-1)**2
     .                                           +CJP(J-2)**2,1d-10)
        CJ2BRHTOAA=CJ2BRHTOAA+(CJ(J)**2+CJP(J)**2)*BRHHH(J,1)
       ELSE
        MH=((CJ(J-1)**2+CJP(J-1)**2)*dsqrt(MH0(J-1))
     .    +(CJ(J)**2+CJP(J)**2)*dsqrt(MH0(J)))
     .            /Max(CJ(J)**2+CJ(J-1)**2+CJP(J)**2+CJP(J-1)**2,1d-10)
        CJ2BRHTOAA=CJ2BRHTOAA+(CJ(J)**2+CJP(J)**2)*BRHHH(J,1)
       ENDIF
      ELSE
       MH=dsqrt(MH0(J))
       CJ2BRHTOAA=(CJ(J)**2+CJP(J)**2)*BRHHH(J,1)
      ENDIF

      I=1
      DOWHILE(MH.GE.MHSM(1)
     .        .AND.MHSM(I).LE.MH.AND.I.LT.NM)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. MH.LT.MHSM(14))THEN
       SMXS_J_8TeV=SMXS_8TeV(I-1)
     .    +(MH-MHSM(I-1))/(MHSM(I)-MHSM(I-1))
     .         *(SMXS_8TeV(I)-SMXS_8TeV(I-1))
       SMXS_J_8TeV=SMXS_J_8TeV*1d3         ! converting in fb
      ENDIF


      I=1
      DOWHILE(MH.GE.HAAMUS4(1,1)
     .        .AND. HAAMUS4(I,1).LE.MH .AND. I.LT.NHAAMUS4)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. MH.LT.HAAMUS4(NHAAMUS4,1))THEN
       LOWBOUND=100d0
       HIGHBOUND=100d0
       IF(dsqrt(MH0(1)).ge.0.25d0 .and. dsqrt(MH0(1)).lt.2d0)THEN
        LOWBOUND=HAAMUS4(I-1,2)+(dsqrt(MH0(1))-0.25d0)/(2d0-0.25d0)
     .         *(HAAMUS4(I-1,3)-HAAMUS4(I-1,2))
        HIGHBOUND=HAAMUS4(I,2)+(dsqrt(MH0(1))-0.25d0)/(2d0-0.25d0)
     .         *(HAAMUS4(I,3)-HAAMUS4(I,2))
       ELSEIF(dsqrt(MH0(1)).ge.2d0 .and. dsqrt(MH0(1)).lt.3.55d0)THEN
        LOWBOUND=HAAMUS4(I-1,3)+(dsqrt(MH0(1))-2d0)/(3.55d0-2d0)
     .         *(HAAMUS4(I-1,4)-HAAMUS4(I-1,3))
        HIGHBOUND=HAAMUS4(I,3)+(dsqrt(MH0(1))-2d0)/(3.55d0-2d0)
     .         *(HAAMUS4(I,4)-HAAMUS4(I,3))
       ENDIF

       LIMIT=LOWBOUND
     .    +(MH-HAAMUS4(I-1,1))/(HAAMUS4(I,1)-HAAMUS4(I-1,1))
     .         *(HIGHBOUND-LOWBOUND)
       PROB(52)=PROB(52)
     .          +DDIM(SMXS_J_8TeV*CJ2BRHTOAA*BRMM(1)**2/LIMIT,1d0)
      ENDIF

        ENDDO

      RETURN

      END
