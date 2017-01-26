      SUBROUTINE runpar_CPV(PAR,IFAIL)

      IMPLICIT NONE

      CHARACTER*256 FILENAME,EXPCON_PATH,catpath

      INTEGER IFAIL,OMGFLAG,MAFLAG,MOFLAG,Q2FIX

      DOUBLE PRECISION PAR(*),MYPHASES(16),fSF1
      DOUBLE PRECISION Yt,Yb,Ytau,alsmt,runmb
      DOUBLE PRECISION mu,M1,M2,M3,mhs2,mas2
      DOUBLE PRECISION QSUSY,QSL,QSQ,l,k,IAL,IAk,IXIS
      DOUBLE PRECISION Pi,EPS,NMB0,NMB1,MA2,MP2
      DOUBLE PRECISION ALSMA,DLA,DLQA,F1,F2,HTMA

      DOUBLE PRECISION Q2MIN,Q2,QSTSB
      DOUBLE PRECISION G1Q,G2Q,GQ,ALSQ
      DOUBLE PRECISION ZHU,ZHD,ZS,vuq,vdq,TANBQ
      DOUBLE PRECISION XIF,XIS,MUP,MSP,M3H
      DOUBLE PRECISION XIFQ,XISQ,MUPQ,MSPQ,M3HQ
      DOUBLE PRECISION lq,kq,Alcos1,Akcos2,muq,nuq
      DOUBLE PRECISION Ytq,Ybq,MTOPQ,MBOTQ
      DOUBLE PRECISION PhiAT,PhiAB,aux1,aux2

      DOUBLE PRECISION AT,AB,ANOMQSTSB,RTOP,RBOT
      DOUBLE PRECISION MQ3,MU3,MD3,ML3,ME3,MQ,MU1,MD,ML,ME
      DOUBLE PRECISION XI,HT2,HB2,L2,K2,HTAU2,MH1,MH2
      DOUBLE PRECISION SP,SIG1,SIG2,SIG3,F20,F21,F22
      DOUBLE PRECISION MHuS,MHdS,MSS

      DOUBLE PRECISION muH2
      DOUBLE PRECISION GF,MZ,MW,g1,g2,alSMZ,S2TW,ALEMMZ
      DOUBLE PRECISION mt,mb,mtau,mmu,mel,MS,MC,MBP,MPI,MSTRANGE
      DOUBLE PRECISION tanb,cosb,sinb,vu,vd
      DOUBLE PRECISION ALSMZC,ALEMMZC,GFC,g1C,g2C,S2TWC
      DOUBLE PRECISION MSC,MCC,MBC,MBPC,MTC,MTAUC,MMUC,MZC,MWC
      DOUBLE PRECISION MPIC,MELC,MSTRANGEC
      DOUBLE PRECISION phi01,phi02,phi0,phiM1,phiM2,phiM3,
     .              phiATQ,phiABQ,phiATAU,phiAC,phiAS,phiAMU
      DOUBLE PRECISION mur,M1r,M2r,msi
      DOUBLE PRECISION MSL3,MSE3,MSL1,MSE1,ATAU,AMU
      DOUBLE PRECISION MSQ1,MSU1,MSD1
      DOUBLE PRECISION Drv
      DOUBLE PRECISION phiF,phiP,phi3,phiS,phiSP,phi3q,phiSq,phiSPq,
     .              mupsi,ks2si
      DOUBLE PRECISION MSQ3,MSU3,MSD3,ATP,ABP
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS

      COMMON/HIGGSMS/muH2
      COMMON/RENSCALE/Q2
      COMMON/STSBSCALE/QSTSB
      COMMON/Q2FIX/Q2MIN,Q2FIX
      COMMON/QGAUGE/G1Q,G2Q,GQ,ALSQ
      COMMON/QHIGGS/ZHU,ZHD,ZS,vuq,vdq,TANBQ
      COMMON/SUSYEXT/XIF,XIS,MUP,MSP,M3H
      COMMON/QEXT/XIFQ,XISQ,MUPQ,MSPQ,M3HQ
      COMMON/QPAR/lq,kq,Alcos1,Akcos2,muq,NUQ
      COMMON/QQUARK/Ytq,Ybq,MTOPQ,MBOTQ
      COMMON/EWPAR/GF,MZ,MW,g1,g2,alSMZ,S2TW,ALEMMZ
      COMMON/SMFERM/mt,mb,mtau,mmu,mel,MS,MC,MBP,MPI,MSTRANGE
      COMMON/GAUGE/ALSMZC,ALEMMZC,GFC,g1C,g2C,S2TWC
      COMMON/SMSPEC/MSC,MCC,MBC,MBPC,MTC,MTAUC,MMUC,MZC,MWC
      COMMON/SMEXT/MPIC,MELC,MSTRANGEC
      COMMON/TBPAR/tanb,cosb,sinb,vu,vd
      COMMON/PHASES/phi01,phi02,phi0,phiM1,phiM2,phiM3,
     .              phiATQ,phiABQ,phiATAU,phiAC,phiAS,phiAMU
      COMMON/GAUGINOPAR/mur,M1r,M2r,msi
      COMMON/SLEPPAR/MSL3,MSE3,MSL1,MSE1,ATAU,AMU
      COMMON/SQUPAR/MSQ1,MSU1,MSD1
      COMMON/VEVSHIFT/Drv
      COMMON/FLAGS/OMGFLAG,MAFLAG,MOFLAG
      COMMON/Z3VAUX/phiF,phiP,phi3,phiS,phiSP,phi3q,phiSq,phiSPq,
     .              mupsi,ks2si
      COMMON/MYPHASES/MYPHASES
      COMMON/MH2TREE/MHuS,MHdS,MSS
      COMMON/IMALAK/IAL,IAK,IXIS
      COMMON/RADCOR2/MSQ3,MSU3,MSD3,ATP,ABP
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS


      EPS=1d-16
      PI=4d0*DATAN(1d0)

c      1) Gauge (at MZ)
      GF=GFC
      MZ=MZC
      MW=MWC
      g1=g1C
      g2=g2C
      ALSMZ=ALSMZC
      S2TW=S2TWC
      ALEMMZ=ALEMMZC

c      2) tanb
      tanb=PAR(3)
      cosb=1d0/dsqrt(1d0+tanb**2)
      sinb=tanb*cosb
      Drv=0d0
      vu=sinb/dsqrt(2d0*dsqrt(2d0)*GF*(1d0-Drv))
      vd=cosb/dsqrt(2d0*dsqrt(2d0)*GF*(1d0-Drv))

c      3) Phases
      phi01=MYPHASES(1) !phil+phi1
      phi02=MYPHASES(2) !phik+phi2
      phi0=phi01-phi02
      phiM1=MYPHASES(3)
      phiM2=MYPHASES(4)
      phiM3=MYPHASES(5)
      phiAT=MYPHASES(6)
      phiAB=MYPHASES(7)
      phiATAU=MYPHASES(8)
      phiAC=MYPHASES(9)
      phiAS=MYPHASES(10)
      phiAMU=MYPHASES(11)
      IF(PAR(2).ne.0d0)THEN
       phiS=MYPHASES(12)
      ELSE
       phiS=0d0
      ENDIF
      phiSP=MYPHASES(13)
      phiF=MYPHASES(14)
      phiP=MYPHASES(15)
      phi3=MYPHASES(16)

c      4) Yukawa (at mt) and SM fermions
      mt=MTC
      mb=MBC            ! MSbar, scale MB
      mtau=MTAUC
      mmu=MMUC
      mel=MELC

      MS=MSC
      MC=MCC
      MBP=MBPC

      MSTRANGE=MSTRANGEC
      MPI=MPIC
      alsmt=ALSMZ/(1.d0+23.d0/12.d0/Pi*ALSMZ*dlog((mt/MZ)**2))
      Yt=mt/vu
      Yt=Yt/(1.d0+4.d0/3.d0/Pi*alsmt+11.d0/Pi**2*alsmt**2)
c      Yb=mb/vd*(1.d0-23.d0/12.d0/Pi*alsmt
c     c   *dlog(mt**2/MZ**2))**(12d0/23d0)
      Yb=RUNMB(mt)/vd
      Ytau=mtau/vd

c      4) Sfermion parameters
      MSL3=PAR(10)
      MSE3=PAR(11)
      MSL1=PAR(18)
      MSE1=PAR(19)
      ATAU=PAR(14)
      AMU=PAR(25)
      MSQ1=PAR(15)
      MSU1=PAR(16)
      MSD1=PAR(17)

*   Definition of the SUSY scale Q2, unless defined initially
      IF(Q2FIX.EQ.0)THEN
       Q2=MAX((2d0*PAR(15)+PAR(16)+PAR(17))/4d0,Q2MIN)
      ENDIF
*   Definition of the scale QSTSB
      QSTSB=DSQRT(MAX(PAR(7)*PAR(8),Q2MIN**2))

      QSUSY=dsqrt(QSTSB)
      QSL=dsqrt((2.d0*MSL3+MSE3+2.d0*MSE1+4.d0*MSL1)/9.d0)
      QSQ=dsqrt((2d0*PAR(7)+PAR(8)+PAR(9)
     .       +4d0*PAR(15)+2d0*PAR(16)+2d0*PAR(17))/12.d0)

c      5) Higgs parameters
       l=PAR(1)
       k=PAR(2)
       mu=PAR(4)
       M1=PAR(20)
       IF(M1.ge.0d0)THEN
        M1=Max(1d0,M1)
       ELSE
        M1=-Max(1d0,dabs(M1))
       ENDIF
       M2=PAR(21)
       IF(M2.ge.0d0)THEN
        M2=Max(1d0,M2)
       ELSE
        M2=-Max(1d0,dabs(M2))
       ENDIF
       M3=PAR(22)
       IF(M3.ge.0d0)THEN
        M3=Max(1d0,M3)
       ELSE
        M3=-Max(1d0,dabs(M3))
       ENDIF
       msi=Max(1d0,
     . dsqrt(4d0*(k/l*mu)**2+muP**2+4d0*k/l*muP*mu*dcos(phi02-phiP)))

      IF(MAFLAG.LT.0)THEN
       Alcos1=PAR(5)
       Akcos2=PAR(6)
       MA2=PAR(23)**2
      ELSE
       IF(MOD(MAFLAG,3).EQ.0)THEN
        Alcos1=PAR(5)
        MA2=MAX(((Alcos1+k*mu/l*dcos(Phi0))*mu+M3H*dcos(phi3)
     .           +mu*mup*dcos(phi01-phip)+l*xif*dcos(phi01-phiF))
     .          *(tanb+1d0/tanb),1d0)
        PAR(23)=DSQRT(MAX(MA2,1d0))
       ELSEIF(MOD(MAFLAG,3).EQ.1)THEN
        MA2=PAR(23)**2
        Alcos1=(MA2*TANB/(1d0+TANB**2)-M3H*dcos(phi3)
     .          -l*xif*dcos(phi01-phiF))/mu-MUP*dcos(phi01-phip)
     .          -k*mu/l*dcos(Phi0)
        PAR(5)=Alcos1
       ELSE
        Alcos1=PAR(5)
        MA2=PAR(23)**2
        aux2=XIF
        aux1=(MA2*TANB/(1d0+TANB**2)-(Alcos1+k*mu/l*dcos(Phi0))*mu
     .  -M3H*dcos(phi3)-mu*mup*dcos(phi01-phip))/(l*dcos(phi01))
     .  -datan(phi01)*XIF
        IF(aux1.eq.0d0)aux1=1d-20
        phiF=datan(aux2/aux1)
        XIF=aux1/DABS(aux1)*DSQRT(aux1**2+aux2**2)
       ENDIF
       IF(MAFLAG/3.EQ.0)THEN 
        Akcos2=PAR(6)
        MP2=-k/l*mu*(3.d0*Akcos2+mup*dcos(phi02-phip))
     . +l**2*vu*vd/mu*(Alcos1+muP*dcos(phiP-phi01)
     .                              +4.d0*k/l*mu*dcos(phi01-phi02))
     . -l/mu*(xiS*dcos(phiS)+XIF*mup*dcos(phiP-phiF))
     . -2d0*MSP*dcos(phiSP)-4d0*k*XIF*dcos(phi02-phiF)
        PAR(24)=DSQRT(MAX(MP2,1d0))
       ELSEIF(MAFLAG/3.EQ.1)THEN
        MP2=PAR(24)**2
        IF(K.EQ.0d0)THEN
         Akcos2=0d0
         PAR(6)=0d0
        ELSE
         Akcos2=(l**2*(Alcos1+4d0*k*mu/l*dcos(Phi0)
     .                        +mup*dcos(phi01-phip))*vu*vd/mu
     .  -l*(XIF*mup*dcos(phip-phiF)+xiS*dcos(phiS))/mu
     .  -mup*k*mu/l*dcos(phi02-phip)-4d0*k*XIF*dcos(phi02-phiF)
     .  -2d0*MSP*dcos(phiSP)-MP2)/(3d0*k*mu/l*dcos(Phi0))
         PAR(6)=Akcos2
        ENDIF
       ELSE
        MP2=PAR(24)**2
        Akcos2=PAR(6)
        aux2=XIS
        aux1=(l*(Alcos1+4d0*k*mu/l*dcos(Phi0)+mup*dcos(phi01-phip))
     .      *vu*vd-(k*mu/l*(3d0*Akcos2+mup*dcos(phi02-phip))
     .       +4d0*k*XIF*dcos(phi02-phiF)+2d0*MSP*dcos(phiSP)+MP2)*mu/l
     .      -XIF*muP*dcos(phip-phiF))
        IF(k.ne.0d0)THEN
        IF(aux1.eq.0d0)aux1=1d-20
         phiS=datan(aux2/aux1)
         XIS=aux1/DABS(aux1)*DSQRT(aux1**2+aux2**2)
        ELSE
         XIS=aux1
        ENDIF
       ENDIF
      ENDIF

       IAl=-k/l*mu*dsin(phi0)-(M3H*dsin(phi3)+l*XIF*dsin(phi01-phiF)
     .                              +mu*mup*dsin(phi01-phip))/mu
       IF(k.ne.0d0)THEN
         IAk=l**2/mu**2/k*(l*vu*vd*(IAl-muP*dsin(phi01-phiP)
     .                                     -2d0*k/l*mu*dsin(phi0))
     .      -xiS*dsin(phiS)-mu/l*MSP*dsin(phiSP)
     .      -muP*k*(mu/l)**2*dsin(phi02-phiP)
     .      -XIF*(muP*dsin(phiP-phiF)+2d0*k/l*mu*dsin(phi02-phiF)))
       ELSE
         IAk=0d0
         aux1=XIS
         aux2=l*vu*vd*(IAl-muP*dsin(phi01-phiP))-mu/l*MSP*dsin(phiSP)
     .          -XIF*muP*dsin(phiP-phiF)
         IF(aux1.eq.0d0)aux1=1d-20
         phiS=datan(aux2/aux1)
         XIS=aux1/DABS(aux1)*dsqrt(aux1**2+aux2**2)
       ENDIF
         IXIS=XIS*dsin(phiS)

c      6) Running gauge couplings

      ALSQ=ALSMT/(1d0+ALSMT/(4d0*PI)*(7d0*DLOG(MAX(QSTSB,MT**2)/MT**2)
     .         -2d0*DLOG(MAX(QSTSB,PAR(22)**2)/MAX(PAR(22)**2,MT**2))))

*   g_2**2 including the Higgs and sparticle thresholds:
      g2q=g2/(1d0+g2/16d0/Pi**2*(DLOG(QSTSB/MZ**2)*19d0/6d0
     .    -DLOG(QSTSB/MIN(QSTSB,MAX(MA2,MZ**2)))/6d0
     .    -DLOG(QSTSB/MIN(QSTSB,MAX(MU**2,MZ**2)))*2d0/3d0
     .    -DLOG(QSTSB/MIN(QSTSB,MAX(PAR(7),MZ**2)))/2d0
     .    -DLOG(QSTSB/MIN(QSTSB,MAX(PAR(10),MZ**2)))/6d0
     .    -DLOG(QSTSB/MIN(QSTSB,MAX(PAR(15),MZ**2)))
     .    -DLOG(QSTSB/MIN(QSTSB,MAX(PAR(18),MZ**2)))/3d0
     .    -DLOG(QSTSB/MIN(QSTSB,MAX(M2**2,MZ**2)))*4d0/3d0))

*   g_1**2 including the top, Higgs and sparticle thresholds:
      g1q=g1/(1d0-g1/16d0/Pi**2*(DLOG(QSTSB/MZ**2)*53d0/9d0
     .    +DLOG(QSTSB/MT**2)*17d0/18d0
     .    +DLOG(QSTSB/MIN(QSTSB,MAX(MA2,MZ**2)))/6d0
     .    +DLOG(QSTSB/MIN(QSTSB,MAX(MU**2,MZ**2)))*2d0/3d0
     .    +DLOG(QSTSB/MIN(QSTSB,MAX(PAR(7),MZ**2)))/18d0
     .    +DLOG(QSTSB/MIN(QSTSB,MAX(PAR(8),MZ**2)))*4d0/9d0
     .    +DLOG(QSTSB/MIN(QSTSB,MAX(PAR(9),MZ**2)))/9d0
     .    +DLOG(QSTSB/MIN(QSTSB,MAX(PAR(10),MZ**2)))/6d0
     .    +DLOG(QSTSB/MIN(QSTSB,MAX(PAR(11),MZ**2)))/3d0
     .    +DLOG(QSTSB/MIN(QSTSB,MAX(PAR(15),MZ**2)))/9d0
     .    +DLOG(QSTSB/MIN(QSTSB,MAX(PAR(16),MZ**2)))*8d0/9d0
     .    +DLOG(QSTSB/MIN(QSTSB,MAX(PAR(17),MZ**2)))*2d0/9d0
     .    +DLOG(QSTSB/MIN(QSTSB,MAX(PAR(18),MZ**2)))/3d0
     .    +DLOG(QSTSB/MIN(QSTSB,MAX(PAR(19),MZ**2)))*2d0/3d0))

      gq=(g1q+g2q)/2d0

c      7) Running Yukawa couplings

      ALSMA=ALSMT/(1d0+ALSMT/(4d0*PI)
     . *(7d0*DLOG(MIN(MAX(MA2,MT**2),QSTSB)/MT**2)
     .   -2d0*DLOG(MAX(MA2,PAR(22)**2)/MAX(PAR(22)**2,MT**2))))
      DLA=(ALSMA/ALSMT)**(1d0/7d0)
      DLQA=(ALSQ/ALSMA)**(1d0/7d0)
      F1=DABS(1d0-9d0*sinb**2*Yt**2*(1d0-DLA)/(8d0*PI*ALSMT))
      HTMA=YT*DLA**4/DSQRT(F1)
      F2=DABS(1d0-9d0*HTMA**2*(1d0-DLQA)/(8d0*PI*ALSMA))

*   Top/Bottom Yukawas at QSTSB
*   (Note: RUNMB(Q) includes QCD corrections only)
*   including electroweak contributions, Logs ht^2*LQT resummed:
*   + Conversion to DR_bar:
      Ytq=Yt/DSQRT(F1*F2)*(1d0+7d0/(4d0*PI)*ALSMT
     .                    *DLOG(MAX(QSTSB,MT**2)/MT**2))**(-4d0/7d0)
     .   *(1d0+((-17d0/6d0*g1q-9d0/2d0*g2q+Yb**2)
     .                    *DLOG(MAX(QSTSB,MT**2)/MT**2)
     .   +((3d0*cosb**2-1d0)*Yb**2+2d0*Ytau**2*cosb**2)
     .                    *DLOG(MIN(MAX(MA2,MT**2),QSTSB)/MT**2)
     .   -2d0*l**2*DLOG(MIN(MAX(mu**2,msi**2,MZ**2),QSTSB)/QSTSB)
     .   -G1Q*DLOG(MIN(MAX(M1**2,mu**2,MZ**2),QSTSB)/QSTSB)
     .   -3d0*G2Q*DLOG(MIN(MAX(M2**2,mu**2,MZ**2),QSTSB)/QSTSB))
     .                                     /64d0/Pi**2)
     .   *(1d0-ALSQ/(3d0*PI)+g2q*3d0/128d0/Pi**2)

      Ybq=RUNMB(DSQRT(QSTSB))/vd*F1**(-1d0/6d0)
     .   *(1d0-3d0*HTMA**2*(1d0-DLQA)/(8d0*PI*ALSMA))**(-1d0/6d0)
     .   *(1d0+((-5d0/6d0*g1q-9d0/2d0*g2q+9d0*Yb**2+2d0*Ytau**2)
     .                     *DLOG(MAX(QSTSB,MT**2)/MT**2)
     .   +(-9d0*Yb**2-2d0*Ytau**2)*sinb**2
     .                     *DLOG(MIN(MAX(MA2,MT**2),QSTSB)/MT**2)
     .   -2d0*L**2*DLOG(MIN(MAX(mu**2,msi**2,MZ**2),QSTSB)/QSTSB)
     .   -G1Q*DLOG(MIN(MAX(M1**2,mu**2,MZ**2),QSTSB)/QSTSB)
     .   -3d0*G2Q*DLOG(MIN(MAX(M2**2,mu**2,MZ**2),QSTSB)/QSTSB))
     .                                     /64d0/Pi**2)
     .   *(1d0-ALSQ/(3d0*PI)+g2q*3d0/128d0/Pi**2)

c      8) Higgs parameters at the SUSY scale

      LQ=l*(1d0+(-G1Q-3d0*G2Q+4d0*L**2+2d0*K**2
     .           +3d0*Ytq**2+3d0*Ybq**2+Ytau**2)/32d0/Pi**2
     .                                      *DLOG(QSTSB/Q2))
      KQ=K*(1d0+3d0/16d0/Pi**2*(L**2+K**2)*DLOG(QSTSB/Q2))
      MUQ=MU*(1d0+(3d0*Ytq**2+3d0*Ybq**2+Ytau**2+2d0*l**2
     .             -G1Q-3d0*G2Q)/32d0/Pi**2*DLOG(QSTSB/Q2))
      NUQ=MUQ*KQ/LQ
      MUPQ=MUP*(1d0+(L**2+K**2)/8d0/Pi**2*DLOG(QSTSB/Q2))
      XIFQ=XIF*(1d0+(L**2+K**2)/16d0/Pi**2*DLOG(QSTSB/Q2))

c      XISQ=XIS+(L**2*(XIS+2d0*XIF
c     .  *(Alcos1*dcos(phiF-phi01-phiS)-IAl*dsin(phiF-phi01-phiS)))
c     . +K**2*(XIS+2d0*XIF
c     .  *(Akcos2*dcos(phiF-phi02-phiS)-IAk*dsin(phiF-phi02-phiS)))
c     . +2d0*L*M3H*(Alcos1*dcos(phi3+phiS)+IAl*dsin(phi3+phiS)
c     .                        +MUP*dcos(phi3+phiP-phi01-phiS))
c     . +K*MSP*(Akcos2*dcos(phiSP+phiS)+IAk*dsin(phiSP+phiS)
c     .                +MUP*dcos(phip+phiSP-phi02-phiS)))
c     .                      /16d0/Pi**2*DLOG(QSTSB/Q2)
      MSPQ=MSP+(2d0*L**2*(MSP+2d0*MUP
     .      *(Alcos1*dcos(phip-phi01-phiSP)-IAl*dsin(phip-phi01-phiSP)))
     . +4d0*K**2*(MSP+2d0*MUP
     .      *(Akcos2*dcos(phip-phi02-phiSP)-IAk*dsin(phip-phi02-phiSP)))
     . +4d0*L*K*M3H*dcos(phi3-phiSP-phi0))/16d0/Pi**2*DLOG(QSTSB/Q2)
      M3HQ=M3H+((3d0*Ytq**2+3d0*Ybq**2+Ytau**2+6d0*l**2-G1Q-3d0*G2Q)
     .    *M3H+2d0*L*K*MSP*dcos(phiSP-phi3+phi0))/32d0/Pi**2
     .                                *DLOG(QSTSB/Q2)
      phiSq=phiS+(2d0*XIF*L**2
     .  *(Alcos1*dsin(phiF-phi01-phiS)+IAl*dcos(phiF-phi01-phiS))
     . +2d0*XIF*K**2
     .  *(Akcos2*dsin(phiF-phi02-phiS)+IAk*dcos(phiF-phi02-phiS))
     . +2d0*L*M3H*(IAl*dcos(phi3+phiS)-Alcos1*dsin(phi3+phiS)
     .                        +MUP*dsin(phi3+phiP-phi01-phiS))
     . +K*MSP*(IAk*dcos(phiSP-phiS)-Akcos2*dsin(phiSP-phiS)
     .         +MUP*dsin(phip+phiSP-phi02-phiS)))
     .       /16d0/Pi**2/Max(xiS,1d-4)*DLOG(QSTSB/Q2)
      phiSPq=phiSP+(2d0*L**2*(IAl*dcos(phip-phi01-phiSP)
     .                         +Alcos1*dsin(phip-phi01-phiSP))
     . +4d0*K**2*2d0*MUP*(IAk*dcos(phip-phi02-phiSP)
     .                     +Akcos2*dsin(phip-phi02-phiSP))
     . +4d0*L*K*M3H*dsin(phi3-phiSP-phi0))
     .    /Max(MSP,1d-4)/16d0/Pi**2*DLOG(QSTSB/Q2)
      phi3q=phi3+l*k*MSP/Max(M3H,1d-4)*dsin(phiSP-phi3+phi0)
     .                        /8d0/Pi**2*DLOG(QSTSB/Q2)

*  In runpar.f, the definition of Al / MA / Alq is different when input=MA:
*             Alq is defined from MA and Al is deduced. 
*             Here, Al is derived from MA and Alq is deduced.
*  Similar discrepancies in PAR(5), PAR(24), etc.

c       Alcos1=Alcos1+(4d0*l**2*Alcos1
c     .    +2d0*k**2*(Akcos2*dcos(phi0)-IAk*dsin(phi0))
c     .    +3d0*Ytq**2*PAR(12)*dcos(phiAt)
c     .    +3d0*Ybq**2*PAR(13)*dcos(phiAb)
c     .    +Ytau**2*PAR(14)*dcos(phiAtau)
c     .    +g1q*M1*dcos(phiM1)
c     .    +3d0*g2q*M2*dcos(phiM2))*DLOG(QSTSB/Q2)/16d0/Pi**2
c       Akcos2=Akcos2+6d0*(l**2*(PAR(5)*dcos(phi0)+IAl*dsin(phi0))
c     .     +k**2*Akcos2)*DLOG(QSTSB/Q2)/16d0/Pi**2


      IF(MAFLAG.LT.0)THEN
       Alcos1=Alcos1+(4d0*l**2*Alcos1
     .    +2d0*k**2*(Akcos2*dcos(phi0)-IAk*dsin(phi0))
     .    +3d0*Ytq**2*PAR(12)*dcos(phiAt)
     .    +3d0*Ybq**2*PAR(13)*dcos(phiAb)
     .    +Ytau**2*PAR(14)*dcos(phiAtau)
     .    +g1q*M1*dcos(phiM1)
     .    +3d0*g2q*M2*dcos(phiM2))*DLOG(QSTSB/Q2)/16d0/Pi**2
       Akcos2=Akcos2+6d0*(l**2*(PAR(5)*dcos(phi0)+IAl*dsin(phi0))
     .     +k**2*Akcos2)*DLOG(QSTSB/Q2)/16d0/Pi**2
       XIFQ=XIF*(1d0+(L**2+K**2)/16d0/Pi**2*DLOG(QSTSB/Q2))
       XISQ=XIS+(L**2*(XIS+2d0*XIF
     .  *(Alcos1*dcos(phiF-phi01-phiS)-IAl*dsin(phiF-phi01-phiS)))
     . +K**2*(XIS+2d0*XIF
     .  *(Akcos2*dcos(phiF-phi02-phiS)-IAk*dsin(phiF-phi02-phiS)))
     . +2d0*L*M3H*(Alcos1*dcos(phi3+phiS)+IAl*dsin(phi3+phiS)
     .                        +MUP*dcos(phi3+phiP-phi01-phiS))
     . +K*MSP*(Akcos2*dcos(phiSP+phiS)+IAk*dsin(phiSP+phiS)
     .                +MUP*dcos(phip+phiSP-phi02-phiS)))
     .                      /16d0/Pi**2*DLOG(QSTSB/Q2)
      ELSE
       IF(MOD(MAFLAG,3).EQ.0)THEN
        Alcos1=Alcos1+(4d0*l**2*Alcos1
     .    +2d0*k**2*(Akcos2*dcos(phi0)-IAk*dsin(phi0))
     .    +3d0*Ytq**2*PAR(12)*dcos(phiAt)
     .    +3d0*Ybq**2*PAR(13)*dcos(phiAb)
     .    +Ytau**2*PAR(14)*dcos(phiAtau)
     .    +g1q*M1*dcos(phiM1)
     .    +3d0*g2q*M2*dcos(phiM2))*DLOG(QSTSB/Q2)/16d0/Pi**2
        XIFQ=XIF*(1d0+(L**2+K**2)/16d0/Pi**2*DLOG(QSTSB/Q2))
        MA2=((Alcos1+NUQ*dcos(Phi0)+MUPQ*dcos(PhiP-Phi01))*MUQ
     . +M3HQ*dcos(phi3)+LQ*XIFQ*dcos(PhIF-Phi01))*(TANB+1d0/TANB)
        PAR(23)=DSQRT(MAX(MA2,1d0))
       ELSEIF(MOD(MAFLAG,3).EQ.1)THEN
        XIFQ=XIF*(1d0+(L**2+K**2)/16d0/Pi**2*DLOG(QSTSB/Q2))
        Alcos1=(MA2*TANB/(1d0+TANB**2)-M3HQ*dcos(phi3)
     . -LQ*XIFQ*dcos(PhIF-Phi01))/MUQ-MUPQ*dcos(PhiP-Phi01)
     . -NUQ*dcos(Phi0)
        PAR(5)=Alcos1-(4d0*l**2*Alcos1
     .    +2d0*k**2*(Akcos2*dcos(phi0)-IAk*dsin(phi0))
     .    +3d0*Ytq**2*PAR(12)*dcos(phiAt)
     .    +3d0*Ybq**2*PAR(13)*dcos(phiAb)
     .    +Ytau**2*PAR(14)*dcos(phiAtau)
     .    +g1q*M1*dcos(phiM1)
     .    +3d0*g2q*M2*dcos(phiM2))*DLOG(QSTSB/Q2)/16d0/Pi**2
       ELSE
        Alcos1=Alcos1+(4d0*l**2*Alcos1
     .    +2d0*k**2*(Akcos2*dcos(phi0)-IAk*dsin(phi0))
     .    +3d0*Ytq**2*PAR(12)*dcos(phiAt)
     .    +3d0*Ybq**2*PAR(13)*dcos(phiAb)
     .    +Ytau**2*PAR(14)*dcos(phiAtau)
     .    +g1q*M1*dcos(phiM1)
     .    +3d0*g2q*M2*dcos(phiM2))*DLOG(QSTSB/Q2)/16d0/Pi**2
        XIFQ=((MA2*TANB/(1d0+TANB**2)
     .   -(Alcos1+NUQ*dcos(Phi0))*MUQ-M3HQ*dcos(Phi3)
     .   -MUQ*MUPQ*dcos(PhiP-Phi01))/LQ)/Max(dcos(PhIF-Phi01),1d-4)
        XIF=XIFQ*(1d0-(L**2+K**2)/16d0/Pi**2*DLOG(QSTSB/Q2))
       ENDIF
       IF(MAFLAG/3.EQ.0)THEN
        IF(k.ne.0d0)then
         Akcos2=Akcos2+6d0*(l**2*(PAR(5)*dcos(phi0)+IAl*dsin(phi0))
     .         +k**2*Akcos2)*DLOG(QSTSB/Q2)/16d0/Pi**2
        ELSE
         Akcos2=0d0
        ENDIF
        XISQ=XIS+(L**2*(XIS+2d0*XIF
     .  *(Alcos1*dcos(phiF-phi01-phiS)-IAl*dsin(phiF-phi01-phiS)))
     . +K**2*(XIS+2d0*XIF
     .  *(Akcos2*dcos(phiF-phi02-phiS)-IAk*dsin(phiF-phi02-phiS)))
     . +2d0*L*M3H*(Alcos1*dcos(phi3+phiS)+IAl*dsin(phi3+phiS)
     .                        +MUP*dcos(phi3+phiP-phi01-phiS))
     . +K*MSP*(Akcos2*dcos(phiSP+phiS)+IAk*dsin(phiSP+phiS)
     .                +MUP*dcos(phip+phiSP-phi02-phiS)))
     .                      /16d0/Pi**2*DLOG(QSTSB/Q2)
        MP2=LQ**2*(Alcos1+4d0*NUQ+MUPQ)*vu*vd/MUQ-3d0*Akcos2*NUQ
     .    -LQ*(XIFQ*MUPQ+XISQ)/MUQ-MUPQ*NUQ-4d0*KQ*XIFQ-2d0*MSPQ
        PAR(24)=DSQRT(MAX(MP2,1d0))
       ELSEIF(MAFLAG/3.EQ.1)THEN
        XISQ=XIS+(L**2*(XIS+2d0*XIF
     .  *(Alcos1*dcos(phiF-phi01-phiS)-IAl*dsin(phiF-phi01-phiS)))
     . +K**2*(XIS+2d0*XIF
     .  *(Akcos2*dcos(phiF-phi02-phiS)-IAk*dsin(phiF-phi02-phiS)))
     . +2d0*L*M3H*(Alcos1*dcos(phi3+phiS)+IAl*dsin(phi3+phiS)
     .                        +MUP*dcos(phi3+phiP-phi01-phiS))
     . +K*MSP*(Akcos2*dcos(phiSP+phiS)+IAk*dsin(phiSP+phiS)
     .                +MUP*dcos(phip+phiSP-phi02-phiS)))
     .                      /16d0/Pi**2*DLOG(QSTSB/Q2)
        IF(K.EQ.0d0)THEN
         Akcos2=0d0
         PAR(6)=0d0
        ELSE
         Akcos2=(LQ**2*(Alcos1+4d0*NUQ+MUPQ)*vu*vd/MUQ
     .      -LQ*(XIFQ*MUPQ+XISQ)/MUQ-MUPQ*NUQ
     .      -4d0*KQ*XIFQ-2d0*MSPQ-MP2)/(3d0*NUQ)
         PAR(6)=Akcos2-6d0*(l**2*(PAR(5)*dcos(phi0)+IAl*dsin(phi0))
     .     +k**2*Akcos2)*DLOG(QSTSB/Q2)/16d0/Pi**2
        ENDIF
       ELSE
        Akcos2=Akcos2+6d0*(l**2*(PAR(5)*dcos(phi0)+IAl*dsin(phi0))
     .     +k**2*Akcos2)*DLOG(QSTSB/Q2)/16d0/Pi**2
        XISQ=(LQ*(Alcos1+4d0*NUQ+MUPQ)*vu*vd-(3d0*Akcos2*NUQ+MUPQ*NUQ
     .    +4d0*KQ*XIFQ+2d0*MSPQ+MP2)*MUQ/LQ-XIFQ*MUPQ)
     .     /Max(dcos(phiSq),1d-4)
        XIS=XISQ-(L**2*(XIS+2d0*XIF
     .  *(Alcos1*dcos(phiF-phi01-phiS)-IAl*dsin(phiF-phi01-phiS)))
     . +K**2*(XIS+2d0*XIF
     .  *(Akcos2*dcos(phiF-phi02-phiS)-IAk*dsin(phiF-phi02-phiS)))
     . +2d0*L*M3H*(Alcos1*dcos(phi3+phiS)+IAl*dsin(phi3+phiS)
     .                        +MUP*dcos(phi3+phiP-phi01-phiS))
     . +K*MSP*(Akcos2*dcos(phiSP+phiS)+IAk*dsin(phiSP+phiS)
     .                +MUP*dcos(phip+phiSP-phi02-phiS)))
     .                      /16d0/Pi**2*DLOG(QSTSB/Q2)
       ENDIF
      ENDIF

c      9) Gaugino parameters

      mhs2=Max(l**2*vu*vd/mu*(PAR(5)+mup*dcos(phip-phi01))
     . +k/l*mu*(PAR(6)+3d0*mup*dcos(phiP-phi02)+4.d0*k/l*mu)
     . -l/mu*(xiS*dcos(phiS)+XIF*mup*dcos(phiP-phiF)),1d0)

      mas2=Max(-k/l*mu*(3.d0*PAR(6)+mup*dcos(phi02-phip))
     . +l**2*vu*vd/mu*(PAR(5)+muP*dcos(phiP-phi01)
     .                              +4.d0*k/l*mu*dcos(phi01-phi02))
     . -l/mu*(xiS*dcos(phiS)+XIF*mup*dcos(phiP-phiF))
     . -2d0*MSP*dcos(phiSP)-4d0*k*XIF*dcos(phi02-phiF),1d0)

      mur=mu*(1d0-1d0/(64d0*Pi**2)*(
     .       12d0*(Ytq**2+Ybq**2)*NMB1(mu**2,0.d0,QSUSY**2,Q2)
     . +3d0*g2*(NMB1(mu**2,M2**2,PAR(23)**2,Q2)
     .    +NMB1(MU**2,M2**2,MZ**2,Q2)+2d0*NMB1(MU**2,MU**2,MZ**2,Q2)
     .    -4d0*NMB0(MU**2,MU**2,MZ**2,Q2)
     .     +2d0*sinb*cosb*mu/M2*(NMB0(MU**2,M2**2,PAR(23)**2,Q2)
     .     -NMB0(MU**2,M2**2,MZ**2,Q2))*dcos(PhiM2+Phi01))
     . +g1*(NMB1(mu**2,M1**2,PAR(23)**2,Q2)
     .    +NMB1(MU**2,M1**2,MZ**2,Q2)+2d0*NMB1(MU**2,MU**2,MZ**2,Q2)
     .    -4d0*NMB0(MU**2,MU**2,MZ**2,Q2)
     .   +2d0*sinb*cosb*M1/mu*(NMB0(MU**2,M1**2,PAR(23)**2,Q2)
     .     -NMB0(MU**2,M1**2,MZ**2,Q2))*dcos(PhiM1+Phi01))
     . +2.d0*l**2*(NMB1(mu**2,msi**2,PAR(23)**2,Q2)
     .    +NMB1(MU**2,msi**2,MZ**2,Q2)
     .    -2d0*sinb*cosb*2d0*k/l*(NMB0(MU**2,msi**2,PAR(23)**2,Q2)
     .     -NMB0(MU**2,msi**2,MZ**2,Q2))*dcos(Phi01-Phi02))))

      M1r=M1*(1d0-g1/(16d0*Pi**2)*(11d0*NMB1(M1**2,0d0,QSQ**2,Q2)
     .     +9d0* NMB1(M1**2,0d0,QSL**2,Q2)
     .     +2d0*sinb*cosb*mu/M1*(NMB0(M1**2,MU**2,PAR(23)**2,Q2)
     .     -NMB0(M1**2,MU**2,MZ**2,Q2))*dcos(PhiM1+Phi01)
     .     +NMB1(M1**2,MU**2,PAR(23)**2,Q2)+NMB1(M1**2,MU**2,MZ**2,Q2)))

      M2r=M2*(1d0-g2/(16d0*Pi**2)*(9d0*NMB1(M2**2,0d0,QSQ**2,Q2)
     .     +3d0* NMB1(M2**2,0d0,QSL**2,Q2)
     .     +2d0*sinb*cosb*mu/M2*(NMB0(M2**2,MU**2,PAR(23)**2,Q2)
     .     -NMB0(M2**2,MU**2,MZ**2,Q2))*dcos(PhiM2+Phi01)
     .     +NMB1(M2**2,MU**2,PAR(23)**2,Q2)+NMB1(M2**2,MU**2,MZ**2,Q2)
     .  -8d0*NMB0(M2**2,M2**2,MW**2,Q2)+4d0*NMB1(M2**2,M2**2,MW**2,Q2)))

      msi=msi*(1d0-2d0*l**2/(16d0*Pi**2)*(
     .   NMB1(msi**2,MU**2,PAR(23)**2,Q2)+NMB1(msi**2,MU**2,MZ**2,Q2)))

      msi=msi-dsqrt(4d0*(k/l*mu)**2+muP**2
     .   +4d0*k/l*muP*mu*dcos(phi02-phiP))*2d0*k**2/(16d0*Pi**2)*(
     .   -(NMB0(Max(1d0,msi**2),msi**2,mhs2,Q2)
     .     -NMB0(Max(1d0,msi**2),msi**2,mas2,Q2))
     .   +NMB1(Max(1d0,msi**2),msi**2,mhs2,Q2)
     .   +NMB1(Max(1d0,msi**2),msi**2,mas2,Q2))


      mupsi=mup*msi/Max(dsqrt(4d0*(k/l*mu)**2+muP**2
     .                         +4d0*k/l*muP*mu*dcos(phi02-phiP)),1d0)

      ks2si=2d0*k/l*mu*msi/Max(dsqrt(4d0*(k/l*mu)**2+muP**2
     .                         +4d0*k/l*muP*mu*dcos(phi02-phiP)),1d0)

      ZHU=1.d0+1.d0/16.d0/Pi**2*(3.d0*Ytq**2
cUE     c                         *NMB0(muH2,mtopq**2,mtopq**2,QSTSB)
     .                         *NMB0(muH2,mt**2,mt**2,QSTSB)
     .    -(g1q+g2q)/2.d0*NMB0(muH2,MZ**2,MZ**2,QSTSB)*sinb**2 !PAR(23)**2
     .    -g2q*NMB0(muH2,MW**2,MW**2,QSTSB)*sinb**2 !PAR(23)**2
     .    +g1q/2.d0*NMB0(muH2,mur**2,M1r**2,QSTSB)
     .    +3.d0*g2q/2.d0*NMB0(muH2,mur**2,M2r**2,QSTSB)
     .    +l**2*NMB0(muH2,mur**2,msi**2,QSTSB))
      ZHD=1.d0+1.d0/16.d0/Pi**2*(3.d0*Ybq**2
cUE     c                         *NMB0(muH2,mbotq**2,mbotq**2,QSTSB)
     .                         *NMB0(muH2,mb**2,mb**2,QSTSB)
     .    +Ytau**2*NMB0(muH2,mtau**2,mtau**2,QSTSB)
c     c    -(g1q+g2q)/2.d0*NMB0(muH2,MZ**2,MZ**2,QSTSB) !PAR(23)**2
c     c    -g2q*NMB0(muH2,MW**2,MW**2,QSTSB) !PAR(23)**2
     .    +g1q/2.d0*NMB0(muH2,mur**2,M1r**2,QSTSB)
     .    +3.d0*g2q/2.d0*NMB0(muH2,mur**2,M2r**2,QSTSB)
     .    +l**2*NMB0(muH2,mur**2,msi**2,QSTSB))
      ZS=1.d0+1.d0/8.d0/Pi**2
     .    *(l**2*NMB0(muH2,mur**2,mur**2,QSTSB)
     .     +k**2*NMB0(muH2,msi**2,msi**2,QSTSB))

      vuq=vu/dsqrt(ZHU)
      vdq=vd/dsqrt(ZHD)
      tanbq=vuq/vdq
      MTOPQ=YTQ*vuQ
      MBOTQ=YBQ*vdQ
      muq=muq/dsqrt(ZS)

c      10) Squark and slepton
*   In order to run the squark masses from Q2 to QSTSB
*   including all thresholds, tree level expressions for the soft
*   Higgs masses and various logarithms are needed:

      AT=PAR(12)
      AB=PAR(13)

      MHuS=(MUQ*(Alcos1+MUPQ*dcos(phi01-phiP)+kq/lq*muq*dcos(Phi0))
     .         +M3HQ+LQ*XIFQ*dcos(phi01-phiF))/TANBQ
     .      -LQ**2*VDQ**2-MUQ**2+GQ/2d0*(VDQ**2-VUQ**2)
      MHdS=(MUQ*(Alcos1+MUPQ*dcos(phi01-phiP)+kq/lq*muq*dcos(Phi0))
     .         +M3HQ+LQ*XIFQ*dcos(phi01-phiF))*TANBQ
     .      -LQ**2*VUQ**2-MUQ**2+GQ/2d0*(VUQ**2-VDQ**2)
      MSS=LQ**2*vuq*vdq/muq
     .        *(Alcos1+2d0*kq/lq*muq*dcos(phi0)+mupq*dcos(Phi01-PhiP))
     . -kq/lq*muq*(Akcos2+2d0*kq/lq*muq)-lq**2*(vuq**2+vdq**2)
     .  -XIFQ*(2d0*KQ*dcos(phi02-PhiF)+LQ*MUPQ*dcos(PhiP-PhiF)/MUQ)
     .  -MUPQ**2-3d0*MUPQ*KQ/LQ*MUQ*dcos(Phi02-PhiP)-MSPQ*dcos(PhiSPq)
     .  -LQ*XISQ/MUQ*dcos(PhiSq)

      MQ3=PAR(7)
      MU3=PAR(8)
      MD3=PAR(9)
      ML3=PAR(10)
      ME3=PAR(11)
      MQ=PAR(15)
      MU1=PAR(16)
      MD=PAR(17)
      ML=PAR(18)
      ME=PAR(19)

*   Running stop/sbottom masses squared and mixings:
*   On input, they are defined at the scale Q2.
*   First, they have to be evaluated at the scale QSTSB
*   (as used in the calculation of the pole masses)

c*********Old soft sfermion running
c*   Taking all possible thresholds into account
c*   (with the exception of the gluino threshold that
c*   is taken into account in the calculation of the pole masses).
c      ANOMQSTSB=G1Q*((MHuS-MHdS)*DLOG(Q2/QSTSB)
c     .              +MQ3*DLOG(Q2/MAX(MQ3,QSTSB))
c     .              -2d0*MU3*DLOG(Q2/MAX(MU3,QSTSB))
c     .              +MD3*DLOG(Q2/MAX(MD3,QSTSB))
c     .              +ME3*DLOG(Q2/MAX(ME3,QSTSB))
c     .              -ML3*DLOG(Q2/MAX(ML3,QSTSB))
c     .       +2d0*(MQ*DLOG(Q2/MAX(MQ,QSTSB))
c     .            -2d0*MU1*DLOG(Q2/MAX(MU1,QSTSB))
c     .            +MD*DLOG(Q2/MAX(MD,QSTSB))
c     .            +ME*DLOG(Q2/MAX(ME,QSTSB))
c     .            -ML*DLOG(Q2/MAX(ML,QSTSB))))

c        RTOP=YTQ**2*(MHuS*DLOG(Q2/QSTSB)
c     .  +MQ3*DLOG(Q2/MAX(MQ3,QSTSB))+MU3*DLOG(Q2/MAX(MU3,QSTSB))
c     .  +AT**2*DLOG(Q2/QSTSB))

c        RBOT=YBQ**2*(MHdS*DLOG(Q2/QSTSB)
c     .  +MQ3*DLOG(Q2/MAX(MQ3,QSTSB))+MD3*DLOG(Q2/MAX(MD3,QSTSB))
c     .  +AB**2*DLOG(Q2/QSTSB))

c        MSQ3=MQ3-(RTOP+RBOT
c     .       -G1Q*M1**2*DLOG(Q2/MAX(M1**2,QSTSB))/9d0
c     .       -3d0*G2Q*M2**2*DLOG(Q2/MAX(M2**2,QSTSB))
c     .       +ANOMQSTSB/3d0)/16d0/Pi**2

c        MSU3=MU3-(2d0*RTOP
c     .       -16d0/9d0*G1Q*M1**2*DLOG(Q2/MAX(M1**2,QSTSB))
c     .       -2d0/3d0*ANOMQSTSB)/16d0/Pi**2
c        print*,RTOP,RBOT,G1Q*M1**2*DLOG(Q2/MAX(M1**2,QSTSB)),
c     .       3d0*G2Q*M2**2*DLOG(Q2/MAX(M2**2,QSTSB)),ANOMQSTSB
c        MSD3=MD3-(2d0*RBOT
c     .       -4d0/9d0*G1Q*M1**2*DLOG(Q2/MAX(M1**2,QSTSB))
c     .       +1d0/3d0*ANOMQSTSB)/16d0/Pi**2
c
c*  Integrate the trilinears from Q2 to QSTSB:
c      aux1=AT*dcos(PhiAT)+
c     .   ((6d0*YTQ**2*AT*dcos(PhiAt)+YBQ**2*AB*dcos(PhiAB)
c     .      +LQ**2*(Alcos1*dcos(phi01)+IAl*dsin(Phi01)))
c     .                                          *dlog(QSTSB/Q2)
c     .    +13d0*G1Q/9d0*M1*dcos(PhiM1)*dlog(MAX(M1**2,QSTSB)/Q2)
c     .    +3d0*G2Q*M2*dcos(PhiM2)*dlog(MAX(M2**2,QSTSB)/Q2)
c     .    +64d0*PI*ALSQ/3d0*M3*dcos(PhiM3)*dlog(MAX(M3**2,QSTSB)/Q2))
c     .             /16d0/Pi**2
     
c      aux2=AT*dsin(PhiAT)+
c     .   ((6d0*YTQ**2*AT*dsin(PhiAt)+YBQ**2*AB*dsin(PhiAB)
c     .      +LQ**2*(IAl*dcos(phi01)-Alcos1*dsin(Phi01)))
c     .                                          *dlog(QSTSB/Q2)
c     .    +13d0*G1Q/9d0*M1*dsin(PhiM1)*dlog(MAX(M1**2,QSTSB)/Q2)
c     .    +3d0*G2Q*M2*dsin(PhiM2)*dlog(MAX(M2**2,QSTSB)/Q2)
c     .    +64d0*PI*ALSQ/3d0*M3*dsin(PhiM3)*dlog(MAX(M3**2,QSTSB)/Q2))
c     .             /16d0/Pi**2

c      IF(dabs(aux1).ge.dabs(aux2)*1.d-10)THEN
c       phiATQ=datan(aux2/aux1)
c       IF(aux1.le.0.d0)PhiATQ=PhiATQ+Pi
c      ELSEIF(aux2.ge.0.d0)THEN
c       phiATQ=Pi/2.d0
c      ELSE
c       phiATQ=-Pi/2.d0
c      ENDIF

c        ATP=dsqrt(aux1**2+aux2**2)

c      aux1=AB*dcos(PhiAB)+
c     .   ((6d0*YBQ**2*AB*dcos(PhiAB)+YTQ**2*AT*dcos(PhiAT)
c     .      +LQ**2*(Alcos1*dcos(phi01)+IAl*dsin(Phi01)))
c     .                                          *dlog(QSTSB/Q2)
c     .    +7d0*G1Q/9d0*M1*dcos(PhiM1)*dlog(MAX(M1**2,QSTSB)/Q2)
c     .    +3d0*G2Q*M2*dcos(PhiM2)*dlog(MAX(M2**2,QSTSB)/Q2)
c     .    +64d0*PI*ALSQ/3d0*M3*dcos(PhiM3)*dlog(MAX(M3**2,QSTSB)/Q2))
c     .             /16d0/Pi**2
     
c      aux2=AB*dsin(PhiAB)+
c     .   ((6d0*YBQ**2*AB*dsin(PhiAB)+YTQ**2*AT*dsin(PhiAT)
c     .      +LQ**2*(IAl*dcos(phi01)-Alcos1*dsin(Phi01)))
c     .                                          *dlog(QSTSB/Q2)
c     .    +7d0*G1Q/9d0*M1*dsin(PhiM1)*dlog(MAX(M1**2,QSTSB)/Q2)
c     .    +3d0*G2Q*M2*dsin(PhiM2)*dlog(MAX(M2**2,QSTSB)/Q2)
c     .    +64d0*PI*ALSQ/3d0*M3*dsin(PhiM3)*dlog(MAX(M3**2,QSTSB)/Q2))
c     .             /16d0/Pi**2

c      IF(dabs(aux1).ge.dabs(aux2)*1.d-10)THEN
c       phiABQ=datan(aux2/aux1)
c       IF(aux1.le.0.d0)PhiABQ=PhiABQ+Pi
c      ELSEIF(aux2.ge.0.d0)THEN
c       phiABQ=Pi/2.d0
c      ELSE
c       phiABQ=-Pi/2.d0
c      ENDIF

c        ABP=dsqrt(aux1**2+aux2**2)


c*********New soft sfermion running************************
       CALL GETSUSYCOUP(PAR,IFAIL)

       XI= g1S*(MHuS-MHdS+MQ3+2d0*MQ-2d0*(MU3+2d0*MU1)+MD3+2d0*MD
     .   +ME3+2d0*ME-ML3+2d0*ML)
       HT2=HTOPS**2
       HB2=HBOTS**2
       HTAU2=HTAUS**2
       MH1=MHuS
       MH2=MHdS
       L2=PAR(1)**2
       K2=PAR(2)**2

      SP= HT2*(-3d0*MH1 - MQ3 + 4d0*MU3)
     .  + HB2*(3d0*MH2 - MQ3 - 2d0*MD3)
     .  + HTAU2*(MH2 + ML3 - 2d0*ME3) + L2*(MH2 - MH1)
     .  + (G1S/18d0 + 3d0/2d0*G2S + 8d0/3d0*G3S)*(MQ3+2d0*MQ)
     .  - (16d0/9d0*G1S + 16d0/3d0*G3S)*(MU3+2d0*MU1)
     .  + (2d0/9d0*G1S + 8d0/3d0*G3S)*(MD3+2d0*MD)
     .  + (G1S/2d0 + 3d0/2d0*G2S)*(MH1-MH2-(ML3+2d0*ML)) 
     .  + 2d0*G1S*(ME3+2d0*ME)
      SIG1= G1S*(MH1 + MH2 + (MQ3+2d0*MQ)/3d0 + 8d0/3d0*(MU3+2d0*MU1) 
     .    + 2d0/3d0*(MD3+2d0*MD)+ ML3+2d0*ML + 2d0*(ME3+2d0*ME))
      SIG2=G2S*(MH1 + MH2 + 3d0*(MQ3+2d0*MQ) + ML3+2d0*ML)
      SIG3=G3S*(2d0*(MQ3+2d0*MQ) + MU3+2d0*MU1 + MD3+2d0*MD)

      F20= HT2*(MH1+MQ3+MU3+AT**2) + HB2*(MH2+MQ3+MD3+AB**2)
     .       - g1S*M1**2/9d0 - 3d0*g2S*M2**2
     .       - 16d0/3d0*g3S*M3**2 + XI/6d0
     .       + (-10d0*HT2**2*(MH1+MQ3+MU3+2d0*AT**2)
     .       - 10d0*HB2**2*(MH2+MQ3+MD3+2d0*AB**2)
     .       - 10d0*HB2**4*(MH2+MQ3+MD3+2d0*AB**2)
     .       - HB2*HTAU2*(2d0*MH2+MQ3+MD3+ML3+ME3+AB**2+ATAU**2
     .                         +2d0*AB*ATAU*dcos(PhiAb-PhiAtau))
     .       - L2*HT2*(2d0*MH1+MH2+MSS+MQ3+MU3+AT**2+Alcos1**2+IAl**2
     .                    +2d0*AT*(Alcos1*dcos(PhiAT)+IAl*dsin(PhiAT)))
     .       - L2*HB2*(MH1+2d0*MH2+MSS+MQ3+MD3+AB**2+Alcos1**2+IAl**2
     .                    +2d0*AB*(Alcos1*dcos(PhiAB)+IAl*dsin(PhiAB)))
     .       + 4d0/3d0*G1S*HT2*(MH1+MQ3+MU3+AT**2+2d0*M1**2
     .                          -2d0*M1*AT*dcos(PhiAT-PhiM1))
     .       + 2d0/3d0*G1S*HB2*(MH2+MQ3+MD3+AB**2+2d0*M1**2
     .                          -2d0*M1*AB*dcos(PhiAB-PhiM1))
     .       + 199d0/54d0*G1S**2*M1**2 + 33d0/2d0*G2S**2*M2**2
     .     - 64d0/3d0*G3S**2*M3**2
     .     + 1d0/3d0*G1S*G2S*(M1**2+M2**2+M1*M2*dcos(PhiM1-PhiM2))
     .       + 16d0/27d0*G1S*G3S*(M1**2+M3**2+M1*M3*dcos(PhiM1-PhiM3))
     .     + 16d0*G2S*G3S*(M2**2+M3**2+M2*M3*dcos(PhiM1-PhiM2))
     .       + 1d0/3d0*G1S*SP + 1d0/18d0*G1S*SIG1
     .       + 3d0/2d0*G2S*SIG2 + 8d0/3d0*G3S*SIG3)/16d0/Pi**2

      F21= 2d0*HT2*(MH1+MQ3+MU3+AT**2)
     .       - 16d0/9d0*g1S*M1**2 - 16d0/3d0*g3S*M3**2 - 2d0*XI/3d0
     .       + (-16d0*HT2**2*(MH1+MQ3+MU3+2d0*AT**2)
     .       - 2d0*HT2*HB2*(MH1+MH2+2d0*MQ3+MU3+MD3+AT**2+AB**2
     .                           +2d0*AT*AB*dcos(PhiAT-PhiAB))
     .       - 2d0*HT2*L2*(2d0*MH1+MH2+MSS+MQ3+MU3+AT**2+Alcos1**2
     .             +IAL**2+2d0*AT*(Alcos1*dcos(PhiAT)+IAl*dsin(PhiAT)))
     .       - 2d0/3d0*G1S*HT2*(MH1+MQ3+MU3+AT**2+2d0*M1**2
     .                          -2d0*AT*M1*dcos(PhiAT-PhiM1))
     .       + 6d0*G2S*HT2*(MH1+MQ3+MU3+AT**2+2d0*M2**2
     .                          -2d0*AT*M2*dcos(PhiAT-PhiM2))
     .       + 1712d0/27d0*G1S**2*M1**2 - 64d0/3d0*G3S**2*M3**2
     .       + 256d0/27d0*G1S*G3S*(M1**2+M3**2+M1*M3*dcos(PhiM1-PhiM3))
     .     - 4d0/3d0*G1S*SP + 8d0/9d0*G1S*SIG1 + 8d0/3d0*G3S*SIG3)
     .                                /16d0/Pi**2

      F22= 2d0*HB2*(MH2+MQ3+MD3+AB**2)
     .       - 4d0/9d0*g1S*M1**2 - 16d0/3d0*g3S*M3**2 + XI/3d0
     .       + (-16d0*HB2**2*(MH2+MQ3+MD3+2d0*AB**2)
     .       - 2d0*HB2*HT2*(MH1+MH2+2d0*MQ3+MU3+MD3+AB**2+AT**2
     .                  +2d0*AT*AB*dcos(PhiAB-PhiAT))
     .       - 2d0*HB2*HTAU2*(2d0*MH2+MQ3+MD3+ML3+ME3+AB**2+ATAU**2
     .                  +2d0*AB*ATAU*dcos(PhiAB-PhiAtau))
     .       - 2d0*HB2*L2*(MH1+2d0*MH2+MSS+MQ3+MD3+AB**2+Alcos1**2
     .            +IAL**2+2d0*AB*(Alcos1*dcos(PhiAB)+IAl*dsin(PhiAB)))
     .       + 2d0/3d0*G1S*HB2*(MH2+MQ3+MD3+AB**2+2d0*M1**2
     .                            -2d0*M1*AB*dcos(PhiAB-PhiM1))
     .       + 6d0*G2S*HB2*(MH2+MQ3+MD3+AB**2+2d0*M2**2
     .                            -2d0*M2*AB*dcos(PhiAB-PhiM2))
     .       + 404d0/27d0*G1S**2*M1**2 - 64d0/3d0*G3S**2*M3**2
     .       + 64d0/27d0*G1S*G3S*(M1**2+M3**2+M1*M3*dcos(PhiM1-PhiM3))
     .     + 2d0/3d0*G1S*SP + 2d0/9d0*G1S*SIG1 + 8d0/3d0*G3S*SIG3)
     .                                /16d0/Pi**2

        MSQ3=MQ3-dlog(Q2/QSTSB)*F20/16d0/Pi**2
        MSU3=MU3-dlog(Q2/QSTSB)*F21/16d0/Pi**2
        MSD3=MD3-dlog(Q2/QSTSB)*F22/16d0/Pi**2

      aux1=AT*dcos(PhiAT)+
     .   (6d0*HT2*AT*dcos(PhiAt)+HB2*AB*dcos(PhiAb)
     .     +L2*(Alcos1*dcos(phi01)+IAl*dsin(Phi01))
     .     +13d0/9d0*g1S*M1*dcos(PhiM1)+3d0*g2S*M2*dcos(PhiM2)
     .     +16d0/3d0*g3S*M3*dcos(PhiM3)
     .     +(-44d0*HT2**2*AT*dcos(PhiAt)
     .       -5d0*HT2*HB2*(AT*dcos(PhiAt)+AB*dcos(PhiAb))
     .       -3d0*HT2*L2
     .         *(AT*dcos(PhiAt)+Alcos1*dcos(phi01)+IAl*dsin(Phi01))
     .       -10d0*HB2**2*AB*dcos(PhiAb)
     .       -HB2*HTAU2*(AB*dcos(PhiAb)+ATAU*dcos(PhiAtau))
     .       -4d0*HB2*L2
     .         *(AB*dcos(PhiAb)+Alcos1*dcos(phi01)+IAl*dsin(Phi01))
     .       -6d0*L2**2*(Alcos1*dcos(phi01)+IAl*dsin(Phi01))
     .       -L2*HTAU2
     .         *(Alcos1*dcos(phi01)+IAl*dsin(Phi01)+ATAU*dcos(PhiAtau))
     .       -2d0*L2*K2*(Alcos1*dcos(phi01)+IAl*dsin(Phi01)
     .                     +Akcos2*dcos(phi02)+IAk*dsin(Phi02))
     .       +2d0*g1S*HT2*(AT*dcos(PhiAt)-M1*dcos(PhiM1))
     .       +6d0*g2S*HT2*(AT*dcos(PhiAt)-M2*dcos(PhiM2))
     .       +16d0*g3S*HT2*(AT*dcos(PhiAt)-M3*dcos(PhiM3))
     .       +2d0/3d0*g1S*HB2*(AB*dcos(PhiAb)-M1*dcos(PhiM1))
     .       -2743d0/81d0*g1S**2*M1*dcos(PhiM1)
     .       -5d0/3d0*g1S*g2S*(M1*dcos(PhiM1)+M2*dcos(PhiM2))
     .       -136d0/27d0*g1S*g3S*(M1*dcos(PhiM1)+M3*dcos(PhiM3))
     .       -15d0*g2S**2*M2*dcos(PhiM2)
     .       -8d0*g2S*g3S*(M2*dcos(PhiM2)+M3*dcos(PhiM3))
     .       +32d0/9d0*g3S**2*M3*dcos(PhiM3))/16d0/Pi**2)
     .             /16d0/Pi**2*dlog(QSTSB/Q2)
     
      aux2=AT*dsin(PhiAT)+
     .   (6d0*HT2*AT*dsin(PhiAt)+HB2*AB*dsin(PhiAB)
     .     +L2*(IAl*dcos(phi01)-Alcos1*dsin(Phi01))                                   
     .     +13d0/9d0*g1S*M1*dsin(PhiM1)+3d0*g2S*M2*dsin(PhiM2)
     .     +16d0/3d0*g3S*M3*dsin(PhiM3)
     .     +(-44d0*HT2**2*AT*dsin(PhiAt)
     .       -5d0*HT2*HB2*(AT*dsin(PhiAt)+AB*dsin(PhiAb))
     .       -3d0*HT2*L2
     .         *(AT*dsin(PhiAt)+IAl*dcos(phi01)-Alcos1*dsin(Phi01))
     .       -10d0*HB2**2*AB*dsin(PhiAb)
     .       -HB2*HTAU2*(AB*dsin(PhiAb)+ATAU*dsin(PhiAtau))
     .       -4d0*HB2*L2
     .         *(AB*dsin(PhiAb)+IAl*dcos(phi01)-Alcos1*dsin(Phi01))
     .       -6d0*L2**2*(IAl*dcos(phi01)-Alcos1*dsin(Phi01))
     .       -L2*HTAU2
     .         *(IAl*dcos(phi01)-Alcos1*dsin(Phi01)+ATAU*dsin(PhiAtau))
     .       -2d0*L2*K2*(IAl*dcos(phi01)-Alcos1*dsin(Phi01)
     .                     +IAk*dcos(phi02)-Akcos2*dsin(Phi02))
     .       +2d0*g1S*HT2*(AT*dsin(PhiAt)-M1*dsin(PhiM1))
     .       +6d0*g2S*HT2*(AT*dsin(PhiAt)-M2*dsin(PhiM2))
     .       +16d0*g3S*HT2*(AT*dsin(PhiAt)-M3*dsin(PhiM3))
     .       +2d0/3d0*g1S*HB2*(AB*dsin(PhiAb)-M1*dsin(PhiM1))
     .       -2743d0/81d0*g1S**2*M1*dsin(PhiM1)
     .       -5d0/3d0*g1S*g2S*(M1*dsin(PhiM1)+M2*dsin(PhiM2))
     .       -136d0/27d0*g1S*g3S*(M1*dsin(PhiM1)+M3*dsin(PhiM3))
     .       -15d0*g2S**2*M2*dsin(PhiM2)
     .       -8d0*g2S*g3S*(M2*dsin(PhiM2)+M3*dsin(PhiM3))
     .       +32d0/9d0*g3S**2*M3*dsin(PhiM3))/16d0/Pi**2)
     .             /16d0/Pi**2*dlog(QSTSB/Q2)

      IF(dabs(aux1).ge.dabs(aux2)*1.d-10)THEN
       phiATQ=datan(aux2/aux1)
       IF(aux1.le.0.d0)PhiATQ=PhiATQ+Pi
      ELSEIF(aux2.ge.0.d0)THEN
       phiATQ=Pi/2.d0
      ELSE
       phiATQ=-Pi/2.d0
      ENDIF

        ATP=dsqrt(aux1**2+aux2**2)

      aux1=AB*dcos(PhiAB)+
     .   (6d0*HB2**2*AB*dcos(PhiAB)+HT2*AT*dcos(PhiAT)
     .      +L2*(Alcos1*dcos(phi01)+IAl*dsin(Phi01))
     .      +7d0/9d0*G1S*M1*dcos(PhiM1)+3d0*G2S*M2*dcos(PhiM2)
     .      +16d0/3d0*g3S*M3*dcos(PhiM3)
     .      +(-44d0*HB2**2*AB*dcos(PhiAB)
     .      -5d0*HT2*HB2*(AT*dcos(PhiAt)+AB*dcos(PhiAB))
     .      -3d0*HB2*HTAU2*(AB*dcos(PhiAB)+ATAU*dcos(PhiAtau))
     .      -3d0*HB2*L2
     .           *(AB*dcos(PhiAB)+Alcos1*dcos(phi01)+IAl*dsin(Phi01))
     .      -10d0*HT2**2*AT*dcos(PhiAT)
     .      -4d0*HT2*L2
     .           *(AT*dcos(PhiAT)+Alcos1*dcos(phi01)+IAl*dsin(Phi01))
     .      -6d0*HTAU2**2*ATAU*dcos(PhiAtau)
     .      -6d0*L2**2*(Alcos1*dcos(phi01)+IAl*dsin(Phi01))
     .      -2d0*L2*K2*(Alcos1*dcos(phi01)+IAl*dsin(Phi01)
     .                     +Akcos2*dcos(phi02)+IAk*dsin(Phi02))
     .      +2d0/3d0*g1S*HB2*(AB*dcos(PhiAB)-M1*dcos(PhiM1))
     .      +6d0*g2S*HB2*(AB*dcos(PhiAB)-M2*dcos(PhiM2))
     .      +16d0*g3S*HB2*(AB*dcos(PhiAB)-M3*dcos(PhiM3))
     .      +4d0/3d0*g1S*HT2*(AT*dcos(PhiAT)-M1*dcos(PhiM1))
     .      +2d0*HTAU2*g1S*(ATAU*dcos(PhiAtau)-M1*dcos(PhiM1))
     .      -1435d0/81d0*g1S**2*M1*dcos(PhiM1)
     .      -5d0/3d0*g1S*g2S*(M1*dcos(PhiM1)+M2*dcos(PhiM2))
     .      -40d0/27d0*g1S*g3S*(M1*dcos(PhiM1)+M3*dcos(PhiM3))
     .      -15d0*g2S**2*M2*dcos(PhiM2)
     .      -8d0*g2S*g3S*(M2*dcos(PhiM2)+M3*dcos(PhiM3))
     .      +32d0/9d0*g3S**2*M3*dcos(PhiM3))/16d0/Pi**2)
     .             /16d0/Pi**2*dlog(QSTSB/Q2)
     
      aux2=AB*dsin(PhiAB)+
     .   (6d0*HB2*AB*dsin(PhiAB)+HT2*AT*dsin(PhiAT)
     .      +L2*(IAl*dcos(phi01)-Alcos1*dsin(Phi01))                                  
     .      +7d0/9d0*g1S*M1*dsin(PhiM1)+3d0*g2S*M2*dsin(PhiM2)
     .      +16d0/3d0*g3S*M3*dsin(PhiM3)
     .      +(-44d0*HB2**2*AB*dsin(PhiAB)
     .      -5d0*HT2*HB2*(AT*dsin(PhiAt)+AB*dsin(PhiAB))
     .      -3d0*HB2*HTAU2*(AB*dsin(PhiAB)+ATAU*dsin(PhiAtau))
     .      -3d0*HB2*L2
     .           *(AB*dsin(PhiAB)+IAl*dcos(phi01)-Alcos1*dsin(Phi01))
     .      -10d0*HT2**2*AT*dsin(PhiAT)
     .      -4d0*HT2*L2
     .           *(AT*dsin(PhiAT)+IAl*dcos(phi01)-Alcos1*dsin(Phi01))
     .      -6d0*HTAU2**2*ATAU*dsin(PhiAtau)
     .      -6d0*L2**2*(IAl*dcos(phi01)-Alcos1*dsin(Phi01))
     .      -2d0*L2*K2*(IAl*dcos(phi01)-Alcos1*dsin(Phi01)
     .                     +IAk*dcos(phi02)-Akcos2*dsin(Phi02))
     .      +2d0/3d0*g1S*HB2*(AB*dsin(PhiAB)-M1*dsin(PhiM1))
     .      +6d0*g2S*HB2*(AB*dsin(PhiAB)-M2*dsin(PhiM2))
     .      +16d0*g3S*HB2*(AB*dsin(PhiAB)-M3*dsin(PhiM3))
     .      +4d0/3d0*g1S*HT2*(AT*dsin(PhiAT)-M1*dsin(PhiM1))
     .      +2d0*HTAU2*g1S*(ATAU*dsin(PhiAtau)-M1*dsin(PhiM1))
     .      -1435d0/81d0*g1S**2*M1*dsin(PhiM1)
     .      -5d0/3d0*g1S*g2S*(M1*dsin(PhiM1)+M2*dsin(PhiM2))
     .      -40d0/27d0*g1S*g3S*(M1*dsin(PhiM1)+M3*dsin(PhiM3))
     .      -15d0*g2S**2*M2*dsin(PhiM2)
     .      -8d0*g2S*g3S*(M2*dsin(PhiM2)+M3*dsin(PhiM3))
     .      +32d0/9d0*g3S**2*M3*dsin(PhiM3))/16d0/Pi**2)
     .             /16d0/Pi**2*dlog(QSTSB/Q2)

      IF(dabs(aux1).ge.dabs(aux2)*1.d-10)THEN
       phiABQ=datan(aux2/aux1)
       IF(aux1.le.0.d0)PhiABQ=PhiABQ+Pi
      ELSEIF(aux2.ge.0.d0)THEN
       phiABQ=Pi/2.d0
      ELSE
       phiABQ=-Pi/2.d0
      ENDIF

        ABP=dsqrt(aux1**2+aux2**2)

!      write(0,*) "init_CPV,bottom:"
!      write(0,*) "RENSCALE",Q2,"STSBSCALE",QSTSB
!      write(0,*) "GAUGE",ALSMZ,ALEMMZ,GF,g1,g2,S2TW
!      write(0,*) "QGAUGE",G1Q,G2Q,GQ,ALSQ
!      write(0,*) "QHIGGS",ZHU,ZHD,ZS,vuq,vdq,TANBQ
!      write(0,*) "QQUARK",Ytq,Ybq,MTOPQ,MBOTQ
!      write(0,*) "QPAR",lq,kq,Alcos1,Akcos2,muq,NUQ
!      write(0,*) Ytq**2,muH2,mtopq**2,QSTSB
!      write(0,*) g1q+g2q,MZ**2,sinb**2
!      write(0,*) g2q,MW**2
!      write(0,*) g1q,mur**2,M1r**2,mur**2,M2r**2
!      write(0,*) YTQ,vuQ

      RETURN
      END
