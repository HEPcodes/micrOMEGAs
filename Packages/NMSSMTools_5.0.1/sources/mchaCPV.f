      SUBROUTINE MCHA_CPV(PAR,IFAIL)

c      Chargino masses: diag(MCH(I)) = U*.MCH.V+
c       The squared chargino masses MCH2(i) (i=1..2) and the rotation
c       matrices U(i,j,k), V(i,j,k) (i,j=1..2; k=1..2: real/imaginary part)
c       are stored in the common CHASPEC.
c       In case of a negative MCH2(i), IFAIL=6.

      IMPLICIT NONE
      
      INTEGER IFAIL

      DOUBLE PRECISION PAR(*),Pi
      DOUBLE PRECISION MCHI2_11,MCHI2_22,MCHI4_12,Tr,det,aux1,aux2,aux
      DOUBLE PRECISION thetaU,thetaV,PhiCHV1,PhiCHU1,PhiCHV2,PhiCHU2

      DOUBLE PRECISION GF,MZ,MW,g1,g2,alSMZ,S2TW,ALEMMZ
      DOUBLE PRECISION tanb,cosb,sinb,vu,vd
      DOUBLE PRECISION phi01,phi02,phi0,phiM1,phiM2,phiM3,      
     .              phiAT,phiAB,phiATAU,phiAC,phiAS,phiAMU
      DOUBLE PRECISION mur,M1r,M2r,msi
      DOUBLE PRECISION MCH2(2),U(2,2,2),V(2,2,2)

      COMMON/EWPAR/GF,MZ,MW,g1,g2,alSMZ,S2TW,ALEMMZ
      COMMON/TBPAR/tanb,cosb,sinb,vu,vd
      COMMON/PHASES/phi01,phi02,phi0,phiM1,phiM2,phiM3,      
     .              phiAT,phiAB,phiATAU,phiAC,phiAS,phiAMU
      COMMON/GAUGINOPAR/mur,M1r,M2r,msi
      COMMON/CHASPEC/MCH2,U,V


      PI=4d0*DATAN(1d0)

c           1) Squared mass matrix
      MCHI2_11=M2r**2+g2*vd**2
      MCHI2_22=mur**2+g2*vu**2
      MCHI4_12=g2*(M2r**2*vu**2+mur**2*vd**2+2.d0*mur*M2r*vu*vd      
     .                                            *dcos(phiM2+phi01))

c           2) Eigenvalues: squared masses
      Tr=MCHI2_11+MCHI2_22
      det=(MCHI2_11-MCHI2_22)**2+4.d0*MCHI4_12
      MCH2(1)=(Tr-dsqrt(det))/2.d0
      MCH2(2)=(Tr+dsqrt(det))/2.d0
      IF(MCH2(1).le.0.d0)THEN
c       print*,'chargino mass problem'
       IFAIL=6
       GOTO 616
      ENDIF

c           3) Rotating matrix V0 (special unitary)
      aux1=-M2r*vu*dsin(phiM2)+mur*vd*dsin(phi01)
      aux2=M2r*vu*dcos(phiM2)+mur*vd*dcos(phi01)
      IF(dabs(aux2).ge.dabs(aux1)*1.d-10)THEN
        phiCHV1=datan(aux1/aux2)
        IF(aux2.le.0.d0)PhiCHV1=PhiCHV1+Pi
      ELSEIF(aux1.ge.0.d0)THEN
       phiCHV1=Pi/2.d0
      ELSE
       phiCHV1=-Pi/2.d0
      ENDIF

      aux1=M2r**2-mur**2+g2*(vd**2-vu**2)+dsqrt(det)
      aux2=dsqrt(4.d0*MCHI4_12)
      IF(dabs(aux2).ge.dabs(aux1)*1.d-10)THEN
       thetaV=datan(aux1/aux2)
       IF(aux2.le.0.d0)thetaV=thetaV+Pi
      ELSEIF(aux1.ge.0.d0)THEN
       thetaV=Pi/2.d0
      ELSE
       thetaV=-Pi/2.d0
      ENDIF

c           3) Rotating matrix U0 (special unitary)
      aux=g2*(M2r**2*vd**2+mur**2*vu**2+2.d0*mur*M2r*vu*vd      
     .                                            *dcos(phiM2+phi01))
      aux1=-M2r*vd*dsin(phiM2)+mur*vu*dsin(phi01)
      aux2=M2r*vd*dcos(phiM2)+mur*vu*dcos(phi01)
      IF(dabs(aux2).ge.dabs(aux1)*1.d-10)THEN
       phiCHU1=datan(aux1/aux2)
       IF(aux2.le.0.d0)PhiCHU1=PhiCHU1+Pi
      ELSEIF(aux1.ge.0.d0)THEN
       phiCHU1=Pi/2.d0
      ELSE
       phiCHU1=-Pi/2.d0
      ENDIF

      aux1=M2r**2-mur**2+g2*(vu**2-vd**2)+dsqrt(det)
      aux2=dsqrt(4.d0*aux)
      IF(dabs(aux2).ge.dabs(aux1)*1.d-10)THEN
       thetaU=datan(aux1/aux2)
       IF(aux2.le.0.d0)thetaU=thetaU+Pi
      ELSEIF(aux1.ge.0.d0)THEN
       thetaU=Pi/2.d0
      ELSE
       thetaU=-Pi/2.d0
      ENDIF

c           4) Rotating matrices U,V (phase-shIFted)
      aux1=M2r*dsin(thetaU)*dsin(thetaV)*dsin(phiCHV1+phiCHU1+phiM2)      
     .      +mur*dcos(thetaU)*dcos(thetaV)*dsin(phi01)      
     .      +dsqrt(g2)*vu*dsin(thetaU)*dcos(thetaV)*dsin(phiCHU1)      
     .      +dsqrt(g2)*vd*dsin(thetaV)*dcos(thetaU)*dsin(phiCHV1)
      aux2=M2r*dsin(thetaU)*dsin(thetaV)*dcos(phiCHV1+phiCHU1+phiM2)      
     .      +mur*dcos(thetaU)*dcos(thetaV)*dcos(phi01)      
     .      +dsqrt(g2)*vu*dsin(thetaU)*dcos(thetaV)*dcos(phiCHU1)      
     .      +dsqrt(g2)*vd*dsin(thetaV)*dcos(thetaU)*dcos(phiCHV1)
      IF(dabs(aux2).ge.dabs(aux1)*1.d-10)THEN
       phiCHV2=datan(aux1/aux2)
        IF(aux2.le.0.d0)PhiCHV2=PhiCHV2+Pi
      ELSEIF(aux1.ge.0.d0)THEN
       phiCHV2=Pi/2.d0
      ELSE
       phiCHV2=-Pi/2.d0
      ENDIF

      aux1=M2r*dcos(thetaU)*dcos(thetaV)*dsin(phiM2)      
     .    +mur*dsin(thetaU)*dsin(thetaV)*dsin(phi01-phiCHV1-phiCHU1)      
     .    +dsqrt(g2)*vu*dcos(thetaU)*dsin(thetaV)*dsin(phiCHV1)      
     .    +dsqrt(g2)*vd*dcos(thetaV)*dsin(thetaU)*dsin(phiCHU1)
      aux2=M2r*dcos(thetaU)*dcos(thetaV)*dcos(phiM2)      
     .    +mur*dsin(thetaU)*dsin(thetaV)*dcos(phiCHV1+phiCHU1-phi01)      
     .    -dsqrt(g2)*vu*dcos(thetaU)*dsin(thetaV)*dcos(phiCHV1)      
     .    -dsqrt(g2)*vd*dcos(thetaV)*dsin(thetaU)*dcos(phiCHU1)
      IF(dabs(aux2).ge.dabs(aux1)*1.d-10)THEN
       phiCHU2=datan(aux1/aux2)
        IF(aux2.le.0.d0)PhiCHU2=PhiCHU2+Pi
      ELSEIF(aux1.ge.0.d0)THEN
       phiCHU2=Pi/2.d0
      ELSE
       phiCHU2=-Pi/2.d0
      ENDIF

      U(1,1,1)=dcos(thetaU)*dcos(PhiCHU2)
      U(1,1,2)=dcos(thetaU)*dsin(PhiCHU2)
      U(1,2,1)=-dsin(thetaU)*dcos(PhiCHU2+PhiCHU1)
      U(1,2,2)=-dsin(thetaU)*dsin(PhiCHU2+PhiCHU1)
      U(2,1,1)=dsin(thetaU)*dcos(PhiCHU1)
      U(2,1,2)=-dsin(thetaU)*dsin(PhiCHU1)
      U(2,2,1)=dcos(thetaU)
      U(2,2,2)=0.d0
c      UCOMP(1,1)=DCMPLX(U(1,1,1),U(1,1,2))
c      UCOMP(1,2)=DCMPLX(U(1,2,1),U(1,2,2))
c      UCOMP(2,1)=DCMPLX(U(2,1,1),U(2,1,2))
c      UCOMP(2,2)=DCMPLX(U(2,2,1),U(2,2,2))

      V(1,1,1)=dcos(thetaV)
      V(1,1,2)=0.d0
      V(1,2,1)=-dsin(thetaV)*dcos(PhiCHV1)
      V(1,2,2)=-dsin(thetaV)*dsin(PhiCHV1)
      V(2,1,1)=dsin(thetaV)*dcos(PhiCHV2-PhiCHV1)
      V(2,1,2)=dsin(thetaV)*dsin(PhiCHV2-PhiCHV1)
      V(2,2,1)=dcos(thetaV)*dcos(PhiCHV2)
      V(2,2,2)=dcos(thetaV)*dsin(PhiCHV2)
c      VCOMP(1,1)=DCMPLX(V(1,1,1),V(1,1,2))
c      VCOMP(1,2)=DCMPLX(V(1,2,1),V(1,2,2))
c      VCOMP(2,1)=DCMPLX(V(2,1,1),V(2,1,2))
c      VCOMP(2,2)=DCMPLX(V(2,2,1),V(2,2,2))

c      print*,'mcha',dsqrt(MCH2(1)),dsqrt(MCH2(2))
c      print*,'U11',U(1,1,1),U(1,1,2)
c      print*,'U12',U(1,2,1),U(1,2,2)
c      print*,'U21',U(2,1,1),U(2,1,2)
c      print*,'U22',U(2,2,1),U(2,2,2)
c      print*,'V11',V(1,1,1),V(1,1,2)
c      print*,'V12',V(1,2,1),V(1,2,2)
c      print*,'V21',V(2,1,1),V(2,1,2)
c      print*,'V22',V(2,2,1),V(2,2,2)

 616  RETURN
      END
