#if 0
* vert_quu5.F
* form factors for {F[6, {1}], F[6, {1}]} -> {V[1], V[1]}
* generated by FormCalc 8.0 on 6-Apr-2015 1:04
#endif

#include "vars.h"


	subroutine vert_quu5
	implicit none

#include "vars.h"


#include "inline.h"

        Cloop(1) = Cloop(1) + 
     &    1/Pi**2*(AAABR(373)*
     &        ((4*AAABR(185) - AAABR(371))*
     &           (-(1/18.D0*(Sub39*pave5(cc00,quu5))) - 
     &             1/36.D0*(Sub21*pave5(cc12,quu5)) + 
     &             1/72.D0*((Sub43 + MNE1*Sub44)*pave5(cc2,quu5)) + 
     &             1/36.D0*(Sub35*pave8(bb0,quu5) + 
     &                Sub38*pave8(bb1,quu5))) - 
     &          1/6.D0*(Eps1*Sub20*MTR103(quu5)*pave5(cc0,quu5)*
     &             quuMass(quu5))))/(-MZ2 + S)

#ifdef DEBUG
	DEB("vert_quu5: Cloop(1) =", Cloop(1))
#endif

	end
