#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: PelForcingForBen
!
! DESCRIPTION
!   !

!   This file is generated directly from OpenSesame model code, using a code
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine PelForcingForBenDynamics(mode)
!
! !USES:
  ! The following Pelagic-states are used (NOT in fluxes): R6c, R6n, R6p, R6s, &
  ! N1p, N3n, N4n, N5s, N6r, O2o
  ! The following box states are used (NOT in fluxes): PhytoPlankton
  ! The following global scalar vars are used: &
  ! BoxNumber, BoxNumberXY
  ! The following Pelagic 1-d global boxvars  are used: ETW, Depth
  ! The following Benthic 1-d global boxvars are modified : PI_Benc, &
  ! PI_Benn, PI_Benp, PI_Bens, ZI_Fc, ZI_Fn, ZI_Fp
  ! The following Benthic 1-d global boxvars got a value: RI_Fc, &
  ! RI_Fn, RI_Fp, RI_Fs, N1p_Ben, N3n_Ben, N4n_Ben, N5s_Ben, N6r_Ben, O2o_Ben, &
  ! ETW_Ben, Depth_Ben
  ! The following groupmember vars  are used: iiPhytoPlankton, iiP1
  ! The following constituent constants  are used: iiC, iiN, iiP, iiS
  ! The following global constants are used: RLEN

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN,LOGUNIT,ZERO,NZERO,DONE
  use mem, ONLY: RZc,R6c, R6n, R6p, R6s, N1p, N3n, N4n, N5s, N6r, O2o, &
   iiPhytoPlankton, ppPhytoPlankton, PhytoPlankton, Z2c, D2STATE, &
   iiMicroZooPlankton,ppMicroZooPlankton,MicroZooPlankton,iiSuspensionFeeders,&
   iiMesoZooPlankton,ppMesoZooPlankton,MesoZooPlankton
  use mem, ONLY: BoxNumber,BoxNumberXY,NO_BOXES_XY,NO_BOXES_Z,  &
    dry_z, &
    puP6Y3,EIR,ETW, ESW, ERHO, ETW_Ben, ESW_Ben, ERHO_Ben,xEPS, &
    ZE_Benc, ZE_Benn, ZE_Benp, PI_Benc, PI_Benn, PI_Benp, PI_Bens,PI_Benl, &
    Nun,Ru_Benn,R2c,R2n,R3c, &
    sediPI_Ben,sediR2_Ben,sediR6_Ben, sediPI, sediR6,sediZE_Ben,sediMeZ, &
    sediR2, &
    Depth, RI_Fc, ZI_Fc, ZI_Fn, ZI_Fp, RI_Fn, RI_Fp,sediRZ, &
    RI_Fs, N1p_Ben, N3n_Ben, N4n_Ben, N5s_Ben, N6r_Ben, O2o_Ben, ETW_Ben, &
    Depth_Ben, iiYy3,iiY3, iiC,iiN,iiP,iiS,iiL, iiP1,iiP5,iiP6,PelBoxAbove, &
    efilP6Y3,efilPART,ctfPm2c,ctfZim2c,ctfZem2c, qR2P1, &
    R2_Benc,R2_Benn ,R3_Benc ,R3_Benn ,R3_Benp,EIR_Ben,qnPc,qpPc,cmO2o,cmO2o_Ben
#ifdef INCLUDE_BENCO2
    use mem, ONLY: O3c_Ben,O3c,O3h_Ben,O3h,DIC_Ben,DIC, &
          CO3_Ben,CO3,HCO3_Ben,HCO3,CO2_Ben,CO2, ctO3m2h
#endif
  use mem_MicroZoo, ONLY:p_qnMic=>p_qnc,p_qpMic=>p_qpc
  use mem_MesoZoo, ONLY:p_qnMec=>p_qnc,p_qpMec=>p_qpc
  use mem_FilterFeeder, ONLY:p_vum

  use mem_Param,  ONLY: CalcPhytoPlankton,CalcMesoZooPlankton, &
       CalcMicroZooPlankton,p_dry_ben,p_clDxm
  use mem_Phaeo, ONLY:CALC_LIMIT_FILTERCAP,CALC_FOOD_FILTERFEEDER, &
                            CALC_FOOD_YOUNG_FILTERFEEDER
  use mem_GlobalFun, ONLY: insw,exp_limit_scalar
  use global_interface,ONLY:PhaeocystisCalc_1l

! use mem,  ONLY: Output2d_1,Output2d_2,Output2d_3, Output2d_4
!
!
! !AUTHORS
!   Piet Ruardij
!
!
! !REVISION_HISTORY
!   Created at Wed Jun 16 02:04:44 PM CEST 2004
!
!
! COPYING
!
!   Copyright (C) 2006 P. Ruardij, the mfstep group, the ERSEM team
!   (rua@nioz.nl, vichi@bo.ingv.it)
!
!   This program is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTEABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!EOP
!-------------------------------------------------------------------------!
!BOC
!
!
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  IMPLICIT NONE
  integer,intent(IN)        :: mode

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Vectors used  of group vectors
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN), dimension(:), pointer  ::lcl_Plankton
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer        :: i
  integer        :: j
  real(RLEN)     :: r
  real(RLEN),dimension(NO_BOXES_XY)     :: rr
  real(RLEN)     :: Zc,Pc,Pnp

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! user defined external functions
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN),external:: GetDelta
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  DO BoxNumberXY=1,NO_BOXES_XY
    BoxNumber=PelBoxAbove(BoxNumberXY)
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Compute total phytoplankton conc. used as food for filtereeders
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  ! normal entry: all forcing for benthic is calculated on basis of the status.
    if ( mode.eq.1) then
      ctfPm2c=ZERO
      R2_Benc=ZERO;R2_Benn=ZERO; R3_Benc=ZERO
      do i = 1 , ( iiPhytoPlankton)
        if (CalcPhytoPlankton(i)) then
          lcl_Plankton => PhytoPlankton(i,iiC)
          Pc=lcl_Plankton(BoxNumber); Pc=Pc * insw(Pc-1.0D-10)
          PI_Benc(BoxNumberXY,i)=Pc
          ctfPm2c=ctfPm2c+sum(Depth*lcl_Plankton)
          lcl_Plankton => PhytoPlankton(i,iiN)
          PI_Benn(BoxNumberXY,i) = lcl_Plankton(BoxNumber)
          lcl_Plankton => PhytoPlankton(i,iiP)
          PI_Benp(BoxNumberXY,i) = lcl_Plankton(BoxNumber)
          lcl_Plankton => PhytoPlankton(i,iiL)
          PI_Benl(BoxNumberXY,i) = lcl_Plankton(BoxNumber)
          j=ppPhytoPlankton(i,iiS)
          if ( j > 0 ) then
            lcl_Plankton => PhytoPlankton(i,iiS)
            PI_Bens(BoxNumberXY,i)  = lcl_Plankton(BoxNumber)
          else
            PI_Bens(BoxNumberXY,i)  =   ZERO
          end if
          select case  (i)
          case (iiP1)
            R2_Benc(BoxNumberXY)=R2_Benc(BoxNumberXY) &
                 +min(R2c(BoxNumber),Pc*qR2P1(BoxNumber))
            R2_Benn(BoxNumberXY)=R2_Benn(BoxNumberXY) &
                  +min(R2n(BoxNumber),PI_Benn(BoxNumberXY,i)*qR2P1(BoxNumber))
          case (iiP6)
            PI_Benn(BoxNumberXY,i) = qnPc(BoxNumber,i)*Pc
            PI_Benp(BoxNumberXY,i) = qpPc(BoxNumber,i)*Pc
            R3_Benc(BoxNumberXY)=R3c(BoxNumber)
            lcl_Plankton => PhytoPlankton(i,iiN)
            Pnp=lcl_Plankton(BoxNumber)
            r=max(ZERO,Pnp-PI_Benn(BoxNumberXY,i))
            R3_Benn(BoxNumberXY)= r
            lcl_Plankton => PhytoPlankton(i,iiP)
            Pnp=lcl_Plankton(BoxNumber)
            r=max(ZERO,Pnp-PI_Benp(BoxNumberXY,i))
            R3_Benp(BoxNumberXY)= r
          end select
        else
          PI_Benc(BoxNumberXY,i) = ZERO
          PI_Benn(BoxNumberXY,i) = ZERO
          PI_Benp(BoxNumberXY,i) = ZERO
          PI_Bens(BoxNumberXY,i) = ZERO
          if ( i==iiP6) then
            R3_Benc(BoxNumberXY)= ZERO
            R3_Benn(BoxNumberXY)= ZERO
            R3_Benp(BoxNumberXY)= ZERO
          endif
        endif
      end do
      ! if no Phaeocystis is included infood uptake for filterfeeder there is no limitation of filtering
      ! in the presence of Phaeocystis
      sediPI_Ben(BoxNumberXY,:)  =  sediPI(BoxNumber,:)

      ZI_Fc =   ZERO
      ZI_Fn =   ZERO
      ZI_Fp =   ZERO
      ctfZim2c=ZERO

      do i = 1 , iiMicroZooPlankton
        if (CalcMicroZooPlankton(i)) then
          lcl_Plankton => MicroZooPlankton(i,iiC)
          ZI_Fc(BoxNumberXY,i)  = lcl_Plankton(BoxNumber)
          ctfZim2c(BoxNumberXY)=ctfZim2c(BoxNumberXY)+sum(Depth*lcl_Plankton)
          j = ppMicroZooPlankton(i,iiN)
          if ( j> 0) then
           lcl_Plankton => MicroZooPlankton(i,iiN)
           ZI_Fn(BoxNumberXY,i)=lcl_Plankton(BoxNumber)
          else
           ZI_Fn(BoxNumberXY,i)=lcl_Plankton(BoxNumber)*p_qnMic(i)
          endif
          j = ppMicroZooPlankton(i,iiP)
          if ( j> 0) then
            lcl_Plankton => MicroZooPlankton(i,iiP)
            ZI_Fp(BoxNumberXY,i)=lcl_Plankton(BoxNumber)
          else
           ZI_Fp(BoxNumberXY,i)=lcl_Plankton(BoxNumber)*p_qpMic(i)
          endif
        endif
      enddo

      ctfZem2c=ZERO
      do i = 1 , ( iiMesoZooPlankton)
        if (CalcMesoZooPlankton(i)) then
          lcl_Plankton => MesoZooPlankton(i,iiC)
          Zc=lcl_Plankton(BoxNumber); Zc=Zc * insw(Zc-1.0D-10)
          ZE_Benc(BoxNumberXY,i)=Zc
          ctfZem2c(BoxNumberXY)=ctfZem2c(BoxNumberXY)+sum(Depth*lcl_Plankton)
          j = ppMesoZooPlankton(i,iiP)
          if ( j> 0) then
            lcl_Plankton => MesoZooPlankton(i,iiP)
            ZE_Benp(BoxNumberXY,i)  = lcl_Plankton(BoxNumber)
            lcl_Plankton => MesoZooPlankton(i,iiN)
            ZE_Benn(BoxNumberXY,i)  =  lcl_Plankton(BoxNumber)
          else
            ZE_Benn(BoxNumberXY,i) =  p_qnMec(i)* Zc
            ZE_Benp(BoxNumberXY,i) =  p_qpMec(i)* Zc
          endif
        else
          ZE_Benc(BoxNumberXY,i) = ZERO
          ZE_Benn(BoxNumberXY,i) = ZERO
          ZE_Benp(BoxNumberXY,i) = ZERO
        endif
      end do
      sediZE_Ben(BoxNumberXY,:)  =  sediMeZ(BoxNumber,:)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Compute total detritus conc. used as food for filtereeders
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      RI_Fc(BoxNumberXY)  =   max(ZERO,R6c(BoxNumber)-RZc(BoxNumber))
      r=RI_Fc(BoxNumberXY)/(NZERO+R6c(BoxNumber))
      RI_Fn(BoxNumberXY)  =   R6n(BoxNumber)*r
      RI_Fp(BoxNumberXY)  =   R6p(BoxNumber)*r
      RI_Fs(BoxNumberXY)  =   R6s(BoxNumber)*r
      sediR6_Ben(BoxNumberXY) = max(ZERO,min(sediR6(BoxNumber),&
        (sediR6(BoxNumber)*R6c(BoxNumber)-sediRZ(BoxNumber)*RZc(BoxNumber))/ &
                (NZERO+max(ZERO,R6c(BoxNumber)-RZc(BoxNumber)))))
!QZ   issue 1: TEP near-bed correction used the detritus (silt-bound) sinking rate
!     sediR2_Ben(BoxNumberXY) = sediR6(BoxNumber)
      sediR2_Ben(BoxNumberXY) = sediR2(BoxNumber)

      Ru_Benn(BoxNumberXY) =  Nun(BoxNumber)

      call FilterLimPart(BoxNumber,efilPART(BoxNumberXY))

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Derive Forcing for benthos
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      !Nutrient Forcing:
      N1p_Ben(BoxNumberXY)  =   max(NZERO,N1p(BoxNumber))
      N3n_Ben(BoxNumberXY)  =   max(NZERO,N3n(BoxNumber))
      N4n_Ben(BoxNumberXY)  =   max(NZERO,N4n(BoxNumber))
      N5s_Ben(BoxNumberXY)  =   max(NZERO,N5s(BoxNumber))
      N6r_Ben(BoxNumberXY)  =   max(NZERO,N6r(BoxNumber))

      !Oxygen Forcing:
      cmO2o_Ben(BoxNumberXY)  =    cmO2o(BoxNumber)
      if ( dry_z(BoxNumberXY).gt.p_clDxm.or.(.not.p_dry_ben)) then
        O2o_Ben(BoxNumberXY)  =   max(NZERO,O2o(BoxNumber))
      else
        O2o_Ben(BoxNumberXY)  =    cmO2o(BoxNumber)
      endif

      ! Temperature in the benthos is made equal of the temperature of the
      ! adjacent level (layer) of the pelagic
      ETW_Ben(BoxNumberXY)  =   ETW(BoxNumber)
      ESW_Ben(BoxNumberXY)  =   ESW(BoxNumber)
      ERHO_Ben(BoxNumberXY) =   ERHO(BoxNumber)
      r=DONE;if (p_dry_ben) r=dry_z(BoxNumberXY)
      EIR_Ben(BoxNumberXY)  =   EIR(BoxNumber) &
                    *exp_limit_scalar(-xEPS(BoxNumber)*Depth(BoxNumber)*r)
      if (isnan(EIR_Ben(BoxNumberXY))) then
         write(LOGUNIT,*) 'PFFB:NaN in EIR_Ben EIR(boxnumber),xEPS,Depth*r:', &
                     EIR(Boxnumber),xEPS(BoxNumber),Depth(BoxNumber),r
         EIR_Ben(BoxNumberXY)=NZERO
      endif
      ! Calculate G2_xavail_o: concentration + potential flux from pelagic
      call BenOxygenDynamics(0)
#ifdef INCLUDE_BENCO2
      O3c_Ben(BoxNumberXY)  =   O3c(BoxNumber)
      O3h_Ben(BoxNumberXY)  =   O3h(BoxNumber)
      DIC_Ben(BoxNumberXY)  =   DIC(BoxNumber)
      CO2_Ben(BoxNumberXY)  =   CO2(BoxNumber)
      CO3_Ben(BoxNumberXY)  =   CO3(BoxNumber)
      HCO3_Ben(BoxNumberXY) =   HCO3(BoxNumber)
      ctO3m2h(BoxNumberXY)  =   sum(O3h *Depth)
#endif

      ! depth of the level aboce the sediment
      Depth_Ben(BoxNumberXY)  =   Depth(BoxNumber)

      ! gr of sediment resuspended based on silt whivh is brought in resuspensnion
    endif

  end DO
  do i = 1,iiSuspensionFeeders
     select case (i)
       case (iiYy3)  ; j=CALC_FOOD_YOUNG_FILTERFEEDER
       case (iiY3)   ; j=CALC_FOOD_FILTERFEEDER
     end select
     ! output:puP6Y3 , input:DONE,DONE
     rr=DONE
     call PhaeocystisCalc_1l(j,iiP6, puP6Y3(:,i),rr,DONE)
  end do

  efilP6Y3(:)=DONE
  if (CalcPhytoPlankton(iiP6)) then
    !Calculate filterlimitation due to presence of Phaeocystis Calc
    call PhaeocystisCalc_1l(CALC_LIMIT_FILTERCAP,iiP6, &
          efilP6Y3, PI_Benc(:,iiP6),p_vum(iiY3))
  endif

  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
