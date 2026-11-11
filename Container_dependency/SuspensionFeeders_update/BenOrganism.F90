#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: BenOrganism
!
! DESCRIPTION
!   !    This submodel describes the carbon dynamics and associated
!    nutrient dynamics in benthic organisms
!
!

!   This file is generated directly from OpenSesame model code, using a code
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine BenOrganismDynamics(y,  ppyc, ppyn, ppyp)
!
! !USES:

  ! For the following Benthic-states fluxes are defined: Q6c, Q6n, Q6p, G2o, &
  ! K4n, K1p, D6m, D7m, D8m
  ! The following Benthic-states are used (NOT in fluxes): D1m
  ! For the following Benthic-group-states fluxes are defined: BenOrganisms, &
  ! BenBacteria
  ! The following Benthic 1-d global boxvars are modified :
  ! The following Benthic 1-d global boxvars  are used: ETW_Ben
  ! The following groupmember vars are used: iiBenOrganisms, iiBenBacteria, &
  ! iiY1, iiY2, iiY4, iiY5
  ! The following constituent constants  are used: iiC, iiN, iiP
  ! The following 0-d global parameters are used: p_d_tot
  ! The following global constants are used: RLEN

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN,ZERO,DONE,NZERO
  use constants,ONLY:p_qnUc,INTEGRAL,AVERAGE,MW_C

#ifdef NOPOINTERS
  use mem,  ONLY: D2STATE
#else
  use mem, ONLY: D2STATE, Q6c, Q6n, Q6p, Q6s, D6m, D7m,G2o, &
    Dfm,Dcm,D1m,D2m,D8m,D9m,BenOrganisms, BenBacteria,SuspensionFeeders,BenPhyto
#endif
  use mem, ONLY: ppQ6c, ppQ6n, ppQ6p,ppG3c,ppO3c, ppG2o,ppO2o, ppQ1c,ppQ11c, &
    ppQun,ppQ1un, ppK1p, ppK11p, ppDfm,ppDcm,ppD6m, ppD7m, ppD8m, ppD9m, &
    iiSuspensionFeeders, iiBenBacteria,iiBenPhyto,jmYIc,jmY3c,jugYIc,jrrYIc,&
    iiH1,iiH2,iiH3,iiHN,iiBen,iiPel, iiY1,iiY2,iiY4,iiY5,iiC, iiN, iiP,iiS,iiL,&
    ppBenOrganisms, ppSuspensionFeeders, ppBenBacteria,ppBenPhyto, &
    iiBenOrganisms,max_change_per_step , &
    flux_vector, sourcesink_flux_vector,NO_BOXES_XY, dry_z,Dfm,G2_xavail_o, &
    jPIQ6s,jO2Y2o,jnetYIc,ETW_Ben,O2o_Ben,Depth_Ben,cmO2o_Ben,jBenFishInput
  use mem_Param,  ONLY: p_pe_R1c, p_pe_R1n, p_pe_R1p, &
        p_d_tot,p_poro,CalcBenPhyto,CalcBenBacteria,p_dry_ben,combine_anabac
  use LimitRates, ONLY:LimitChange_vector
  use Constants, ONLY:NEGATIVE,POSITIVE
  use botflux,onlY:addbotflux_vector,openbotflux_vector
  use mem_BenOrganism
  use mem_BenPhyto,ONLY: CalculateBenPhyto_vector
  use mem_Bioturbation,only:sw_irr,p_pu_xoxPelBen,p_xeff_oxuptake
! use mem,  ONLY: Output2d_1,Output2d_2,Output2d_3, Output2d_4
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following vector functions are used:eTq_vector, &
  ! MM_vector, PartQ_vector
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use mem_globalfun,ONLY:eTq_vector,MM_vector,PartQ_vector,insw_vector,exp_limit
  use global_interface, only: FindMidPointExp_vector



  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  IMPLICIT NONE

! !INPUT:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer,intent(IN)  :: y
  integer,intent(IN) :: ppyc
  integer,intent(IN) :: ppyn
  integer,intent(IN) :: ppyp
!
!
! !AUTHORS
!   W. Ebenhoh and C. Kohlmeier.
!
!
! !REVISION_HISTORY
!   !
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
  ! Set up Local Variable for copy of state var. object
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real,external              :: GetDelta
  real(RLEN),dimension(NO_BOXES_XY) :: yc
  real(RLEN),dimension(NO_BOXES_XY) :: yn
  real(RLEN),dimension(NO_BOXES_XY) :: yp
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer  :: i,j
  real(RLEN)                         ::scalar
  real(RLEN),dimension(NO_BOXES_XY,iiBenBacteria)  :: avail_H
  logical,dimension(NO_BOXES_XY,iiBenPhyto)  :: xBPavailable
  real(RLEN),dimension(NO_BOXES_XY)  :: clm,clmp
  real(RLEN),dimension(NO_BOXES_XY)  :: cm,cmp
  real(RLEN),dimension(NO_BOXES_XY)  :: cmm
  real(RLEN),dimension(NO_BOXES_XY)  :: et
  real(RLEN),dimension(NO_BOXES_XY)  :: eO
  real(RLEN),dimension(NO_BOXES_XY)  :: food
  real(RLEN),dimension(NO_BOXES_XY)  :: food_PIc,food_YIc,food_FFc, food_HIc
  real(RLEN),dimension(NO_BOXES_XY)  :: food_src
  real(RLEN),dimension(NO_BOXES_XY)  :: eF
  real(RLEN),dimension(NO_BOXES_XY)  :: sug
  real(RLEN),dimension(NO_BOXES_XY)  :: sug_y
  real(RLEN),dimension(NO_BOXES_XY)  :: sun
  real(RLEN),dimension(NO_BOXES_XY)  :: sunQ6
  real(RLEN),dimension(NO_BOXES_XY)  :: se_u
  real(RLEN),dimension(NO_BOXES_XY)  :: se_uQ6
  real(RLEN),dimension(NO_BOXES_XY)  :: choice,choicel
  real(RLEN),dimension(NO_BOXES_XY)  :: availQ6_c,availQ6_n,availQ6_p
  real(RLEN),dimension(NO_BOXES_XY)  :: rugc,rugn,rugp,rro
  real(RLEN),dimension(NO_BOXES_XY)  :: runc,runn,runp
  real(RLEN),dimension(NO_BOXES_XY)  :: rqt6c,rqt6n,rqt6p,rqt6s
  real(RLEN),dimension(NO_BOXES_XY)  :: rq6c,rq6n,rq6p,rq6s
  real(RLEN),dimension(NO_BOXES_XY)  :: ruYIc,ruYIn,ruYIp
  real(RLEN),dimension(NO_BOXES_XY)  :: ruBIc,ruBIn,ruBIp
  real(RLEN),dimension(NO_BOXES_XY)  :: ruPIn,ruPIp
  real(RLEN),dimension(NO_BOXES_XY)  :: ruQ6c,ruQ6n,ruQ6p
  real(RLEN),dimension(NO_BOXES_XY)  :: rrc
  real(RLEN),dimension(NO_BOXES_XY)  :: sm
  real(RLEN),dimension(NO_BOXES_XY)  :: qnc,qpc,qsc,qlc
  real(RLEN),dimension(NO_BOXES_XY)  :: renc,renn,renp
  real(RLEN),dimension(NO_BOXES_XY)  :: edry
  real(RLEN),dimension(NO_BOXES_XY)  :: pUptake_ox
  real(RLEN),dimension(NO_BOXES_XY)  :: pudiln,pudilp! phytobenthos mass(mgC/m2)
  real(RLEN),dimension(NO_BOXES_XY)  :: xfood_m
  real(RLEN),dimension(NO_BOXES_XY)  :: rx_any,px_any,cx_any,sx_any
  real(RLEN),dimension(NO_BOXES_XY)  :: rx_any_m,mx_any ! any depth/distance
  real(RLEN),dimension(NO_BOXES_XY)  :: x_foodinD1m_c,px_foodinD1
  real(RLEN),dimension(NO_BOXES_XY,iiSuspensionFeeders)  :: pYF
  real(RLEN),dimension(NO_BOXES_XY,iiBenOrganisms)  :: pYY
  real(RLEN),dimension(NO_BOXES_XY,iiBenBacteria)  :: pYH
  real(RLEN),dimension(NO_BOXES_XY,iiBenPhyto)  :: pc
  real(RLEN),dimension(NO_BOXES_XY,iiBenPhyto)  :: ruPIc
  real(RLEN),dimension(NO_BOXES_XY,iiBenPhyto)  :: rq6l

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  Copy  state var. object in local var
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  yc = D2STATE(:,ppyc)
  yn = D2STATE(:,ppyn)
  yp = D2STATE(:,ppyp)

  edry=DONE
  ! p_xdry==0 : edry=1
  ! p_xdry==1 : edry=1.1*dry_z/(0.1_dry_z) (dry_z==0: edry=0, dry_z==1: edry=1)
  if (p_dry_ben) edry=1.1*p_xdry(y)*dry_z/(0.1+dry_z) + (DONE-p_xdry(y))
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Physiological temperature and oxygen response
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  !grazing on other benthic organism is defined according parameter
  !predation-matrix (p_Yn) . With exception for grazing on Y3: there are only the young ones
  !available as food
  cm  =   p_cm(y); avail_H=DONE; !avail_H(:,iiH2:iiH3)=ZERO
  if (combine_anabac ) then
  else
    mx_any=p_d_tot
    clm=max(p_clm(y),D1m)
    cm=min(p_cm(y),D2m)
    where (cm.gt.D1m) avail_H(:,iiH2)= &
      min(DONE, PartQ_vector(D6m(:),clm,cm,      p_d_tot) / &
                PartQ_vector(D6m(:),clm,D2m,  p_d_tot))
    clm=max(p_clm(y),D2m)
    cm=min(p_cm(y),p_d_tot)
    where (cm.gt.D2m) avail_H(:,iiH3)= &
      min(DONE, PartQ_vector(D6m(:),clm,cm,      p_d_tot) / &
                PartQ_vector(D6m(:),clm,mx_any,  p_d_tot))
  endif
  et  =   eTq_vector(  ETW_Ben(:),  p_q10(y))
  !oxygen in lowest layer  determine input to aerobic layer and
  !is a better proxy for oxygen limitation. Moreover that all processes
  !depend on oxygen use much more oxygen than present in the aerobic layer.
  !p_xsenso=0:y is not limited to layer;p_senso=1;y is found oxic layer;p_senso=2:y is found on sediment
  select case (p_xsenso(y))
    case (0);eO = DONE
    case (1);eO = DONE-exp(- max(NZERO,G2o(:)/D1m/p_poro) / p_clO2o(y))
    case (2);eO = DONE-exp(- max(NZERO,O2o_Ben(:)) / p_clO2o(y))
  end select

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculate total food
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  food  =   NZERO
  x_foodinD1m_c =   NZERO

  rugc  =   ZERO
  rugn  =   ZERO
  rugp  =   ZERO

  rqt6c  =   ZERO
  rqt6n  =   ZERO
  rqt6p  =   ZERO
  rqt6s  =   ZERO

  ! For other benthic organisms:


  clmp=p_clm(y);cmp =p_cm(y) ;mx_any=cmp-clmp
  if (p_xsenso(y)>=1) then
    cmp=min(D1m,cmp)
    where(clmp>cmp)clmp=max(ZERO,cmp-mx_any)
  endif
  if (y.eq.iiY2) then
    mx_any=p_d_tot
    cmm=FindMidPointExp_vector(D6m,d1_from=D1m,d1_to=mx_any)
    cmp=4.0D+00*cmm-clmp
  endif
  xfood_m=DONE;if (p_vum(y)>ZERO) xfood_m=cmp-clmp

  pYY=ZERO; food_YIc=ZERO

  do i = 1 ,iiBenOrganisms
    clm  =   p_clm(i)
    cm  =   p_cm(i); if (p_xsenso(i).ge.1)cm=min(cm,D1m)
    food_src  = max(ZERO,BenOrganisms(i,iiC))
    px_any=max(ZERO,min(cmp,cm)-max(clmp,clm))/(cm-clm)*p_Yn(y,i)
    pYY(:,i)= px_any*p_Yn(y,i)* MM_vector( food_src* px_any,  p_clu(y))
    food_YIc  =  food_YIc + food_src*pYY(:,i)/xfood_m
    x_foodinD1m_c   =  x_foodinD1m_c  + food_src*pYY(:,i)/xfood_m &
                             *max(ZERO,min(cmp,D1m)-max(clmp,ZERO))/(cm-clm)
  enddo
  food  =   food_YIc+ food

  food_FFc=ZERO
  do i = 1 , ( iiSuspensionFeeders)
    cx_any  = max(ZERO,SuspensionFeeders(i,iiC))
    pYF(:,i)=p_Sn(y,i)* MM_vector(  cx_any*p_Sn(y,i),  p_clu(y))
    food_FFc  =   food_FFc+ cx_any/xfood_m*pYF(:,i)
    x_foodinD1m_c  =    x_foodinD1m_c+ cx_any/xfood_m*pYF(:,i)
  end do
  food  =   food_FFc+ food

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Bacteria:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  food_HIc=ZERO;pYH=ZERO
  do i = 1,iiBenBacteria
    if (CalcBenBacteria(i)) then
      scalar=min(p_clu(y),p_cluH)
      select case (i)
        case (iiH1,iiHN);clm=ZERO;cm=D1m
        case (iiH2);     clm=D1m;cm=D2m ;if (combine_anabac) cm=p_d_tot
        case (iiH3);     clm=D2m;cm=p_d_tot
      end select
      if (i.ne.iiH3.or.(i.eq.iiH3.and.(.not.combine_anabac) ) ) then
        cx_any  =   max(ZERO,BenBacteria(i,iiC))
        if (y.eq.iiY2.and.i.eq.iiH2 .and. combine_anabac) then
          px_any=(DONE-exp_limit(-max(ZERO,xfood_m-D1m)/D6m))/ &
                     (DONE-exp_limit(-(cm-D1m)/D6m))
!         write(LOGUNIT,*) px_any,xfood_m,cm,D1m,D6m
        else
          px_any=max(ZERO,min(cmp,cm)-max(clmp,clm))/(cm-clm)
        endif
        pYH(:,i)  =px_any/xfood_m*p_Hn(y,i)*MM_vector(cx_any*px_any, scalar )
        food_HIc  =   food_HIc +cx_any* pYH(:,i)
        x_foodinD1m_c=x_foodinD1m_c +  &
          cx_any*pYH(:,i)*max(ZERO,min(cmp,D1m)-max(clmp,ZERO))/(cm-clm)
      endif
    endif
  end do
  food      =  food+ food_HIc

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Benthic Phytoplankton:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  food_PIc=ZERO;pc=ZERO

  xBPavailable= .false.
  do i = 1 ,  iiBenPhyto
    if (CalcBenPhyto(i)) then
      cx_any  =   BenPhyto(i,iiC)
      xBPavailable(:,i)= cx_any > 1.0D-6.and.p_Pn(y,i)>NZERO
      where (xBPavailable(:,i))
         px_any= CalculateBenPhyto_vector(iiC,INTEGRAL,clmp,cmp)
         pc(:,i) =p_Pn(y,i)*MM_vector(px_any*p_Pn(y,i)*cx_any,  p_clu(y))
         food_PIc  =   food_PIc + px_any*cx_any*pc(:,i)/xfood_m
         x_foodinD1m_c  =    x_foodinD1m_c + px_any*cx_any*pc(:,i)/xfood_m
      endwhere
    endif
  end do
  food      =  food+ food_PIc

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Detritus (if eaten) First calculate the available portion
  ! and then add it to the total amount of food
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  if ( p_puQ6(y)> ZERO) then
    availQ6_c  =   Q6c(:)* PartQ_vector(  D6m(:),  clm,  cm,  p_d_tot)
    availQ6_n  =   Q6n(:)* PartQ_vector(  D7m(:),  clm,  cm,  p_d_tot)
    availQ6_p  =   Q6p(:)* PartQ_vector(  D8m(:),  clm,  cm,  p_d_tot)

    food_src  =   p_puQ6(y)* availQ6_c
    food  =   food+ food_src/xfood_m* MM_vector(  food_src,  p_clu(y))
  else
    availQ6_c  =   ZERO
    availQ6_n  =   ZERO
    availQ6_p  =   ZERO
  end if

  px_foodinD1=min(DONE,x_foodinD1m_c/(NZERO+food))

  if (p_vum(y)==ZERO) then
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Simple monod....
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Correct for too much food...
    eF  =   MM_vector( food,  p_chu(y))
    ! Growth rate at actual amount with correction of growth rate for environmental factors
    sug_y  =   p_su(y)* et* eO* eF *edry

  else
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    !  Holling...........
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    sug_y= eO*et*edry* p_su(y)* MM_vector( p_vum(y)* food,p_su(y))
  endif

  !if maintenance costs for searching is defined there is only eaten
  !when the gain of searching is larger than the costs:
  if ( p_sra(y) > ZERO) then
     sug_y=sug_y*insw_vector(sug_y*( DONE- p_pue(y))*(DONE-p_pur(y))-p_sra(y))

  endif

  !Correction for low oxygen......
  rro=yc*(sug_y*(DONE-p_pue(y))*p_pur(y)+p_sra(y))*insw_vector(sug_y)/MW_C
  !-----------------------------------------------------------------------------
  if (isnan(rro(1))) then
    write(LOGUNIT,*) 'y,yc,rro', y,yc,rro,food,sug_y,eO,eF,edry,et
    write(LOGUNIT,*) 'food_PIc,food_YIc,food_FFc food_HIc', &
                                 food_PIc,food_YIc,food_FFc, food_HIc
  endif
  !-----------------------------------------------------------------------------
  !reset of var p,.
  cx_any=G2_xavail_o
  if (p_xsenso(y)==2) cx_any=O2o_Ben*Depth_Ben
  call LimitChange_vector(POSITIVE,rro,cx_any,max_change_per_step,px_any)
  pUptake_ox=ZERO
  if ((y==iiY2.or.y==iiY4).and.sw_irr>1 )then
    !The parameters in this if-endif part are defined in bioturbation.nml
    pUptake_ox=min(p_pu_xoxPelBen*px_any,DONE-px_any)* &
                                            min(DONE,O2o_Ben/cmO2o_Ben)
    px_any=px_any+p_xeff_oxuptake*pUptake_ox
  endif
  sug_y=max(ZERO,min(sug_y, &
      (px_any*rro*MW_C/(NZERO+yc)-p_sra(y))/((DONE-p_pue(y))*p_pur(y))))

  sug  =max(ZERO,( sug_y* yc/xfood_m)/ ( NZERO + food))
  jugYIc(:,y)=sug_y*yc

  ! Net uptake:
  sun    =   sug*( DONE- p_pue(y))
  sunQ6  =   sug*( DONE- p_pueQ6(y))

  ! Execreted part:
  se_u    =   sug- sun
  se_uQ6  =   sug- sunQ6


  if (isnan(sug(1))) then
    write(LOGUNIT,*) 'y,sug', y,sug,food,sug_y,sun,sunQ6,px_any,G2_xavail_o
  endif
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! Calculation of uptake rate:
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   pudiln=(DONE-p_pe_R1n)/(DONE-p_pe_R1c)
   pudilp=(DONE-p_pe_R1p)/(DONE-p_pe_R1c)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Benthic organisms:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  do i = 1 , iiBenOrganisms
    choice = pYY(:,i)
    cx_any  = BenOrganisms(i,iiC)
    qnc  = max(ZERO,BenOrganisms(i,iiN))/ (NZERO+cx_any)
    qpc  = max(ZERO,BenOrganisms(i,iiP))/ (NZERO+cx_any)
    ruYIc  = cx_any* sug* choice
    ruYIn  = qnc* ruYIc
    ruYIp  = qpc* ruYIc

    ! In case of cannibalism rate of change in state = zero!
    if ( i/= y) then
      call flux_vector( iiBen, ppBenOrganisms(i,iiC),ppyc, ruYIc )
      call flux_vector( iiBen, ppBenOrganisms(i,iiN),ppyn, ruYIn )
      call flux_vector( iiBen, ppBenOrganisms(i,iiP),ppyp, ruYIp )
    elseif ( y==iiY1 .or. y==iiY5) then
     ! grazing within a functional group is considered as 'closure' of the model
     ! This grazing can be seen therefor seen as food available for the higher
     ! trophic levels.  The higher trophic levels (=demersel fish) eatonly
     ! of the 2 top predators Y1 and Y5
      jBenFishInput=jBenFishInput + ruYIc
    end if

    rugc  =   rugc+ ruYIc
    rugn  =   rugn+ ruYIn
    rugp  =   rugp+ ruYIp
    jmYIc(:,i)=jmYIc(:,i)+ruYIc

    rq6c  =   cx_any* se_u* choice
    rq6n  =   rq6c * qnc * pudiln
    rq6p  =   rq6c * qpc * pudilp

    rqt6c  =   rqt6c+ rq6c
    rqt6n  =   rqt6n+ rq6n
    rqt6p  =   rqt6p+ rq6p

  end do

  do i = 1, iiSuspensionFeeders
    choice = pYF(:,i)
    cx_any  = SuspensionFeeders(i,iiC)
    qnc =max(ZERO, SuspensionFeeders(i,iiN)/ (NZERO+cx_any))
    qpc =max(ZERO, SuspensionFeeders(i,iiP)/ (NZERO+cx_any))
    ruYIc  = cx_any* sug* choice
    ruYIn  = qnc* ruYIc
    ruYIp  = qpc* ruYIc
    jmY3c(:,i)=jmY3c(:,i)+ruYIc

    call flux_vector( iiBen, ppSuspensionFeeders(i,iiC),ppyc, ruYIc )
    call flux_vector( iiBen, ppSuspensionFeeders(i,iiN),ppyn, ruYIn )
    call flux_vector( iiBen, ppSuspensionFeeders(i,iiP),ppyp, ruYIp )

    rugc  =   rugc+ ruYIc
    rugn  =   rugn+ ruYIn
    rugp  =   rugp+ ruYIp

    rq6c  =   cx_any* se_u* choice
    rq6n  =   rq6c * qnc * pudiln
    rq6p  =   rq6c * qpc * pudilp

    rqt6c  =   rqt6c+ rq6c
    rqt6n  =   rqt6n+ rq6n
    rqt6p  =   rqt6p+ rq6p
  end do

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Bacteria:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  do i = 1,iiBenBacteria
    if (CalcBenBacteria(i)) then
      choice  = pYH(:,i)
      cx_any=BenBacteria(i,iiC)
      qnc =max(ZERO, BenBacteria(i,iiN)/ (NZERO+cx_any))
      qpc =max(ZERO, BenBacteria(i,iiP)/ (NZERO+cx_any))
      ruBIc  =   max(ZERO,cx_any)* sug* choice
      call LimitChange_vector(POSITIVE,ruBIc,cx_any,max_change_per_step,px_any)
      ruBIc=ruBIc*px_any
      choice=choice*px_any
      ruBIn  = qnc* ruBIc
      ruBIp  = qpc* ruBIc

      rugc  =   rugc+ ruBIc
      rugn  =   rugn+ ruBIn
      rugp  =   rugp+ ruBIp

      rq6c  =   cx_any* se_u* choice
      rq6n  =   rq6c * qnc * pudiln
      rq6p  =   rq6c * qpc * pudilp

      rqt6c  =   rqt6c+ rq6c
      rqt6n  =   rqt6n+ rq6n
      rqt6p  =   rqt6p+ rq6p
    else
      ruBIc=ZERO;ruBIn=ZERO;ruBIp=ZERO;
    endif
    call flux_vector( iiBen, ppBenBacteria(i,iiC),ppyc, ruBIc )
    call flux_vector( iiBen, ppBenBacteria(i,iiN),ppyn, ruBIn )
    call flux_vector( iiBen, ppBenBacteria(i,iiP),ppyp, ruBIp )
  end do
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Benthic Phytoplankton
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  rx_any=ZERO;choice=ZERO;choicel=ZERO
  do i =1,iiBenPhyto
    where (xBPavailable(:,i))
      choice =  pc(:,i)* CalculateBenPhyto_vector(iiC,INTEGRAL,clmp,cmp)
      choicel=  pc(:,i)* CalculateBenPhyto_vector(iiL,INTEGRAL,clmp,cmp)
      cx_any=BenPhyto(i,iiC)
      qnc =max(ZERO, BenPhyto(i,iiN)/ (NZERO+cx_any))*choicel/(NZERO+choice)
      qpc =max(ZERO, BenPhyto(i,iiP)/ (NZERO+cx_any))*choicel/(NZERO+choice)
      qsc =max(ZERO, BenPhyto(i,iiS)/ (NZERO+cx_any))*choicel/(NZERO+choice)
      qlc =max(ZERO, BenPhyto(i,iiL)/ (NZERO+cx_any))*choicel/(NZERO+choice)
      rx_any= max(ZERO,cx_any)* sug* choice
    endwhere
    call LimitChange_vector(POSITIVE,rx_any,cx_any, &
                                       max_change_per_step,px_any)
    choice=choice*px_any
    where (xBPavailable(:,i))
      ruPIc(:,i)=   rx_any*px_any
      ruPIn     =   ruPIc(:,i)*qnc
      ruPIp     =   ruPIc(:,i)*qpc
      rq6s      =   ruPIc(:,i)*qsc
      rq6l(:,i) =   ruPIc(:,i)*qlc
    elsewhere
     ruPIc(:,i)=ZERO;ruPIn=ZERO;ruPIp=ZERO;rq6l(:,i)=ZERO;rq6s=ZERO
    endwhere
    jPIQ6s(:,i)=jPIQ6s(:,i)+ rq6s

    rugc  =   rugc+ ruPIc(:,i)
    rugn  =   rugn+ ruPIn
    rugp  =   rugp+ ruPIp

    rq6c  =   cx_any* se_u* choice
    rq6n  =   rq6c * qnc * pudiln
    rq6p  =   rq6c * qpc * pudilp

    rqt6c  =   rqt6c+ rq6c
    rqt6n  =   rqt6n+ rq6n
    rqt6p  =   rqt6p+ rq6p
    rqt6s  =   rqt6s+ rq6s

    call flux_vector( iiBen, ppBenPhyto(i,iiC),ppyc, ruPIc(:,i) )
    call flux_vector( iiBen, ppBenPhyto(i,iiN),ppyn, ruPIn )
    call flux_vector( iiBen, ppBenPhyto(i,iiP),ppyp, ruPIp )
    j=ppBenPhyto(i,iiL)
    call flux_vector( iiBen, j,j, -rq6l(:,i) )
  end do

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! detritus
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  select case ( p_puQ6(y)> ZERO)

    case( .TRUE. )
      choice  =  p_puQ6(y)* MM_vector(  p_puQ6(y)* availQ6_c/xfood_m,  p_clu(y))
      ruQ6c  =   sug* choice* availQ6_c
      ruQ6n  =   sug* choice* availQ6_n
      ruQ6p  =   sug* choice* availQ6_p

      call flux_vector( iiBen, ppQ6c,ppyc, ruQ6c )
      call flux_vector( iiBen, ppQ6n,ppyn, ruQ6n )
      call flux_vector( iiBen, ppQ6p,ppyp, ruQ6p )
      rugc  =   rugc+ ruQ6c
      rugn  =   rugn+ ruQ6n
      rugp  =   rugp+ ruQ6p

      rq6c  =   se_uQ6* choice* availQ6_c
      rq6n  =   se_uQ6* choice* availQ6_n* pudiln
      rq6p  =   se_uQ6* choice* availQ6_p* pudilp

      rqt6c  =   rqt6c+ rq6c
      rqt6n  =   rqt6n+ rq6n
      rqt6p  =   rqt6p+ rq6p

    case( .FALSE. )
      ruQ6c  =   ZERO
      ruQ6n  =   ZERO
      ruQ6p  =   ZERO
  end select


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculation of respiration:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  jrrYIc(:,y)=p_sr(y)* max(ZERO,yc)* et
  rrc  =jrrYIc(:,y)   + (p_pur(y)*( rugc- rqt6c)  &
                                       +p_sra(y)*yc)*insw_vector(sug)
  if (y==iiY2.or.y==iiY4) jO2Y2o=jO2Y2o+pUptake_ox*rrc/MW_C

  if (p_xsenso(y)<2) then 
    call sourcesink_flux_vector( iiBen, ppyc, ppG3c, rrc )
    call flux_vector(iiBen, ppG2o,ppG2o,-rrc/MW_C )
  else
     call addbotflux_vector(POSITIVE,iiBen,ppyc,iiPel,ppO3c,rrc)
     call openbotflux_vector(NEGATIVE,iiPel,ppO2o,-rrc/MW_C)
  endif

  ! in case of a negative value of one of the following values there is a &
  ! situation
  ! of startvation and very low biomass values. Check on quota in the food is &
  ! out of order
  runc  =   max(  ZERO,  rugc- rqt6c -rrc)
  runn  =   max(  ZERO,  rugn- rqt6n)
  runp  =   max(  ZERO,  rugp- rqt6p)

  jnetYIc(:,y)=runc
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculation of nutrient release and correction of C:N:P
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

   renc  =  max(ZERO,(-runn/( NZERO+ runc)/ p_qn(y)+DONE),   &
                    (-runp/( NZERO+ runc)/ p_qp(y)+DONE))* runc

  !Correct excretion of renn for the fact that for ecretion some C is needed
  ! for excretion as urea.
  runc=runc-renc
  renn=max(ZERO,runn-p_qn(y)*runc)/(DONE-p_qn(y)/p_qnUc)
  runc=runc-renn/p_qnUc
  renp=max(ZERO,runp-p_qp(y)*runc)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Correction for cases where initial conditions deviate strongly from
  ! Redfield C:N:P. In this way the C:N:P does not become too extreme
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  renn = min( max( ZERO, renn), max( ZERO, renn-( p_qn(y)* yc- yn)))
  renp = min( max( ZERO, renp), max( ZERO, renp-( p_qp(y)* yc- yp)))
  !urea excretion
  call flux_vector( iiBen, ppyn,ppQun, renn*px_foodinD1 )
  call flux_vector( iiBen, ppyc,ppQ1c, renn/p_qnUc*px_foodinD1 )
  call flux_vector( iiBen, ppyn,ppQ1un,renn*(DONE-px_foodinD1))
  call flux_vector( iiBen, ppyc,ppQ11c,renn/p_qnUc*(DONE-px_foodinD1))
  ! dissolved phoshpate excretion
  call flux_vector( iiBen, ppyp,ppK1p, renp*px_foodinD1 )
  call flux_vector( iiBen, ppyp,ppK11p,renp*(DONE-px_foodinD1))

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculation of mortality
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  !Only mortality if they are hungry + density dependent mortality
  !Densitiy dependent mortality is a closure  and controlled by grazers
  !hence biological process : corrected for temperature.
  ! Density p_vum>0 food_density per m3  therefor mortality also per m3
  sm  = max(ZERO,  p_sd(y)* insw_vector(NZERO-sug_y)+  &
          et*p_sd2(y)*yc/xfood_m)
  cx_any=DONE
  call LimitChange_vector(POSITIVE,sm,cx_any, max_change_per_step)

  rq6c  =   max(ZERO,yc)* sm
  rq6n  =   max(ZERO,yn)* sm
  rq6p  =   max(ZERO,yp)* sm

  rqt6c  =   rqt6c+ rq6c
  rqt6n  =   rqt6n+ rq6n
  rqt6p  =   rqt6p+ rq6p

  if ( y==iiY1 .or. y==iiY4)  jBenFishInput=jBenFishInput + rq6c

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculate total flux of Y to Q6:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  call flux_vector( iiBen, ppyc,ppQ6c, rqt6c +renc )
  call flux_vector( iiBen, ppyn,ppQ6n, rqt6n )
  call flux_vector( iiBen, ppyp,ppQ6p, rqt6p )

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Assign organism-dependent parameters for benphyt
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rx_any_m=ZERO
  do i = 1,iiBenPhyto
    where (xBPavailable(:,i))
      cx_any  =   NZERO+BenPhyto(i,iiC)
      mx_any= CalculateBenPhyto_vector(iiC,AVERAGE,clmp,cmp)
      rx_any_m=-(mx_any-Dfm(:))* ruPIc(:,i)/cx_any
    endwhere
    call LimitChange_vector(NEGATIVE,rx_any_m,Dfm,max_change_per_step)
    call flux_vector(iiBen, ppDfm,ppDfm,rx_any_m)
    where (xBPavailable(:,i))
      cx_any  =   NZERO+BenPhyto(i,iiL)
      mx_any= CalculateBenPhyto_vector(iiL,AVERAGE,clmp,cmp)
      rx_any_m=-(mx_any-Dcm(:))* rq6l(:,i)/cx_any
    endwhere
    call LimitChange_vector(NEGATIVE,rx_any_m,Dcm,max_change_per_step)
    call flux_vector(iiBen, ppDcm,ppDcm,rx_any_m)
  end do

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Assign organism-dependent parameters for detritus
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  select case (y)
    case ( iiY1 ) ; cmm  =  ZERO
    case ( iiY4 ) ; cmm  = (cmp+clmp )/ 2.0D+00
    case ( iiY5 ) ; cmm  =  D1m(:)/ 2.0D+00
  end select

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculation of changes due to uptake of detritus in distribution of
  ! state variables (Dx.m is an undetermined source).
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  if (y==iiY2) cmm=FindMidPointExp_vector(D6m,p_from=p_clm(y),d1_to=cmp)
  sx_any=(cmm- abs(D6m)) *(rqt6c-ruQ6c)/(NZERO+Q6c(:))
  call LimitChange_vector(NEGATIVE,sx_any,abs(D6m),max_change_per_step)
  call flux_vector(iiBen,ppD6m,ppD6m,sx_any*sign(DONE,D6m))

  if (y==iiY2) cmm=FindMidPointExp_vector(D7m,p_from=p_clm(y),d1_to=cmp)
  sx_any=(cmm- abs(D7m)) *(rqt6n-ruQ6n)/(NZERO+Q6n(:))
  call LimitChange_vector(NEGATIVE,sx_any,abs(D7m),max_change_per_step)
  call flux_vector(iiBen,ppD7m,ppD7m,sx_any*sign(DONE,D7m))

  if (y==iiY2) cmm=FindMidPointExp_vector(D8m,p_from=p_clm(y),d1_to=cmp)
  sx_any= (cmm- abs(D8m)) *(rqt6p-ruQ6p)/(NZERO+Q6p(:))
  call LimitChange_vector(NEGATIVE,sx_any,abs(D8m),max_change_per_step)
  call flux_vector(iiBen,ppD8m,ppD8m,sx_any*sign(DONE,D8m))

  if (y==iiY2) cmm=FindMidPointExp_vector(D9m,p_from=p_clm(y),d1_to=cmp)
  sx_any=(cmm- abs(D9m)) *(rqt6s)     /(NZERO+Q6s(:))
  call LimitChange_vector(NEGATIVE,sx_any,abs(D9m),max_change_per_step)
  call flux_vector(iiBen,ppD9m,ppD9m,sx_any*sign(DONE,D9m))

  end

! end function
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
