#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: FilterFeeder
!
! DESCRIPTION
!   This process describes the carbon dynamics and associated
!   nutrient dynamics in benthic organism Y3 (suspension feeders)
!   Y3 is handled separately because it also feeds from the water column.
!

!   This file is generated directly from OpenSesame model code, using a code
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
!JM  subroutine FilterFeederDynamics(y,ppyc,ppyn,ppyp)
  subroutine FilterFeederDynamics(y_arg,ppyc_arg,ppyn_arg,ppyp_arg)
!
! !USES:

  ! For the following Benthic-states fluxes are defined: yc, yn, yp, &
  ! G2o, K4n, K1p, D9m
  ! The following Benthic-states are used (NOT in fluxes): D1m
  ! The following Benthic 1-d global boxvars are modified : &
  ! The following Benthic 1-d global boxvars got a value: jPIY3c, jZIY3c, &
  ! jY3RIc, jY3RIn, jY3RIp, jY3RIs
  ! The following Benthic 1-d global boxvars are used: ETW_Ben,PI_Benc,RI_Fc, &
  ! ZI_Fc, PI_Benn, PI_Benp, PI_Bens, ZI_Fn, ZI_Fp, RI_Fn, RI_Fp, RI_Fs
  ! The following global constants are used: RLEN

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN,ZERO,DONE,NZERO,LOGUNIT
  use constants,  ONLY:p_qnUc,MW_C
!JM double #ifdef NOPOINTERS
!  use mem,  ONLY: D2STATE
!#endif
  use mem,  ONLY: D2STATE,Y3c,Yy3c
  use mem, ONLY: ppY3c,ppY3n,ppY3p, ppYy3c, ppQ6c, ppQ6n, ppQ6p, ppG2o, &
    ppYs3c,Ys3c,ppQ1c,ppQ1p,ppQ1n,ppK1p, iiPel,iiBen, ppG3c,iiY3,iiYy3, &
    ppQ6c,ppQ6p,ppQ6n,ppQ6s,ppQun,sediR2_Ben,dry_z,TauBed, &
    ppYy3c,ppYy3n,ppYy3p, flux_vector,max_change_per_step,LocalDelta,   &
    jPIY3c,jP6Y3c,jZEY3c, jZIY3c, jY3RIc,jPLO3c, &
    jspaY3c,jnetY3c,jnetYy3c,jCaCO3Y3c, sediPI_Ben,sediZE_Ben,sediR6_Ben, &
    jRIQIc,jRIQIs,jY3N1p,jY3N4n,jrrY3c,rugY3c, &
    iiC,iiN,iiP,iiL,iiS,&
    O2o_Ben, ZE_Benc,ZE_Benn,ZE_Benp, PI_Benc,PI_Benn,PI_Benp,PI_Bens,PI_Benl,&
    ZI_Fc, ETW_Ben, RI_Fc, ZI_Fn, ZI_Fp, RI_Fn, RI_Fp, RI_Fs, jmY3c,&
    R2_Benn,R2_Benc,R3_Benc,R3_Benn,R3_Benp,puPIY3,puZEY3, &
    puP6Y3, efilP6Y3,efilPART,NO_BOXES_XY, Depth_Ben,  &
    ctfPm2c,ctfZim2c,ctfZem2c,OcDepth,efsatY3
  !apoiter to pelagic
  use mem,only:iiMesoZooPlankton,ppMesoZooPlankton,iiPhytoPlankton,&
    ppPhytoPlankton,iiP1,iiP5,iiP6,iiMicroZooPlankton,ppMicroZooPlankton,ppO2o,&
    ppO3c,ppN1p,ppR2c,ppR2n,ppR6c,ppR6n,ppR6p,ppR6s,ppR3c,ppR1c,ppR1n,ppR1p
! Input: Calculated in PelForcingFoBen-------------------------------------
!   ETW_Ben : temperature
!   PI_Benc : potential food phytoplankton
!   PI_Benn :
!   PI_Benp :
!   PI_Bens :
!   ZI_Fc   : potential food zooplankton
!   ZI_Fn   :
!   ZI_Fp   :
!   RI_Fc   : potential food detritus
!   RI_Fn   :
!   RI_Fp   :
!   RI_Fs   :
!  R3_Benc   : potential food TEP
!  R3_Benn   :
!  R3_Benp   :
!  Depth_Ben :depth(height) lowest layer
! sediPI_Ben   : sedimentation rate of Phyto
! sediR6_Ben   : sedimetation rate of detritus
!  efilP6Y3  : filter limitaion due too high Phaeocystis concentration
!  ctfPm2c   : total amount of Phyto in pelagic per m2 ( used to limit rates dependng on time step)
!  ctfZim2c   :total amount of zooplankton in pelagic per m2
!  ctfRm2c   :total amount of detritus in pelagic per m2
!  cZ2m2c    :total amount of young larvae in pelagic per m2 (used to caluclate quotient young/total filterfeeders)
! Output: ---- Used in BentoPelCoup.F90
!  pyfoodY3  : part of larvae available for other benthic orgranisms
!  jPIY3c   : flux phyto ->filterfeeder (mgC/m2/d)
!  jZIY3c   : flux zoo ->filterfeeder (mgC/m2/d)
!  jY3RIc   : flux  filterfeeder -> detritus (mgC/m2/d)
!  jY3QIc   : flux  filterfeeder -> benthic detritus (mgC/m2/d)
!  jY3QIs   : flux  filterfeeder -> benthic detritus (mgC/m2/d)
! jCaCO3Y3c : shellformation
! jY3O3c    :  CO2 release to pelagic due to respiration
! jRIQIc    :  pseudofaeces production form detritus to benthic detritus
! jRIQIs   :   pseudofaeces production form detritus to benthic detritus
! jY3N1p   :   P flux to pelagic
! jY3N4n   :   N (=urea) flux to pelagic
! diagnostic:
!  puPIY3   : relative availability for food corrected for power distribution near sediment
!  jnetY3c  : nett production of filterfeeders
!  jnetYy3c : nett production of filterfeeders
!  rugYIc  :  rate uptake gros of filterfeeders


  use mem_Param,ONLY:p_pe_R1c,p_pe_R1n,p_pe_R1p,p_peZ_R1c,p_peZ_R1n,p_peZ_R1p, &
    CalcSuspensionFeeders,p_dry_ben
  use LimitRates, ONLY:LimitChange_vector
  use mem_FilterFeeder
  use constants,only:NEGATIVE,POSITIVE
  use botflux,onlY:addbotflux_vector,openbotflux_vector
  use global_interface,only:CorrectConcNearBed_vector,CalcPelMassInM2
  use mem_Phaeo,ONLY:CALC_GRAZING_FILTERFEEDER, CALC_GRAZING_YOUNG_FILTERFEEDER
  use mem_Phyto,only:p_qnRc,p_qpRc,p_qsRc,p_xqn,p_xqp,p_xqs,p_qnR2c



! use mem,  ONLY: Output2d_1,Output2d_2,Output2d_3, Output2d_4,jbotR6c

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following vector functions are used:eTq_vector, eramp_vector, &
  ! MM_vector, PartQ_vector
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use mem_globalfun,   ONLY: eTq_vector, MM_vector, insw_vector,insw

  IMPLICIT NONE
! !INPUT:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!JM  integer,intent(IN)  :: y
!JM  integer,intent(IN) :: ppyc
!JM  integer,intent(IN) :: ppyn
!JM  integer,intent(IN) :: ppyp
  integer,intent(IN)  :: y_arg
  integer,intent(IN) :: ppyc_arg
  integer,intent(IN) :: ppyn_arg
  integer,intent(IN) :: ppyp_arg
!

!
! !AUTHORS
!   P.Ruardij ,        W. Ebenhoeh and C. Kohlmeier
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
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

!
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Set up Local Variable for copy of state var.3object
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN),dimension(NO_BOXES_XY) :: yc
  real(RLEN),dimension(NO_BOXES_XY) :: yn
  real(RLEN),dimension(NO_BOXES_XY) :: yp

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer  :: i,j,k

  integer  :: y
  integer :: ppyc
  integer :: ppyn
  integer :: ppyp

  real(RLEN) :: clu
  real(RLEN) :: corr_max,corr_not
  real(RLEN),dimension(NO_BOXES_XY)  :: R6_corr !correction for power distitustion of R6 ln layer near sediment
  real(RLEN),dimension(NO_BOXES_XY)  :: et    !temperature correction
  real(RLEN),dimension(NO_BOXES_XY)  :: eO    !correction for low oxygen
  real(RLEN),dimension(NO_BOXES_XY)  :: edry  !on tidal flat when it is nearly food uptake is stopped.
  real(RLEN),dimension(NO_BOXES_XY)  :: foodpm2
  real(RLEN),dimension(NO_BOXES_XY)  :: food
  real(RLEN),dimension(NO_BOXES_XY,iiPhytoPlankton):: limit_PIc !calc. food
  real(RLEN),dimension(NO_BOXES_XY,iiPhytoPlankton):: food_PIc !calc. food
  real(RLEN),dimension(NO_BOXES_XY,iiMesoZooPlankton)::food_ZEc ! calc. food
  real(RLEN),dimension(NO_BOXES_XY)  :: food_PT ! calculated phyto food
  real(RLEN),dimension(NO_BOXES_XY)  :: food_ZE ! calculated phyto food
  real(RLEN),dimension(NO_BOXES_XY)  :: food_ZI ! calculated zoo food
  real(RLEN),dimension(NO_BOXES_XY)  :: food_RI ! calculated detritus food
  real(RLEN),dimension(NO_BOXES_XY)  :: food_R3 ! calculated TEP food included in Phaeo colonies
  real(RLEN),dimension(NO_BOXES_XY)  :: R2x_food_c,R2x_food_n
  real(RLEN),dimension(NO_BOXES_XY)  :: suf     ! calculated specific uptake (uptake/food)
  real(RLEN),dimension(NO_BOXES_XY)  :: ruc,qx_any
  real(RLEN),dimension(NO_BOXES_XY)  :: rugc         ! food uptake
  real(RLEN),dimension(NO_BOXES_XY)  :: rugc_corr    ! food uptake
  !specific excretions
  real(RLEN),dimension(NO_BOXES_XY)  :: se_uPIR6c,se_uPIR6n, se_uPIR6p
  real(RLEN),dimension(NO_BOXES_XY)  :: se_uPIR1c,se_uPIR1n, se_uPIR1p
  real(RLEN),dimension(NO_BOXES_XY)  :: se_uR2  ! specific excretion
  real(RLEN),dimension(NO_BOXES_XY)  :: se_uR3  ! specific excretion
  real(RLEN),dimension(NO_BOXES_XY)  :: se_uZIR6c,se_uZIR6n, se_uZIR6p
  real(RLEN),dimension(NO_BOXES_XY)  :: se_uZIR1c,se_uZIR1n, se_uZIR1p
  real(RLEN),dimension(NO_BOXES_XY)  :: se_uZER6c,se_uZER6n, se_uZER6p
  real(RLEN),dimension(NO_BOXES_XY)  :: se_uZER1c,se_uZER1n, se_uZER1p
  real(RLEN),dimension(NO_BOXES_XY)  :: se_uR6  ! specific excretion
  real(RLEN),dimension(NO_BOXES_XY)  :: choice
  real(RLEN),dimension(NO_BOXES_XY)  :: lim_tot_column

  real(RLEN),dimension(NO_BOXES_XY)  :: rrmc    ! maintenance rate respiration C
  real(RLEN),dimension(NO_BOXES_XY)  :: rrac    ! activity rate respiration C
  real(RLEN),dimension(NO_BOXES_XY)  :: sm      ! specific mortality

  ! rate excetion net C N P
  real(RLEN),dimension(NO_BOXES_XY)  :: renc,renn,renp
  ! rate uptake net C N P
  real(RLEN),dimension(NO_BOXES_XY)  :: runc,runn,runp
  ! rate uptake of TEP detritus C,N P, Si
  real(RLEN),dimension(NO_BOXES_XY)  :: ruR3c,ruR3n,ruR3p
  ! rate excretion of TEP  C,N P, Si
  real(RLEN),dimension(NO_BOXES_XY)  :: reR3c,reR3n,reR3p
  ! rate excretion of slow degrading pelagic detritus C,N P, Si
  real(RLEN),dimension(NO_BOXES_XY)  :: reR6c,reR6n,reR6p
  ! rate excretion to total slow degrading pelagic detritus C,N P, Si
  real(RLEN),dimension(NO_BOXES_XY)  :: retR6c,retR6n,retR6p
  real(RLEN),dimension(NO_BOXES_XY)  :: retR1c,retR1n,retR1p
  ! rate excretion to slow degrading benthic detritus C,N P, Si
  real(RLEN),dimension(NO_BOXES_XY)  :: reQ6c,reQ6n,reQ6p
  real(RLEN),dimension(NO_BOXES_XY)  :: reQ1c,reQ1n,reQ1p
  ! rate excretion to total slow degrading benthic detritus C,N P, Si
  real(RLEN),dimension(NO_BOXES_XY)  :: retQ6c,retQ6n,retQ6p,retQ6s
  ! rate excretion phytoplankton C N P
  real(RLEN),dimension(NO_BOXES_XY)  :: rePIR6c,rePIR6n,rePIR6p
  real(RLEN),dimension(NO_BOXES_XY)  :: rePIR1c,rePIR1n,rePIR1p
  ! rate of uptake phytoplankton C N P SI
  real(RLEN),dimension(NO_BOXES_XY)  :: ruPIc,ruPIn,ruPIp, ruPIs
  ! rate of uptake microzooplankton C N P
  real(RLEN),dimension(NO_BOXES_XY)  :: ruZIc,ruZIn,ruZIp
  ! rate of uptake mesozooplankton C N P
  real(RLEN),dimension(NO_BOXES_XY)  :: ruZEc,ruZEn,ruZEp
  ! rate of excretion microzooplankton C N P
  real(RLEN),dimension(NO_BOXES_XY)  :: reZIR6c,reZIR6n,reZIR6p
  real(RLEN),dimension(NO_BOXES_XY)  :: reZIR1c,reZIR1n,reZIR1p
  ! rate of excretion mesozooplankton C N P
  real(RLEN),dimension(NO_BOXES_XY)  :: reZER6c,reZER6n,reZER6p
  real(RLEN),dimension(NO_BOXES_XY)  :: reZER1c,reZER1n,reZER1p

  real(RLEN),dimension(NO_BOXES_XY)  :: reR2c,reR2n   ! rate excretion TEP C
  real(RLEN),dimension(NO_BOXES_XY)  :: ruR2c,ruR2n   ! rate uptake TEP C
  real(RLEN),dimension(NO_BOXES_XY)  :: RTc     ! total available R6c for food uptake
 ! rate uptake detritus  rep. C,B,P,Si
  real(RLEN),dimension(NO_BOXES_XY)  :: ruR6c,ruR6n,ruR6p,ruR6s
  real(RLEN),dimension(NO_BOXES_XY)  :: su      ! specific uptake
  real(RLEN),dimension(NO_BOXES_XY)  :: r     ! --help variables---
  real(RLEN),dimension(NO_BOXES_XY)  :: fluc    ! -- help variables--
  real(RLEN),dimension(NO_BOXES_XY)  :: vum     ! volume filtered for uptake maximum (m3 m-2/mgC)
  real(RLEN),dimension(NO_BOXES_XY)  :: erqu    ! ratio between filtering respiration loss  and max. uptake (-)
  real(RLEN),dimension(NO_BOXES_XY)  :: active  ! if enough foood active > 0.0
  real(RLEN),dimension(NO_BOXES_XY)  :: sunmc   ! rate uptake net maximal C
  real(RLEN),dimension(NO_BOXES_XY)  :: jnetYc  ! rate uptake net maximal C
  real(RLEN),dimension(NO_BOXES_XY) :: x_quality! p6/food
  real(RLEN),dimension(NO_BOXES_XY) :: actual_depth! p6/food
  !resp . any rate,limitation(<1),specific rate,concentration
  real(RLEN),dimension(NO_BOXES_XY) :: rx_any,px_any,sx_any,cx_any,px_limit 
  real(RLEN),dimension(NO_BOXES_XY) :: in_arg !JM added

!write(LOGUNIT,*) 'FilterFeederDynamics',y_arg
!stop

  !JM variables made local
  y=y_arg
  ppyc=ppyc_arg
  ppyn=ppyn_arg
  ppyp=ppyp_arg
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  Copy  state var. object in local var
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  yc = D2STATE(:,ppyc)
  yn = D2STATE(:,ppyn)
  yp = D2STATE(:,ppyp)

   !puP6Y3 stands for the correction for Phaeo Colonies: too largen ones are not eaten.
   limit_PIc =DONE; limit_PIc(:,iiP6)=puP6Y3(:,y)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! response on lowtide on tial flat:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  edry=DONE
  actual_depth=Depth_Ben
  if(p_dry_ben) then
    r =OcDepth(1)*dry_z
    edry= r/(r+0.10)
    actual_depth=actual_depth*dry_z
  endif
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Physiological temperature response
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  et  =   eTq_vector(  ETW_Ben(:),  p_q10(y))

  eO=max(ZERO,DONE-exp(-(O2o_Ben(:)-0.5*p_clO2o(y))/(0.5*p_clO2o(y))))

  !all food offered is  avaible per timestep
  !density of food in layer filtered compared with average density in layer
  !However maxm uptake is limited to p_max/GetDelta().
  corr_max=p_max/LocalDelta
  corr_not=1.0D+80
  !max. amount per timestep which can be eaten from whole pelagic.

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculate total food Cfluxes!
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  clu=p_clu+NZERO
  food  =   ZERO

  ! For phytoplankton:

  food_PT=ZERO;R2x_food_c=ZERO;food_R3=ZERO;R2x_food_n=ZERO
  do i=1,iiPhytoPlankton
     cx_any =  limit_PIc(:,i) * PI_Benc(:,i) *p_puPI(y,i)
     px_any= MM_vector(  cx_any,  clu)
     cx_any =  cx_any  *px_any
     call CorrectConcNearBed_vector(actual_depth(:),sediPI_Ben(:,i),p_height(y),&
                                        corr_not,NO_BOXES_XY,  puPIY3(:,i))
     food_PIc(:,i)=cx_any*puPIY3(:,i)
     food_PT(:)  = food_PT(:)+ food_PIc(:,i)
!write(LOGUNIT,*) 'i,food_PT:',i,food_PT(:)
!write(LOGUNIT,*) 'cx_any',cx_any
!write(LOGUNIT,*) 'puPIY3(:,i)',puPIY3(:,i)
!write(LOGUNIT,*) 'limit_PIc(:,i)',limit_PIc(:,i)
!write(LOGUNIT,*) 'p_puPI(y,i)',p_puPI(y,i)
     select  case  (i)
       case (iiP1,iiP5)
         !take up only R2 is proportion to (resuspended Benthic) diatom
         ! correction for sedimentation is doen in one time
         cx_any =  food_PIc(:,i)*R2_Benc(:)/(NZERO+PI_Benc(:,i))*insw(p_pueR2(y))
         R2x_food_c(:)=R2x_food_c(:)+cx_any
         R2x_food_n(:)=R2x_food_n(:)+cx_any*R2_Benn(:)/(NZERO+R2_Benc(:))
       case (iiP6)
        cx_any =  limit_PIc(:,i) * R3_Benc*  px_any
        food_R3=food_R3+cx_any*puPIY3(:,i)
        food=food+food_R3
     end select
  enddo
  food  =   food  + food_PT(:)

  ! For microzooplankton:
  cx_any=sum(ZI_Fc,1)
  food_ZI  =    cx_any * MM_vector(cx_any,  clu)
  food  =   food+ food_ZI

  ! For meso:
  food_ZE=ZERO
  do i=1,iiMesoZooPlankton
     cx_any =  ZE_Benc(:,i) *p_puZE(y,i)
     px_any =   MM_vector(  cx_any,  clu)
     cx_any =  cx_any  * px_any
     call CorrectConcNearBed_vector(actual_depth(:),sediZE_Ben(:,i),p_height(y), &
                                       corr_not,NO_BOXES_XY,  puZEY3(:,i))
     food_ZEc(:,i)=cx_any*puZEY3(:,i)
     food_ZE(:)  = food_ZE(:)+ food_ZEc(:,i)
  enddo
  food  =   food  + food_ZE(:)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Detritus, (if eaten) first calculate available amount
  ! and add it to the total amount of food
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  !uptake of R2 (EPS) depend on uptqke rate of phytoplankton
  call CorrectConcNearBed_vector(actual_depth(:), sediR2_Ben(:), p_height(y), &
                                             corr_not,NO_BOXES_XY, px_limit)
  R2x_food_c=R2x_food_c*px_limit
  food  =   food+ R2x_food_c

  cx_any=   RI_Fc(:)* MM_vector(  RI_Fc(:),  clu)
  call CorrectConcNearBed_vector(actual_depth(:), sediR6_Ben(:), p_height(y), &
                                             corr_not,NO_BOXES_XY, R6_corr)
  RTc=cx_any*R6_corr
  food_RI=RTc*p_R6(y)
  food  =   food+ food_RI

  px_any=DONE
  if (sw_lim==1.or.sw_lim==3) px_any=min(px_any,efilP6Y3)
  if (sw_lim>1) px_any=min(px_any,efilPart)
  vum=p_vum(y)* px_any

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  !  Food uptake as in zooplankton:
  !  using the modifed Holling response equation which tkae in account
  !  the maximum growth rate and the volume filtered.
  !  Further is assumed that the filterfeeder (nearly) stop filtering
  !  as soon as the costs for filtering  are lower than the  profit
  !  For this we solve the next equation in which r is the unknown:
  !    (left side == profit , right side=costs)
  !    r* p_sum* MM_vector(  vum* food,  r* p_sum)* Y3c(:)*sunmc = p_sra *r
  !  If r > 1.0 : there is enough food to grow
  !  if r < 1.0 : there is balance between costs and profit if
  !  r*erqu*su is larger than the rest respiration.
  !
  !  It is assumed that the detritus sedimentation is defined as a netto
  !  process (p_bursel << P_sediR6). Therefor it assumed that
  ! filterfeeders do noet eat benthic detritus (Q6c).
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  erqu = p_sra(y)/p_sum(y)

  sunmc= (DONE-(p_pueR6(y)*food_RI+p_puePI(y)*food_PT &
              +p_pueR2(y)*R2x_food_c  +p_pueR3(y)*food_R3  &
              +p_pueZE(y)*food_ZE+p_pueZI(y)*food_ZI)/food ) * (DONE-p_pur(y))

  active= edry* (sunmc /erqu - p_sum(y)/(NZERO+ vum *food))

  active=max(1.0D-6,active)
  px_any=min(DONE,active)

  ! Calculate relative uptake
  su= DONE/( DONE/(NZERO+ px_any* vum *food * et *eO ) &
                                   + DONE/(NZERO+ p_sum(y) *et *eO ))
  ! The minimal uptake rate is equal to rest respiration.
  ! With filtering the filterfeeder provide himself also with oxygen.
  rugc= min(su,et*eO*p_sum(y)) *max(ZERO,yc(:))
  ! filtering saturation ( high at low , low at hight food)
  efsatY3(:,y)=min(DONE,su/(et*eO*p_sum(y)))
! fsat=su/(NZERO+ et*eO*edry * vum*food)
  ! Calculate cost of energy for filtering based on realized rate of uptake.
  rrmc = max(edry * eO * p_sra(y)*efsatY3(:,y), p_srs(y))* max(ZERO,yc(:))* et

  !check oxygen and limit growth if Oxygen consumption is too high
  rx_any=(rugc*p_pur(y)+rrmc)/MW_C
  call LimitChange_vector(POSITIVE,rx_any,O2o_Ben,max_change_per_step,px_limit)
  rrmc=px_limit*rrmc
  rugc=px_limit*rugc
  efsatY3(:,y)=px_limit*efsatY3(:,y)

  foodpm2 =food*Depth_Ben

  !diagnostic output:
  rugY3c(:,y)=rugc
  jrrY3c(:,y)=p_srs(y)* max(NZERO,yc(:))* et

  ! Relative growth rate corrected for actual amount of food:

  ! suf= food uptake per unit of food
  suf  =   rugc/ foodpm2

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Execreted part:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! is case p_pueR2 is set to 0, this parameter is used as a switch:
  ! no uptake of R2 by FilterFeeders

  se_uPIR6c  =   suf*min(p_puePI(y),DONE-p_pe_R1c)
  se_uPIR1c  =   suf*p_puePI(y)-se_uPIR6c
  se_uPIR6n  =   suf*min(p_puePI(y),DONE-p_pe_R1n)
  se_uPIR1n  =   suf*p_puePI(y)-se_uPIR6n
  se_uPIR6p  =   suf*min(p_puePI(y),DONE-p_pe_R1p)
  se_uPIR1p  =   suf*p_puePI(y)-se_uPIR6p
  se_uR2  =   suf*p_pueR2(y)*insw(p_pueR2(y))
  se_uR3  =   suf*p_pueR3(y)
!JM  se_uZIR6c  =   suf*min(p_puePI(y),DONE-p_peZ_R1c)
!JM  se_uZIR1c  =   suf*p_puePI(y)-se_uZIR6c
!JM  se_uZIR6n  =   suf*min(p_puePI(y),DONE-p_peZ_R1n)
!JM  se_uZIR1n  =   suf*p_puePI(y)-se_uZIR6n
!JM  se_uZIR6p  =   suf*min(p_puePI(y),DONE-p_peZ_R1p)
!JM  se_uZIR1p  =   suf*p_puePI(y)-se_uZIR6p
  se_uZIR6c  =   suf*min(p_pueZI(y),DONE-p_peZ_R1c)
  se_uZIR1c  =   suf*p_pueZI(y)-se_uZIR6c
  se_uZIR6n  =   suf*min(p_pueZI(y),DONE-p_peZ_R1n)
  se_uZIR1n  =   suf*p_pueZI(y)-se_uZIR6n
  se_uZIR6p  =   suf*min(p_pueZI(y),DONE-p_peZ_R1p)
  se_uZIR1p  =   suf*p_pueZI(y)-se_uZIR6p
!JM  se_uZER6c  =   suf*min(p_puePI(y),DONE-p_peZ_R1c)
!JM  se_uZER1c  =   suf*p_puePI(y)-se_uZIR6c
!JM  se_uZER6n  =   suf*min(p_puePI(y),DONE-p_peZ_R1n)
!JM  se_uZER1n  =   suf*p_puePI(y)-se_uZIR6n
!JM  se_uZER6p  =   suf*min(p_puePI(y),DONE-p_peZ_R1p)
!JM  se_uZER1p  =   suf*p_puePI(y)-se_uZIR6p
  se_uZER6c  =   suf*min(p_pueZE(y),DONE-p_peZ_R1c)
  se_uZER1c  =   suf*p_pueZE(y)-se_uZER6c
  se_uZER6n  =   suf*min(p_pueZE(y),DONE-p_peZ_R1n)
  se_uZER1n  =   suf*p_pueZE(y)-se_uZER6n
  se_uZER6p  =   suf*min(p_pueZE(y),DONE-p_peZ_R1p)
  se_uZER1p  =   suf*p_pueZE(y)-se_uZER6p
  se_uR6  =   suf*p_pueR6(y)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! CALCULATION OF UPTAKE RATE:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! calculation of quotum N/C P/C in exretion product slow degrading detritus
  ! (R6)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!write(LOGUNIT,*) 'somewhere halfway'
!stop

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Pelagic Phytoplankton:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  ruPIc  =   ZERO ; rePIR1c=ZERO;  rePIR6c=ZERO;  ruR3c  =   ZERO ; reR3c=ZERO
  ruPIn  =   ZERO ; rePIR1n=ZERO;  rePIR6n=ZERO;  ruR3n  =   ZERO ; reR3n=ZERO
  ruPIp  =   ZERO ; rePIR1p=ZERO;  rePIR6p=ZERO;  ruR3p  =   ZERO ; reR3p=ZERO
  ruPIs  =   ZERO ;  ruR2c  =   ZERO;  reR2c=ZERO;ruR2n  =   ZERO;  reR2n=ZERO

  ! calculate max uptake  of food/m2. Compare uptake pressure
  ! with food in whole water column and limit if necessary
  fluc=food_PT *suf *Depth_Ben ! pot. uptaken food  mgC/m2
!write(LOGUNIT,*) 'fluc1',fluc
!write(LOGUNIT,*) 'food_PT:',food_PT(:)
  call LimitChange_vector(POSITIVE,fluc,ctfPm2c,max_change_per_step, &
                                                         lim_tot_column)
!write(LOGUNIT,*) 'fluc2',fluc
!write(LOGUNIT,*) 'ctfPm2c',ctfPm2c
!write(LOGUNIT,*) 'max_change_per_step',max_change_per_step


  k=CALC_GRAZING_FILTERFEEDER;if (y.eq.iiYy3)k=CALC_GRAZING_YOUNG_FILTERFEEDER
  do i=1,iiPhytoPlankton
    ! calculate rel max uptake  of /m3/d. Check uptake pressure
    ! in lowest layer and limit if necessary
    fluc = suf*food_PIc(:,i)
    call LimitChange_vector(POSITIVE,fluc,PI_Benc(:,i),corr_max,px_limit)
    choice=food_PIc(:,i)* Depth_Ben/(NZERO + PI_Benc(:,i))*min(lim_tot_column,px_limit)
    ! all uptake below are uptakes per square m2  dimension choice :m
    ruc=suf*choice*PI_Benc(:,i)
    jPIY3c(:,i) = jPIY3c(:,i) + ruc
    j= ppPhytoPlankton(i,iiC)
    call addbotflux_vector(POSITIVE,iiPel,j,iiBen,ppyc,ruc)
!if (i==2) then
!if (minval(ruc).le.0.0) then
!write(LOGUNIT,*) 'FilterFeeder: ruc',ruc
!write(LOGUNIT,*) 'i',i
!write(LOGUNIT,*) 'PI_Benc:',PI_Benc(:,i)
!write(LOGUNIT,*) 'suf:',suf
!write(LOGUNIT,*) 'choice:',choice
!write(LOGUNIT,*) 'food_PIc:',food_PIc(:,i)
!write(LOGUNIT,*) 'lim_tot_column,px_limit,corr_max:',lim_tot_column,px_limit,corr_max
!write(LOGUNIT,*) 'fluc:',fluc
!write(LOGUNIT,*) 'Depth_Ben:',Depth_Ben
!endif
!endif
    ruPIc  = ruPIc  +  ruc
    rePIR6c  = rePIR6c  +  PI_Benc(:,i)* se_uPIR6c* choice
    rePIR1c  = rePIR1c  +  PI_Benc(:,i)* se_uPIR1c* choice
    !N fluxes
!   t=suf*choice*PI_Benn(:,i)
    j= ppPhytoPlankton(i,iiN)
    qx_any=PI_Benn(:,i)/(NZERO+PI_Benc(:,i))
    qx_any=min(qx_any,p_qnRc(i)*p_xqn(i))
    rx_any=suf*qx_any*choice*PI_Benc(:,i)
    call addbotflux_vector(POSITIVE,iiPel,j,iiBen,ppyn,rx_any)
    ruPIn  = ruPIn  +  rx_any
!   rePIn  = rePIn  +  rec*qx_any
    rePIR6n  = rePIR6n  +  PI_Benn(:,i)* se_uPIR6n* choice
    rePIR1n  = rePIR1n  +  PI_Benn(:,i)* se_uPIR1n* choice
    !P fluxes
!   t=suf*choice*PI_Benp(:,i)
    j= ppPhytoPlankton(i,iiP)
    qx_any=PI_Benp(:,i)/(NZERO+PI_Benc(:,i))
    qx_any=min(qx_any,p_qpRc(i)*p_xqp(i))
    rx_any=suf*qx_any*choice*PI_Benc(:,i)
    call addbotflux_vector(POSITIVE,iiPel,j,iiBen,ppyp,rx_any)
    ruPIp  = ruPIp  + rx_any
!   rePIp  = rePIp  + rec*qx_any
    rePIR6p  = rePIR6p  +  PI_Benp(:,i)* se_uPIR6p* choice
    rePIR1p  = rePIR1p  +  PI_Benp(:,i)* se_uPIR1p* choice
    !Chl fluxes
    j= ppPhytoPlankton(i,iiL) ;rx_any=PI_Benl(:,i)* suf* choice
    call openbotflux_vector(NEGATIVE,iiPel,j,-rx_any)
    !Si fluxes
    j= ppPhytoPlankton(i,iiS) ; if (j>0) then
      qx_any=PI_Bens(:,i)/(NZERO+PI_Benc(:,i))
      qx_any=min(qx_any,p_qsRc(i)*p_xqs(i))
      ruPIs=suf*qx_any*choice*PI_Benc(:,i)
      call addbotflux_vector(POSITIVE,iiPel,j,iiBen,ppQ6s,ruPIs)
    endif

    if (i==iiP6) then
       !if more than one  filterfeedergroup members eat Phaeocystis and eat
       !from different size
       !classes we need to know how much eachmember
       ! eat in order to correct for size changes in the PhaeoColonies
       ! See alsp in Settling.
       ruc=PI_Benc(:,i)* suf* choice
       jP6Y3c(:,y)= jP6Y3c(:,y)+ruc
       call PhaeocystisCalc_1l(k,i, r,ruc,DONE)

       !Phaeo.uptake: carbon outside cells but in colony
        ruR3c= R3_Benc *suf * choice
        call addbotflux_vector(POSITIVE,iiPel,ppR3c,iiBen,ppyc,ruR3c)
        ruR3n= R3_Benn *suf * choice
        j= ppPhytoPlankton(i,iiN)
        call addbotflux_vector(POSITIVE,iiPel,j,iiBen,ppyn,ruR3n)
        ruR3p= R3_Benp *suf * choice
        j= ppPhytoPlankton(i,iiP)
        call addbotflux_vector(POSITIVE,iiPel,j,iiBen,ppyp,ruR3p)

        reR3c= R3_Benc* se_uR3* choice
        reR3n= R3_Benn* se_uR3* choice
        reR3p= R3_Benp* se_uR3* choice
    endif
  enddo

  ! new food originating for phytoplankton for Y3 is added here to Y3
  ! losses to phytoplankton are set in BenPelCoup using the total
  ! total phytoplankton uptake jPIY3c

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Pelagic MicroZooplankton:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  ! calculate max uptake  of food/m3. Check uptake pressure
  ! in lowest layer and limit if necessary
  fluc=food_ZI *suf ;cx_any=sum(ZI_Fc,1)
  call LimitChange_vector(POSITIVE,fluc,cx_any,corr_max,px_limit)
  ! calculate max uptake  of food/m2. Compare uptake pressure
  ! with food in whole water column and limit if necessary
  fluc=fluc *Depth_Ben
  call LimitChange_vector(POSITIVE,fluc,ctfZim2c,max_change_per_step, &
                                                          lim_tot_column)

  !Calculate relative uptaken food  (m)
  choice  =   min(lim_tot_column,px_limit) *food_ZI*Depth_Ben/(NZERO+cx_any)

  ! mgC/m2/d = mgC/m3  */d * m
  ruZIc=ZERO;ruZIn=ZERO;ruZIp=ZERO
  reZIR6c=ZERO;reZIR6n=ZERO;reZIR6p=ZERO
  reZIR1c=ZERO;reZIR1n=ZERO;reZIR1p=ZERO
  do i=1,iiMicroZooPlankton
    rx_any=ZI_Fc(:,i)* suf* choice; ruZIc=ruZIc+rx_any
    j= ppMicroZooPlankton(i,iiC)
    call addbotflux_vector(POSITIVE,iiPel,j,iiBen,ppyc,rx_any)
    rx_any=ZI_Fn(:,i)* suf* choice; ruZIn=ruZIn+rx_any
    j= ppMicroZooPlankton(i,iiN)
    if (j>0) call addbotflux_vector(POSITIVE,iiPel,j,iiBen,ppyn,rx_any)
    if (j<1)call flux_vector(iiBen,ppyn,ppyn,rx_any)
    rx_any=ZI_Fp(:,i)* suf* choice; ruZIp=ruZIp+rx_any
    j= ppMicroZooPlankton(i,iiP)
    if (j>0) call addbotflux_vector(POSITIVE,iiPel,j,iiBen,ppyp,rx_any)
    if (j<1)call flux_vector(iiBen,ppyp,ppyp,rx_any)

    reZIR6c  =reZIR6c +   ZI_Fc(:,i)* se_uZIR6c* choice
    reZIR1c  =reZIR1c +   ZI_Fc(:,i)* se_uZIR1c* choice
    reZIR6n  =reZIR6n +   ZI_Fn(:,i)* se_uZIR6n* choice
    reZIR1n  =reZIR1n +   ZI_Fn(:,i)* se_uZIR1n* choice
    reZIR6p  =reZIR6p +   ZI_Fp(:,i)* se_uZIR6p* choice
    reZIR1p  =reZIR1p +   ZI_Fp(:,i)* se_uZIR1p* choice
  enddo

  jZIY3c(:)  = jZIY3c(:) +    ruZIc

  ! new food originating fro microzooplankton for Y3 is added here to Y3
  ! losses to microzooplankton are set in BenPelCoup using the total
  ! total microzooplankton uptake jZIY3c


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Pelagic MesoZooplankton:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  ruZEc  =   ZERO ; reZER6c=ZERO; reZER1c=ZERO
  ruZEn  =   ZERO ; reZER6n=ZERO; reZER1n=ZERO
  ruZEp  =   ZERO ; reZER6p=ZERO; reZER1p=ZERO

  ! calculate max uptake  of food/m2. Compare uptake pressure
  ! with food in whole water column and limit if necessary
    fluc=food_ZE *suf *Depth_Ben ! pot. uptaken food  mgC/m2
    call LimitChange_vector(POSITIVE,fluc,ctfZem2c,max_change_per_step, &
                                                              lim_tot_column)

    do i=1,iiMesoZooPlankton
      ! calculate rel max uptake  of /m3/d. Check uptake pressure
      ! in lowest layer and limit if necessary
      fluc = suf*food_ZEc(:,i)
      call LimitChange_vector(POSITIVE,fluc, ZE_Benc(:,i),corr_max,px_limit)
      choice=food_ZEc(:,i)* Depth_Ben/(NZERO + ZE_Benc(:,i))* &
                                                    min(lim_tot_column,px_limit)
      ! all uptake below are uptakes per square m2  dimension choice :m
      jZEY3c(:,i) = jZEY3c(:,i) + ZE_Benc(:,i)* suf* choice
      rx_any=  ZE_Benc(:,i)* suf* choice;ruZEc  = ruZEc  +   rx_any
      j= ppMesoZooPlankton(i,iiC)
      call addbotflux_vector(POSITIVE,iiPel,j,iiBen,ppyc,rx_any)
      rx_any=  ZE_Benn(:,i)* suf* choice;ruZEn  = ruZEn  +   rx_any
      j= ppMesoZooPlankton(i,iiN)
      if (j>0) call addbotflux_vector(POSITIVE,iiPel,j,iiBen,ppyn,rx_any)
      if (j<1)call flux_vector(iiBen,ppyn,ppyn,rx_any)
      rx_any=  ZE_Benp(:,i)* suf* choice;ruZEp  = ruZEp  +   rx_any
      j= ppMesoZooPlankton(i,iiP)
      if (j>0) call addbotflux_vector(POSITIVE,iiPel,j,iiBen,ppyp,rx_any)
      if (j<1)call flux_vector(iiBen,ppyp,ppyp,rx_any)

      reZER6c  = reZER6c  +   ZE_Benc(:,i)* se_uZER6c* choice
      reZER1c  = reZER1c  +   ZE_Benc(:,i)* se_uZER1c* choice
      reZER6n  = reZER6n  +   ZE_Benn(:,i)* se_uZER6n* choice
      reZER1n  = reZER1n  +   ZE_Benn(:,i)* se_uZER1n* choice
      reZER6p  = reZER6p  +   ZE_Benp(:,i)* se_uZER6p* choice
      reZER1p  = reZER1p  +   ZE_Benp(:,i)* se_uZER1p* choice
    enddo

  ! new food originating for mesozooplankton for Y3 is added here to Y3
  ! losses to mesozooplankton are set in BenPelCoup using the total
  ! total mesozooplankton uptake jZEY3c
!write(LOGUNIT,*) 'quite a bit further'
!stop

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Pelagic Detritus
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  call CorrectConcNearBed_vector(Depth_Ben(:),sediR2_Ben(:),p_height(y),&
                                        p_max,NO_BOXES_XY,px_limit)
  ruR2c= R2x_food_c *suf *Depth_Ben*px_limit
  call addbotflux_vector(POSITIVE,iiPel,ppR2c,iiBen,ppyc,ruR2c)
  reR2c= R2x_food_c *se_uR2 *Depth_Ben*px_limit

!write(LOGUNIT,*) 'here1'
!stop

!JM
if (p_qnR2c > NZERO) then
  ruR2n= R2x_food_n *suf *Depth_Ben*px_limit
  call addbotflux_vector(POSITIVE,iiPel,ppR2n,iiBen,ppyn,ruR2n)
  reR2n= R2x_food_n *se_uR2 *Depth_Ben*px_limit
endif

!write(LOGUNIT,*) 'here2'
!stop

  fluc=RI_Fc *suf *Depth_Ben
  !limit uptake to maximum of p_max per step (see at calulation of p_max)
!write(LOGUNIT,*) 'here3'
!stop
  cx_any=CalcPelMassInM2(ppR6c)
!write(LOGUNIT,*) 'here4'
!stop
  call LimitChange_vector(POSITIVE,fluc,cx_any,max_change_per_step, &
                                                        lim_tot_column)
  choice  =   food_RI * Depth_Ben/(NZERO+RI_Fc(:))*lim_tot_column

!write(LOGUNIT,*) 'here'
!stop

  ruR6c      =   RI_Fc(:)* suf* choice
  ruR6n      =   RI_Fn(:)* suf* choice
  ruR6p      =   RI_Fp(:)* suf* choice
  ruR6s      =   RI_Fs(:)* suf* choice
  call addbotflux_vector(POSITIVE,iiPel,ppR6c,iiBen,ppyc,ruR6c)
  call addbotflux_vector(POSITIVE,iiPel,ppR6n,iiBen,ppyn,ruR6n)
  call addbotflux_vector(POSITIVE,iiPel,ppR6p,iiBen,ppyp,ruR6p)
  call addbotflux_vector(POSITIVE,iiPel,ppR6s,iiBen,ppQ6s,ruR6s)

  reR6c  =   RI_Fc(:)* se_uR6* choice
  reR6n  =   RI_Fn(:)* se_uR6* choice
  reR6p  =   RI_Fp(:)* se_uR6* choice

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Book keeping
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  runc  =   ruPIc+ ruR2c + ruR3c + ruZIc+ ruZEc+ ruR6c
  runn  =   ruPIn+ ruR2n + ruR3n + ruZIn+ ruZEn+ ruR6n
  runp  =   ruPIp+         ruR3p + ruZIp+ ruZEp+ ruR6p

  retR6c  =   rePIR6c+ reR2c +reR3c + reZIR6c+ reZER6c+ reR6c
  retR6n  =   rePIR6n+ reR2n +reR3n + reZIR6n+ reZER6n+ reR6n
  retR6p  =   rePIR6p+        reR3p + reZIR6p+ reZER6p+ reR6p

  retR1c  =   rePIR1c+ reZIR1c+ reZER1c
  retR1n  =   rePIR1n+ reZIR1n+ reZER1n
  retR1p  =   rePIR1p+ reZIR1p+ reZER1p

  retQ6c  =   ZERO
  retQ6n  =   ZERO
  retQ6p  =   ZERO
  retQ6s  =   ruPIs + ruR6s

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! net growth rate ( only for output)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!write(LOGUNIT,*) 'net growth rate'
!stop

  rugc_corr=ruPIc+ ruR2c + ruR3c +  ruZIc+ ruZEc+ ruR6c
  rrac= p_pur(y)*( rugc_corr- retR1c-retR6c- retQ6c)
  jnetYc=max(ZERO,ruPIc+ ruR2c + ruR3c +  ruZIc+ ruZEc+ ruR6c &
          -rrac -retQ6c-retR1c-retR6c)

  ! in case of a negative value of one of the following values there is a &
  ! situation
  ! of startvation and very low biomass values. Check on quota in the food is &
  ! out of order

  sx_any=(rrmc+rrac)/(NZERO+yc)
  runc  =   max(  ZERO,  runc -retR1c-retR6c-retQ6c-rrmc-rrac)
  runn  =   max(  ZERO,  runn -retR1n-retR6n-retQ6n+sx_any*yn)
  runp  =   max(  ZERO,  runp -retR1p-retR6p-retQ6p+sx_any*yp)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculation of nutrient release and correction of C:N:P
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  ! Carbon which will be excreted/respired when internal quota (N/C) are
  ! below optimum
  renc  =   max(ZERO,(-runn/( NZERO+ runc)/p_qnc+DONE),   &
                    (-runp/( NZERO+ runc)/p_qpc+DONE))* runc

  ! take in consideration that a (small) part of runC is used to make urea
  runc=runc-renc
  renn=max(ZERO,runn-p_qnc*runc)/(DONE-p_qnc/p_qnUc)
  ! Correct uptake of N for loss to ure production
  runc=runc-renn/p_qnUc
  renp=max(ZERO,runp-p_qpc*runc)

  ! redistribute renc ovr flux to pelagic and benthic
  rx_any=retQ6c;r=retR6c
  retQ6c= retQ6c +renc * rx_any/(NZERO + rx_any+ r)
  retR6c= retR6c +renc * r/(NZERO + rx_any+ r)

  renn = max( ZERO, renn+ yn(:) -DONE * p_qnc* yc(:))
  renp = max( ZERO, renp+ yp(:) -DONE * p_qpc* yc(:))


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculation of respiration:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  r  =  rrmc+ max(ZERO,rrac)

  call flux_vector(               iiBen, ppyc,ppG3c,  r*(DONE-p_pePel(y)) )
  call addbotflux_vector(POSITIVE,iiBen,ppyc,iiPel,ppO3c,r*p_pePel(y))
  jPLO3c=jPLO3c+r*p_pePel(y)
  jrrY3c(:,y)=r

  call flux_vector(iiBen, ppG2o,ppG2o,-( r/ MW_C))
  call addbotflux_vector(POSITIVE,iiPel,ppO2o,iiBen,ppG2o, &
                                               r/MW_C *(DONE-p_pePel(y)))

  ! Here only sink flux defined for Y3
  ! Source fluxes are here only calculated in jY3N4n and JY3N1p
  ! Flux correction for N4n and N1p in BenPelCoup

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! fluxes
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rx_any=(renn+retR1n)*p_pePel(y)
  call addbotflux_vector(POSITIVE,iiBen,ppyn,iiPel,ppR1n,rx_any)
  rx_any=(renn/p_qnUc+retR1c)*p_pePel(y)
  call addbotflux_vector(POSITIVE,iiBen,ppyc,iiPel,ppR1c, rx_any)
  call addbotflux_vector(POSITIVE,iiBen,ppyp,iiPel,ppN1p,renp*p_pePel(y))
  call addbotflux_vector(POSITIVE,iiBen,ppyp,iiPel,ppR1p,retR1p*p_pePel(y))

  jY3N4n(:)=jY3N4n(:)+renn *p_pePel(y)
  jY3N1p(:)=jY3N1p(:)+renp *p_pePel(y)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculation of mortality
  !      first order: mortality; "age"
  !      second order: mortality; "density"
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- =-=-=-=-=-=-=-=-=-=-=-=-=-

  sm  =   p_sd(y)* et  +p_sd2(y) * max(ZERO,Y3c+Yy3c)
  if ( y==iiY3) sm=sm+p_sxresus(y)*TauBed/(p_xtauC(y)+TauBed)&
                       *insw_vector(TauBed-p_xtauC(y))

  rx_any  =   max(ZERO,yc(:))* sm
  reQ1c  =   rx_any*p_peZ_R1c
  reQ6c  =   rx_any-reQ1c
  rx_any  =   max(ZERO,yn(:)) *sm
  reQ1n  =   rx_any*p_peZ_R1n
  reQ6n  =   rx_any-reQ1n
  rx_any  =   max(ZERO,yp(:)) *sm
  reQ1p  =   rx_any*p_peZ_R1p
  reQ6p  =   rx_any-reQ1p

  jmY3c(:,y)=jmY3c(:,y)+reQ6c
  retQ6c  =   retQ6c+ reQ6c
  retQ6n  =   retQ6n+ reQ6n
  retQ6p  =   retQ6p+ reQ6p

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  !Full flux defintion
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  call flux_vector( iiBen, ppyn,ppQ1n,reQ1n+retR1n* (DONE-p_pePel(y)))
  call flux_vector( iiBen, ppyn,ppQun,renn* (DONE-p_pePel(y)))
  call flux_vector( iiBen, ppyc,ppQ1c, &
                              reQ1c+(retR1c+renn/p_qnUc)* (DONE-p_pePel(y)))
  call flux_vector( iiBen, ppyp,ppQ1p,reQ1p+retR1p * (DONE-p_pePel(y)))
  call flux_vector( iiBen, ppyp,ppK1p,renp * (DONE-p_pePel(y)))

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculate total flux from Suspension feeders to Q6:
  ! Full flux:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  call flux_vector( iiBen, ppyc,ppQ6c, retQ6c+ retR6c* (DONE-p_pR6Pel(y)))
  call flux_vector( iiBen, ppyn,ppQ6n, retQ6n+ retR6n* (DONE-p_pR6Pel(y)))
  call flux_vector( iiBen, ppyp,ppQ6p, retQ6p+ retR6p* (DONE-p_pR6Pel(y)))

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculate total flux from Suspension feeders to R6:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  call addbotflux_vector(POSITIVE,iiBen,ppyc,iiPel,ppR6c,retR6c*p_pR6Pel(y))
  call addbotflux_vector(POSITIVE,iiBen,ppyn,iiPel,ppR6n,retR6n*p_pR6Pel(y))
  call addbotflux_vector(POSITIVE,iiBen,ppyp,iiPel,ppR6p,retR6p*p_pR6Pel(y))
  call addbotflux_vector(POSITIVE,iiBen,ppQ6s,iiPel,ppR6s,retQ6s*p_pR6Pel(y))

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculate NET flux from R6 to Suspension feeders :
  ! (can be negative!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  ! The ruR6s which is uptaken is directly relased back to R6: net food flux &
  ! from Y3 to/from R6 is 0
  jY3RIc(:)  =jY3RIc(:)  +   retR6c


!write(LOGUNIT,*) 'pseudo faeces'
!stop


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! pseudo faeces production
  ! This production lead only to a flux to the sediment!
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  !uptake of detritus for pseudofaeces production per m3
  ! Check flux in lowest layer
  ! Check flux with mass of R6 in the whole water-column
  call CorrectConcNearBed_vector(Depth_Ben(:),sediR6_Ben(:),p_height(y),&
                                        corr_not,NO_BOXES_XY,px_limit)
  !flux mg detritus /m3
  fluc =   px_limit*(DONE- p_pR6Pel(y))  &
         * edry * et * eO *  vum * efsatY3(:,y)* yc(:) *RTc
  call LimitChange_vector(POSITIVE,fluc,RTc,corr_max,px_limit)
  !flux mg detritus /m2
  fluc= fluc*Depth_Ben
  sx_any=fluc/(NZERO + RI_Fc(:))
!write(LOGUNIT,*) 'pf1'
!stop
  
  cx_any=CalcPelMassInM2(ppR6c)
  call LimitChange_vector(POSITIVE,fluc,cx_any,corr_max,lim_tot_column)
  rx_any = sx_any * min(lim_tot_column,px_limit)*RI_Fc
  jRIQIc(:)=jRIQIc(:)+max(ZERO,r * RI_Fc(:) -ruR6c)

!write(LOGUNIT,*) 'pf2'
!stop

  call addbotflux_vector(POSITIVE,iiPel,ppR6c,iiBen,ppQ6c, &
                                           +max(ZERO,rx_any-ruR6c))
  cx_any=CalcPelMassInM2(ppR6n)
  rx_any= sx_any*RI_Fn
  call LimitChange_vector(POSITIVE,rx_any,cx_any,corr_max,lim_tot_column)
  rx_any = rx_any * min(lim_tot_column,px_limit)
  call addbotflux_vector(POSITIVE,iiPel,ppR6n,iiBen,ppQ6n, &
                                           +max(ZERO,rx_any-ruR6n))
!write(LOGUNIT,*) 'pf3'
!stop

  cx_any=CalcPelMassInM2(ppR6p)
  rx_any= sx_any*RI_Fp
  call LimitChange_vector(POSITIVE,rx_any,cx_any,corr_max,lim_tot_column)
!write(LOGUNIT,*) 'pf4'
!stop
  rx_any = rx_any * min(lim_tot_column,px_limit)
!write(LOGUNIT,*) 'pf4b'
!stop
  call addbotflux_vector(POSITIVE,iiPel,ppR6p,iiBen,ppQ6p, &
                                           +max(ZERO,rx_any -ruR6p))
!write(LOGUNIT,*) 'pf5'
!stop

  !uptake of detritus for pseudofaeces production per m3
  call CorrectConcNearBed_vector(Depth_Ben(:),sediR6_Ben(:),p_height(y), &
                                        corr_not,NO_BOXES_XY,px_limit)
  fluc =   px_limit*(DONE- p_pR6Pels(y))  &
       * edry * et * eO * vum * efsatY3(:,y)* yc(:) *RTc
  call LimitChange_vector(POSITIVE,fluc,RTc,corr_max,px_limit)
  sx_any=fluc/(NZERO + RI_Fc(:))
!write(LOGUNIT,*) 'pf55'
!stop

  cx_any=CalcPelMassInM2(ppR6s)
  rx_any= sx_any*RI_Fs
  call LimitChange_vector(POSITIVE,rx_any,cx_any,corr_max,lim_tot_column)
!write(LOGUNIT,*) 'pf58'
!stop
  rx_any = rx_any * min(lim_tot_column,px_limit)
  jRIQIs(:)=jRIQIs(:)+max(ZERO,rx_any-ruR6s)
!write(LOGUNIT,*) 'pf59'
!stop
  call addbotflux_vector(POSITIVE,iiPel,ppR6s,iiBen,ppQ6s, &
                                           +max(ZERO,rx_any -ruR6s))

!write(LOGUNIT,*) 'pf6'
!stop
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! larvae spawning is setupped in such a way that it start when food
  ! saturation is reached and filtering-rate can be lowered.
  ! The variable fsat is the controlling factor which trigger the spawning
  ! as soon as fsat is lower than p_fsat_y (default value: 0.7)
  ! The spawning is modelled as an event by dividing the rate by the timestep.
  ! The spawning is modelled in pragmatic way: new spawning events are limited
  ! by the quotient young/total. In this way only 2 per year a spwawing-event .
  ! take place(in spring and sometimes in autumn)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!write(LOGUNIT,*) 'larvae spawning'
!stop

  if (CalcSuspensionFeeders(iiYy3)) then
    select case (y)
      case (iiY3)
        jspaY3c=NZERO
!       cx_any=food_PIc(iiP6,:)*(DONE-efilP6Y3)
!       x_quality=(food -food_PIc(iiP6,:))/(food+cx_any)
        if (p_xspawning>ZERO) then
          !limit spawning if already a high biomass of young ones are present.
          !This is a way to avoid modelling explicit the modelling the reserve
          ! biomass which is excreted as young ones.
          sx_any=max(ZERO,p_xspawning-(Ys3c)/(NZERO+yc))
          ! biomass increase of adult filterfeeders must at least p_lxspawning!
          !only swpawning,when
          ! active: filterfeeder filters because there is enough food present
          ! there is food saturation
          ! quality food is high (low amount of pHaeocsysts)
          ! above p_clxTemp (10 C)
          in_arg=efsatY3(:,y)-p_fsat_y  !JM pass this way
          sx_any=sx_any*insw_vector(active-DONE) &
!JM            *insw_vector(efsatY3(:,y)-p_fsat_y) &
            *insw_vector(in_arg) &
!           *insw_vector(x_quality-0.50) &
            *insw_vector(ETW_Ben-p_clxTemp)
          sx_any=sx_any*insw_vector(sx_any-p_lxspawning) *&
                                                     insw_vector(10.0-TauBed)
          jspaY3c=max(NZERO,sx_any*yc)/LocalDelta
        endif
        call LimitChange_vector(POSITIVE,jspaY3c,yc,max_change_per_step)
        call flux_vector(iiBen,ppYs3c,ppYs3c, &
                          +jspaY3c-Ys3c*p_smYy3c/et*insw_vector(Ys3c) )
        call flux_vector( iiBen, ppY3c,ppYy3c,jspaY3c )
        call flux_vector( iiBen, ppY3n,ppYy3n,jspaY3c *p_qnc)
        call flux_vector( iiBen, ppY3p,ppYy3p,jspaY3c *p_qpc)
        jnetY3c=jnetYc
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! shell formation
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      case (iiYy3)
        jnetYy3c= jnetYc
        cx_any=DONE
        !artifical way to move young ones to old .This to overcome that young
        !ones are still prresent at start of growing season.
        sx_any=-50.0*min(ZERO,jnetYc-rrmc)/(NZERO+yc)
        call LimitChange_vector(POSITIVE,sx_any,cx_any,corr_max)
        call flux_vector(iiBen,ppYs3c,ppYs3c,-Ys3c*sx_any/et*insw_vector(Ys3c) ) 
        sx_any=p_smYy3c+sx_any
        call flux_vector(iiBen, ppYy3c,ppY3c, yc* sx_any/et*insw_vector(yc) )
        call flux_vector(iiBen, ppYy3n,ppY3n, yn* sx_any/et*insw_vector(yn) )
        call flux_vector(iiBen, ppYy3p,ppY3p, yp* sx_any/et*insw_vector(yp) )
        jCaCO3Y3c=jCaCO3Y3c+p_qCaCO3Y3c * max(ZERO,jnetY3c-jnetYy3c)
    end select
  else
    jnetY3c=jnetYc
  endif
!write(LOGUNIT,*) 'end FilterFeederDynamics'
!stop

  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
