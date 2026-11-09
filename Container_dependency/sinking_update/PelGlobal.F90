#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: PelGlobal
!
! DESCRIPTION
!   !
!

!   This file is generated directly from OpenSesame model code, using a code
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine PelGlobalDynamics(mode)
!
! !USES:
  ! The following Pelagic-states are used (NOT in fluxes): R6p, R6c, R6n, R6s, &
  ! P1s, P1c, B1p, B1c, B1n
  ! The following box states are used (NOT in fluxes): &
  ! MicroZooPlankton, MesoZooPlankton, PhytoPlankton
  ! The following Pelagic 1-d global boxvars got a value: flP1R6s, flPTN6r, &
  ! qpR6c, qnR6c, qsR6c, qpB1c, qnB1c, sediR6
  ! The following Pelagic 2-d global boxvars got a value: qp_mz, qn_mz, qpZc, &
  ! qnZc, qpPc, qnPc, qlPc, qsPc, sediPI
  ! The following groupmember vars are used: iiMicroZooPlankton, &
  ! iiMesoZooPlankton, iiPhytoPlankton

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN,LOGUNIT,ZERO,NZERO,DONE
#ifdef NOPOINTERS
  use mem,  ONLY: D3STATE
#else
  use mem, ONLY: N1p,N3n, R1c,R1n,R1p, R2c,R6p, R6c, R6n, R6s,R9x, RZc,O2o, &
    B1p,B1c, B1n,P5c, MicroZooPlankton, MesoZooPlankton, PhytoPlankton,Pcc
#endif
  use mem, ONLY: ppR6c, ppB1c,ppZ6c,ppP3c, &
    ppMicroZooPlankton, ppMesoZooPlankton, ppPhytoPlankton,iiPhytoPlankton, &
    iiMicroZooPlankton,iiMesoZooPlankton,iiP,iiC,iiN,iiL,iiS,iiP2,iiP5,iiP6, &
    ppB1c,B1c,ETW,ESW,iiP1, &
    flN3N4n,jnetPIc,jnetB1c,jnetMeZc, jnetMiZc,jPelFishInput, sediMeZ,sediMiZ,&
    qpR6c,qnR6c,qsR6c, qpB1c,qnB1c, qp_mz,qn_mz,qpZc,qnZc, qpPc,qnPc,qlPc,qsPc,&
    limnuti,PTi,flPTN6r, jZIR6n,jZIR6p,jZIDIn,jZIDIp,rnetPTc,jmY3c,jPLO3c,&
    NO_BOXES_XY,NO_BOXES,NO_BOXES_Z,Depth,OCDepth,PelBoxAbove,BMLd, &
    jnetY3c,jnetYy3c,rml,jPTZTn,jPTZTp,flB1N4n,flPIN4n,flN3N4n,PI_dw,&
    flPIR6n,flPIR1n,flPIR6p,flPIR1p,flPIR6s,flR3R2c,flBaZTc,flR1O3c,flR1N4n, &
    Nun,Nup,Chla, O2o_vr,N1p_vr,N3n_vr,Chla_vr,R6c_vr,jrrPTc,jrrMic,jrrMec, &
    iiProduction,p_xfree_R2,iiR2,qR2P1, &
    fl_xgrazing_PIc,POC,PON,POP,xSizeMA_m,xSizeMA_d

  use mem,ONLY:sediB1,sediR2, sediR6,sediR6s, sediR9, sediRZ,sediPI,&
     CoupledtoBDc
  use mem_PelGlobal
  use BFM_ERROR_MSG,only:set_warning_for_getm
  use constants,  ONLY: p_qnUc,MW_C,MW_N,MW_P,MW_Si,MW_H,MW_Si,MW_O

  use SourceFunctions,only: Source_D3_withstate,Source_D3_withgroup
  use mem_Param,ONLY:CalcPhytoPlankton,p_qnUlc,p_qpPhc,p_qR6cQ9x
  use mem_globalfun,   ONLY: MM_power_vector,exp_limit,insw_vector
  use global_interface,ONLY:PhaeocystisCalc
  use mem_Phyto, ONLY:p_qnRc, p_qpRc, p_qsRc, p_lqnlc, p_qplc, p_qslc, &
    p_xqp,p_xqn,p_qlPlc,p_qchlc,p_iRI,p_thdo,p_xsize_m,p_xsize_c,p_pu_ea,p_pu_ra
  use mem_MesoZoo,ONLY:p_qnMec=>p_qnc,p_qpMec=>p_qpc
  use mem_MicroZoo,ONLY:p_qnMic=>p_qnc,p_qpMic=>p_qpc
  use mem_PelBac,ONLY:p_xsizeB_m=>p_xsize_m,p_xsizeB_c=>p_xsize_c, &
                 p_qnBlc=>p_qnlc,p_qpBlc=>p_qplc,p_qnBc=>p_qnc,p_qpBc=>p_qpc
  use mem_Phaeo,only:CALC_SEDIMENTATION,CALC_REL_RADIUS
  use AttachRatetoTEP,only:RelAttachRateToTEP
  use bio_bfm, ONLY: DeriveFromGotm


!  use mem,ONLY: Output2d_1,Output2d_2,Output2d_3, Output2d_4

!
!
! !AUTHORS
!   Piet Ruardij
!
!
! !REVISION_HISTORY
!   Created at Tue Apr 20 09:11:59 AM CEST 2004
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
  integer,intent(IN)              ::mode

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer                         :: i,j,iout,jout,from,ito
  integer,parameter               :: sw_mode=1
  real(RLEN)                      :: rscalar
  real(RLEN),parameter            :: xlow_mass_c=1.0D-5
  real(RLEN),parameter            :: xlow_mass_c_phyto=3.0D+0
  real(RLEN),dimension(NO_BOXES)  :: cx_any,rx_any,px_any,mx_any,ex_any,sx_any
  real(RLEN),dimension(NO_BOXES)  :: r,h,dummy,rec_qR2P1
  real(RLEN),dimension(NO_BOXES)  :: iNN
  real(RLEN),dimension(NO_BOXES)  :: pxMAStickControl !MacroAggregaat Parameters
  real(RLEN),dimension(NO_BOXES)  :: sxMACatchControl
  real(RLEN),dimension(NO_BOXES)  :: E,buoyancy
  real(RLEN),dimension(NO_BOXES)  :: x_macro_c
  real(RLEN),dimension(NO_BOXES)  :: eo,meet
  real(RLEN),dimension(NO_BOXES)  :: ref_meet
  real(RLEN),dimension(NO_BOXES)  :: sedi_macragg
  real(RLEN),dimension(NO_BOXES)  :: nc,nn,np,ns
!                                    number pf macro.agg,size_macroagg
  real(RLEN),dimension(NO_BOXES)  :: x_macro_n,x_cell_n
  real(RLEN),dimension(NO_BOXES,iiPhytoPlankton)  :: PIx_macr_c,R2x_macr_c

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Vectors used  of group vectors
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN), dimension(:), pointer  ::lcl_Ph

  select case (mode)
  case(1)
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Reset var in  betnhic vars which are used in the pelagic
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    jnetY3c(:) =ZERO
    jnetYy3c(:)=ZERO
    jrrPTc(:)=ZERO
    jrrMic(:)=ZERO
    jrrMec(:)=ZERO
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Reset var in which silica and TEP fluxes are collected:
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    flR3R2c(:)= ZERO
    flPTN6r(:)= ZERO

    fl_xgrazing_PIc=ZERO
    jnetPIc(:,:)=ZERO
    jnetB1c(:)=ZERO
    jnetMeZc(:)=ZERO
    jnetMiZc(:)=ZERO
    jZIR6n=ZERO
    jZIR6p=ZERO
    jZIDIn=ZERO
    jZIDIp=ZERO
    jPTZTn=ZERO
    jPTZTp=ZERO
    jPelFishInput(:)=ZERO
    rnetPTc=ZERO
    flPIR6n(:,:)=ZERO
    flPIR1n(:,:)=ZERO
    flPIR6p(:,:)=ZERO
    flPIR1p(:,:)=ZERO
    flBaZTc(:)=ZERO
    flN3N4n(:)=ZERO
    rml(:)=ZERO
    jPLO3c(:)=ZERO

    flR1O3c=ZERO ! due to urea breakdown
    flR1N4n=ZERO ! due to urea breakdown
    flB1N4n=ZERO ! due to urea breakdown
    flPIN4n=ZERO
    flN3N4n=ZERO !reverse nitrification
    POC=ZERO
    PON=ZERO
    POP=ZERO

    !Benthic flux initialized in pelagic,because mortality of larvae  of
    !bnenthic filterfeeders are added to this rate.
    jmY3c=ZERO

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Comput/Estimate Urea and 'rich' labile organic P
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    rscalar=0.01*sum(Depth*R1c)/sum(Depth)
    nc=max(ZERO,R1c-rscalar)
    rscalar=0.01*sum(Depth*R1n)/sum(Depth)
    nn=max(ZERO,R1n-rscalar)
    r = max(ZERO,min(p_qnUlc*nc,nn(:)))
    Nun(:)=max(ZERO, min(nn(:),nc(:)*p_qnUc)-r)
    r = max(ZERO,min(p_qpPhc*R1c(:),R1p(:)))
    Nup(:)=max(ZERO, min(R1p(:),R1c(:)*p_qpPhc)-r)

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Compute nutrient quota in pelagic detritus
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    qpR6c(:)= R6p(:)/(NZERO+ R6c(:))
    qnR6c(:)= R6n(:)/(NZERO+ R6c(:))
    qsR6c(:)= R6s(:)/(NZERO+ R6c(:))
    POC=POC+R6c
    PON=PON+R6n
    POP=POP+R6p

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Compute nutrient quota in microzooplankton and HNAN
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    do i = 1 , (iiMicroZooPlankton)
      cx_any=MicroZooPlankton(i,iiC)
!write(LOGUNIT,*)'pelglobal microzoo'
      if (ppMicroZooPlankton(i,iiN) > 0 )  then
        qn_mz(:,i)= MicroZooPlankton(i,iiN)/(NZERO+cx_any )
      else
!write(LOGUNIT,*)'pelglobal microzoo 2'
        qn_mz(:,i)= p_qnMic(i)
!write(LOGUNIT,*)'pelglobal microzoo 3'
      endif
      if (ppMicroZooPlankton(i,iiP) > 0 ) then
        qp_mz(:,i)= MicroZooPlankton(i,iiP)/(NZERO+cx_any )
      else
        qp_mz(:,i)= p_qpMic(i)
      endif
      POC=POC+cx_any
      PON=PON+cx_any*qn_mz(:,i)
      POP=POP+cx_any*qp_mz(:,i)
    end do
!write(LOGUNIT,*)'pelglobal mesozoo'
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Compute nutrient quota in omnivorous and herbivorous mesozooplankton
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    do i = 1 , (iiMesoZooPlankton)
      cx_any=MesoZooPlankton(i,iiC)
      if (ppMesoZooPlankton(i,iiN) > 0 ) then
        qnZc(:,i)= MesoZooPlankton(i,iiN)/(NZERO+ cx_any)
      else
        qnZc(:,i)= p_qnMec(i)
      endif
      if (ppMesoZooPlankton(i,iiP) > 0 ) then
        qpZc(:,i)= MesoZooPlankton(i,iiP)/(NZERO+ cx_any)
      else
        qpZc(:,i)= p_qpMec(i)
      endif
      POC=POC+cx_any
      PON=PON+cx_any*qnZc(:,i)
      POP=POP+cx_any*qpZc(:,i)
    end do

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Compute nutrient quota in phytoplankton
    ! Compute light prop.or chl. quota in phytoplankton (dep. on ChlLightFlag)
    ! Max quotum is limited to p_xqn*p_nrc because Phaeocystis may contain more
    ! than the maximum. This N is found in the colony but ouside the cells.
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!write(LOGUNIT,*)'pelglobal phyto'

    do i = 1 , iiPhytoPlankton
      cx_any= PhytoPlankton(i,iiC)
      where (cx_any.gt.NZERO)
        qnPc(:,i)= PhytoPlankton(i,iiN)/cx_any
        qpPc(:,i)= PhytoPlankton(i,iiP)/cx_any
        POC=POC+cx_any
        PON=PON+cx_any*qnPc(:,i)
        POP=POP+cx_any*qpPc(:,i)
      elsewhere
        qnPc(:,i)=p_xqn(i)*p_qnRc(i)
        qpPc(:,i)=p_xqp(i)*p_qpRc(i)
      endwhere
      PI_dw(:,i)=p_xsize_c(i) &
          *(DONE+(MW_O+2.0*MW_H+MW_N*qnPc(:,i)+(MW_P+4.0*MW_O)*qpPc(:,i))/MW_C)
      if (i== iiP6 ) then
        !Limit in an artifical way the quotum.
        !If Phaeocystis contain more N an P than maximum
        !the  N and P is found in the colony but outside the cells.
        qnPc(:,i)=min(qnPc(:,i),p_xqn(i)*p_qnRc(i))
        qpPc(:,i)=min(qpPc(:,i),p_xqp(i)*p_qpRc(i))
      endif
      qlPc(:,i)=min(DONE*p_qchlc(i), max(0.1* p_qlPlc(i), &
                  PhytoPlankton(i,iiL)/(NZERO+ PhytoPlankton(i,iiC))))
      qsPc(:,i)=ZERO
      j=ppPhytoPlankton(i,iiS)
      if (j>0) then
         qsPc(:,i) = PhytoPlankton(i,iiS)/(NZERO+ PhytoPlankton(i,iiC))
         PI_dw(:,i)=PI_dw(:,i)+p_xsize_c(i)*((MW_Si+3.0*MW_O)*qsPc(:,i))/MW_C
      endif
      flPIR6s(:,i) = ZERO
    end do

!write(LOGUNIT,*)'pelglobal bac'

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Compute nutrient quota in Pelagic Bacteria
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    qpB1c(:) = B1p(:)/(NZERO+ B1c(:))
    qnB1c(:) = B1n(:)/(NZERO+ B1c(:))
    call findlarge(qpB1c,NO_BOXES,2.0*p_qpBc,iout)
    call findlarge(qnB1c,NO_BOXES,2.0*p_qnBc,jout)

    if (iout.gt.0.or.jout.gt.0) then
       iout=max(jout,iout)
       if (iout.gt.0)write(LOGUNIT,*) 'qpB1c>2*p_qpc:',qpB1c(iout)
       if (jout.gt.0)write(LOGUNIT,*) 'qnB1c>2*p_qnc:',qnB1c(iout)
    endif
!   qpB1c=max(p_qpBlc,min(p_qpBc,qpB1c))
!   qnB1c=max(p_qnBlc,min(p_qnBc,qnB1c))
    POC=POC+B1c
    PON=PON+B1n
    POP=POP+B1p

    ! Compute sedimentation velocities
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    sediPI= ZERO
    sediR2= ZERO
    sediR6= NZERO
    sediB1= NZERO
    sedi_macragg=ZERO

    !limitation of sedimentation rate at very low concentration (1.0E-6)
!write(LOGUNIT,*)'pelglobal macroagg'

  ! concentration TEP determines sinking rate!
    if (p_raPIm> ZERO)  then

      x_macro_c=ZERO
      !STEP1 : determine the size of the macroAggregates.
      !        macroaggrages are formed in presence of R2 and algae which
      !        produce Carbonhydrates =(diatoms)
      do i=1,iiPhytoPlankton
        !Select all diatoms  types and exclude resuspended benthic diatoms
        if (p_iRI(i).eq.iiR2.and.CoupledtoBDc(i).eq.0.and. &
                                              CalcPhytoPlankton(i)) then
          lcl_Ph => PhytoPlankton(i,iiC)
          !first determine nutrient status.range:0-2  (depleted->rich)
          iNN=min((qpPc(:,i)- p_qplc(i))/(p_qpRc(i)- p_qplc(i)), &
              (qnPc(:,i)- p_lqnlc(i))/(p_qnRc(i)- p_lqnlc(i)))
          if (p_qslc(i).gt.ZERO) &
             iNN=min(iNN,(qsPc(:,i)- p_qslc(i))/(p_qsRc(i)- p_qslc(i)))
          iNN=max(ZERO,iNN)
          ! Check if diatoms are really nutrient limited
          !Calculate fraction of algal type which will be part of MacroAggratate
          !nutrient limitation*(quotient potential macroaggregate/phytoplankton)
          ! (unit: -) value between 0-1 valid for the diatoms

!         qhP1R2=(DONE-p_pu_ea(i)-p_pu_ra(i))/p_pu_ea(i)
          qlR2P1=p_pu_ea(i)/(DONE-p_pu_ea(i)-p_pu_ra(i))
!         meet=max(NZERO,R2c-qhP1R2*lcl_Ph/(NZERO+R2c))
          meet=max(NZERO,R2c-lcl_Ph/qlR2P1/(NZERO+R2c))
          !high stickness when nutrients are depleted.
          pxMAStickControl=p_thdo(i)/(p_thdo(i)+iNN)
          ! potential stickyness:
!         qP1R2=lcl_Ph/(NZERO+R2c)
          qR2P1=R2c/(NZERO+lcl_Ph)
          ! set maximum size to macroaggregates by limiting min proportion 
          ! P1 to R2
!         qP1R2=max(qP1R2,p_qlP1R2)
          qR2P1=min(qR2P1,p_qhR2P1)
          ! atucal stickyness: at low nutrient depletion when pxMAStickControl 
          ! is small qP1P2 is decreased
!         qP1R2=qP1R2/pxMAStickControl
          qR2P1=qR2P1*pxMAStickControl
!         qP1R2=min(qhP1R2,qP1R2)
          qR2P1=max(qlR2P1,qR2P1)
!         rec_qP1R2=(DONE/(NZERO+qP1R2))
          rec_qR2P1=(DONE/(NZERO+qR2P1))
          ! Determine the fraction of R2 which is sticked to the MacroAggrate
          ! amount min(free R2 ,(limitation for diatom)
          !  * (amount R2 based on phytoplankton and fixed relation)
!         R2x_macr_c(:,i)= min(R2c,rec_qP1R2*lcl_Ph)
          R2x_macr_c(:,i)= min(R2c,qR2P1*lcl_Ph)
!         PIx_macr_c(:,i)= min(lcl_Ph,R2c*qP1R2)
          PIx_macr_c(:,i)= min(lcl_Ph,R2c*rec_qR2P1)
          !Im sediPI the value of pxMAStickControl is tempory stored
        endif
      enddo
      !STEP2+3 :check only if sums of R2 is not larger then total of R2
      ! if this is the case correct r
      ! calculate maximal sedimnetations rates in case
      ! if MacroAggregates are on its largest.
      do i=1,iiPhytoPlankton
        !Select all diatoms  types and exlcude resuspended benthic datoms
        if (p_iRI(i).eq.iiR2.and.CoupledtoBDc(i).eq.0.and. &
                              CalcPhytoPlankton(i)) then
          lcl_Ph => PhytoPlankton(i,iiC)
          !Im sediPI the value of pxMAStickControl is tempory stored
           ! MacroAggregaat = phytoplankton+ R2c
!          x_macro_c=(R2x_macr_c(:,i) +PIx_macr_c(:,i) )  &
!            *min(DONE,(NZERO+lcl_Ph)/ &
!             max(NZERO,R2c-lcl_Ph/p_qlP1R2/pxMAStickControl) )
           x_macro_c=(R2x_macr_c(:,i) +PIx_macr_c(:,i) )  &
             *min(DONE,(NZERO+lcl_Ph)/ &
              max(NZERO,R2c-lcl_Ph*p_qhR2P1*pxMAStickControl) )
           ! size of macroaggregaat is determined by a fraction of conc. which 
           ! is a macro-aggrgaat divided by the optimum conentration TEP/
!          xSizeMA_m=p_xsize_m(i)*max(DONE,R2x_macr_c(:,i)/p_chR2R2c)
           ! recalculation in dry weight for use in Stoke equation. Weight of 
           ! one cells is multiplied with tota mass of macroaggregaat in dry 
           ! weight and weight of diatom fraction in macroaggregaat.(Implicitly
           ! assuming that macroaggrgaat consists of 1 diatom cell with
           ! around TEP.
           cx_any=PI_dw(:,i)*PIx_macr_c(:,i)/p_xsize_c(i) &
                           +R2x_macr_c(:,i)*(DONE+(2.0*MW_H+MW_O)/MW_C)
           xSizeMA_d=p_xsize_c(i)*max(DONE,cx_any/(NZERO &
             +max(p_xsize_c(i),PIx_macr_c(:,i)) ) )
           xSizeMA_m=p_xsize_m(i)*max(DONE,xSizeMA_d/ &
                                         PI_dw(:,i))**(0.333333)
           buoyancy=pxMAStickControl
           call CalcSinking(NO_BOXES,xSizeMA_m,xSizeMA_d,ETW,ESW,sedi_macragg)
           sedi_macragg=sedi_macragg*buoyancy*p_fMAsink
           sediPI(:,i)=sedi_macragg* PIx_macr_c(:,i)/(NZERO+lcl_Ph)+p_rrPIm(i)
           sediR2=sediR2 + R2x_macr_c(:,i) *sedi_macragg
         elseif (CalcPhytoPlankton(i)) then
           sediPI(:,i)=p_rrPIm(i)
         endif
         call findnan(sediR2,NO_BOXES,iout)
         if (iout>0) then
            write(logunit,*) 'PelGLobal I in sediR2c layer',iout, &
              sedi_macragg(iout),R2c(iout), R2x_macr_c(iout,i), &
              qR2P1(iout),rec_qR2P1(iout),lcl_Ph(iout),R2x_macr_c(iout,i)
        endif
       enddo
       !STEP5
       if (sum(sedi_macragg).gt.NZERO) then
         call DeriveFromGotm(1,NO_BOXES,E)
         sx_any=sedi_macragg/Depth
         x_macro_n=x_macro_c/xSizeMA_d
         do i=1,iiPhytoPlankton
           if (p_iRI(i).ne.iiR2.and.CoupledtoBDc(i).eq.0 &
                                      .and.CalcPhytoPlankton(i)) then
!write(LOGUNIT,*)'pelglobal pheocystis step5, i',i
             call PhaeocystisCalc(CALC_REL_RADIUS,iiP6,px_any,dummy,p_xsize_m(i))
!write(LOGUNIT,*)'pelglobal pheocystis step5 after call'
             mx_any=px_any*p_xsize_m(i)
             sxMACatchControl=RelAttachRateToTEP(1,NO_BOXES,mx_any,E, &
                radius_macro=xSizeMA_m,number_macros=x_macro_n)
             px_any=meet*pxMAStickControl*sxMACatchControl/(NZERO+sx_any)
             sediPI(:,i)=sediPI(:,i) &
                    +min(DONE,px_any)*sedi_macragg
           endif
         enddo
         mx_any=p_xsizeB_m
         sxMACatchControl=RelAttachRateToTEP(1,NO_BOXES,mx_any,E, &
          radius_macro=xSizeMA_m,number_macros=x_macro_n)
         px_any=meet*pxMAStickControl*sxMACatchControl/(NZERO+sx_any)
         sediB1(:)=sediB1(:) +min(DONE,px_any)*sedi_macragg
         call findnan(sediB1,NO_BOXES,iout)
         if (iout>0) then
          write(logunit,*) 'PelGLobal II in sediB1c layer',iout, &
          sedi_macragg(iout),R2c(iout), R2x_macr_c(iout,i), &
          qR2P1(iout),lcl_Ph(iout),R2x_macr_c(iout,i),&
           px_any(iout),sxMACatchControl(iout),meet(iout),pxMAStickCOntrol(iout)
         endif
       endif
     endif
!    sediPI(iiP5,:)=35.0*P5c**(0.9647)
!    sediPI(iiP5,:)=2.0*sediR9
     sediPI(:,iiP5)=86400.0*0.0022
    ! Calculate sinking rate (dummy is here a dummy variable) on basis of size and wight of colonies.
    ! Compare sinking with siinkin due to sticking to macroaggregates produced by diatoms.
!write(LOGUNIT,*)'pelglobal pheocystis calc_sedimentation',rx_any
!stop
!JM return leads to segmentation fault. Skip for now, solve later
   call PhaeocystisCalc(CALC_SEDIMENTATION,iiP6,rx_any,dummy,p_xsize_m(iiP6))
!write(LOGUNIT,*)'pelglobal sediPI'
!stop
    sediPI(:,iiP6)=max(sediPI(:,iiP6),rx_any)
!write(LOGUNIT,*)'pelglobal sediPI 2'
!stop

     do i = 1 , iiPhytoPlankton
      if (CalcPhytoPlankton(i)) &
        sediPI(2:NO_BOXES,i)= &
                       (sediPI(2:NO_BOXES,i)+sediPI(1:NO_BOXES-1,i))*0.5
     end do

     p_xfree_R2=DONE
     where (R2c(:) > 1.0D-10 )
      p_xfree_R2=max(ZERO,R2c-R2x_macr_c(:,iiP1))/(NZERO+R2c)
      sediR2(:)=(sediR2(:)+max(ZERO,R2c-R2x_macr_c(:,iiP1))*p_raR2m)/(NZERO+R2c)
      sediR2(2:NO_BOXES)=(sediR2(2:NO_BOXES)+sediR2(1:NO_BOXES-1))*0.5
    endwhere
!write(LOGUNIT,*)'pelglobal findnan'
!stop
    call findnan(sediR2,NO_BOXES,iout)
    if (iout>0) then
        write(logunit,*) 'PelGLobal II in sediR2c layer',iout
    endif

    do i = 1, iiMesoZooPlankton
      if (p_rMem(i) < ZERO ) then
        eo= MM_power_vector(max(NZERO,O2o(:)),  p_clMeO2o(i),3)
        sediMeZ(:,i)=p_rMem(i)* (DONE-eo)
        sediMeZ(NO_BOXES_Z,i)=ZERO
      else
        sediMeZ(:,i)=p_rMem(i)
      endif
    end do

    do i = 1, iiMicroZooPlankton
      eo= MM_power_vector(max(NZERO,O2o(:)),  p_clMiO2o(i),3)
      sediMiZ(:,i)=p_rMim(i)* (DONE-eo)
      sediMiZ(NO_BOXES_Z,i)=ZERO
    end do

!write(LOGUNIT,*)'pelglobal sediR6'
!stop
    !Sedmentation of R6 "uncoupled" to silt.
    !we assume that filterfeeders only take up the uncoupled detritus
    sediRZ(:)=p_raRZm*OCDepth(1)/(5.0+OCDepth(1))
    sediR6(:)=p_rrR6m

    !During period of stratification there is above the startified layers no
    ! coupled transport of detritus.
    !Sedimentation of detritus depend on the origin of the  detritus recently
!write(LOGUNIT,*)'pelglobal sourced3 0'
!stop
    r=Source_D3_withstate(ppR6c,ppB1c,iiProduction) &
     +Source_D3_withgroup(ppR6c,ppMicroZooPlankton,iiMicroZooPlankton, &
              iiC,iiProduction)
!write(LOGUNIT,*)'pelglobal sourced3 1'
!stop
    h=r+Source_D3_withgroup(ppR6c,ppPhytoPlankton,iiPhytoPlankton, &
              iiC,iiProduction)
!write(LOGUNIT,*)'pelglobal sourced3 2'
!stop
    r=r +Source_D3_withstate(ppR6c,ppP3c,iiProduction) &
        +Source_D3_withstate(ppR6c,ppZ6c,iiProduction)
!write(LOGUNIT,*)'pelglobal sourced3'
!stop

    do j=1,NO_BOXES_XY
      from=PelBoxAbove(j)
      ito= from-1+NO_BOXES_Z
      !Be sure that sedimentation rate in the 2 lowest layers is not changed
      !upper side of startification is limited to for NO_BOXES_Z-2 layers
      rscalar=min(OCDepth(2),-BMLd(j))
      where (OCDepth(from:ito).lt.rscalar) sediR6(from:ito)= &
         max(ZERO,(NZERO+h(from:ito)-r(from:ito))/(NZERO+h(from:ito)))*p_rrR6m
    enddo

    sediR6(:)= (sediR6(:)*max(ZERO,R6c(:)-RZc(:)) &
             +p_raRZm*(NZERO+RZc(:)))/ (NZERO+ max(RZc(:),R6c(:)))
    sediR6s(:)=(sediR6s(:)*max(ZERO,R6s(:)-RZc(:)*qsR6c) &
             +p_raRZm*(NZERO+RZc(:)*qsR6c))/ (NZERO+ max(RZc(:)*qsR6c,R6s(:)))

    !limit sedimentation at very small values of P2 (avoid problems in GOTM)
    do i = 1 , iiPhytoPlankton
     if (CalcPhytoPlankton(i) ) then
       cx_any=PhytoPlankton(i,iiC)
       rscalar=max(ZERO,sum(cx_any*Depth))/OCDepth(1)
!JM       sediPI(:,i)= sediPI(:,i) *max(ZERO, &
!JM         (rscalar/(rscalar+ xlow_mass_c))- (xlow_mass_c/(rscalar+ xlow_mass_c)))
       sediPI(:,i)= sediPI(:,i) *max(ZERO, &
         (rscalar/(rscalar+ xlow_mass_c_phyto))- (xlow_mass_c_phyto/(rscalar+ xlow_mass_c_phyto)))
     else
       sediPI(:,i)=ZERO
     endif
    enddo
!write(LOGUNIT,*)'pelglobal limit'
!stop

    !limit sedimentation at very small values of R2 (avoid problems in GOTM)
    rscalar=max(ZERO,sum(R2c*Depth))
    sediR2(:)= sediR2(:)*max(ZERO, (rscalar/(rscalar+ xlow_mass_c))- &
                                          (xlow_mass_c/(rscalar+ xlow_mass_c)))
    !limit sedimentation at very small values of R6 (avoid problems in GOTM)
    rscalar=max(ZERO,sum(R6c*Depth))
    sediR6(:)= sediR6(:)*max(ZERO, (rscalar/(rscalar+ xlow_mass_c))- &
                                       (xlow_mass_c/(rscalar+ xlow_mass_c)))
    rscalar=max(ZERO,sum(R6s*Depth))
    sediR6s(:)= sediR6s(:)*max(ZERO, (rscalar/(rscalar+ xlow_mass_c))- &
                                         (xlow_mass_c/(rscalar+ xlow_mass_c)))
  case(2)
    ! Here the rate only for Settling is calulated for which is
    ! assumed that a part of the R6 is fixed to Silt.
    px_any=max(ZERO,min(p_qR6cQ9x*R9x,R6c-RZc)/(NZERO+R6c-RZc))
    sediR6(:)=sediR9*px_any+sediR6*(DONE-px_any)
    sediR6s=sediR6

    rscalar=max(ZERO,sum(R6c*Depth))
    sediR6(:)= sediR6(:)*max(ZERO, (rscalar/(rscalar+ xlow_mass_c))- &
                                       (xlow_mass_c/(rscalar+ xlow_mass_c)))
    rscalar=max(ZERO,sum(R6s*Depth))
    sediR6s(:)= sediR6s(:)*max(ZERO, (rscalar/(rscalar+ xlow_mass_c))- &
                                         (xlow_mass_c/(rscalar+ xlow_mass_c)))
    !----------------------------------------------------
    ! Calculation of nutrient limitation index
    ! 1=n-limitation, 5=p-limitation,3=si-limitations,
    ! 2= si-limitation for diatoms,n-limitation for other
    ! 4= si-limitation for diatoms,p-limitation for other
    !----------------------------------------------------

    nn=ZERO;np=ZERO;ns=ZERO;cx_any=ZERO;h=ZERO
    do i=1,iiPhytoPlankton
        nn=nn+ PhytoPlankton(i,iiN)/p_qnRc(i)
        np=np+ PhytoPlankton(i,iiP)/p_qpRc(i)
        if (ppPhytoPlankton(i,iiS) >0) then
          ns=ns+ PhytoPlankton(i,iiS)/p_qsRc(i)
          h=h+ PhytoPlankton(i,iiC)
        endif
        cx_any=cx_any * PhytoPlankton(i,iiC)
    enddo


    nn=nn/(NZERO+cx_any)
    np=np/(NZERO+cx_any)
    ns=ns/(NZERO+h)
    limnuti=ZERO
    h=ZERO
    where (ns < DONE .and. ns<nn .and. ns<np  ) h=3.0
    where (nn < DONE .and.nn < np )
      limnuti=DONE
    elsewhere (np < DONE .and. np < nn )
      limnuti=5.0
    endwhere
    where (h>ZERO .and. limnuti>ZERO)
      limnuti=(limnuti + h)* 0.5
    elsewhere(h > ZERO)
      limnuti=h
    endwhere

!  kwdratic sums
    O2o_vr=O2o*O2o
    N1p_vr=N1p*N1p
    N3n_vr=N3n*N3n
    Chla_vr=Chla*Chla
    R6c_vr=R6c*R6c
    !----------------------------------------------------
    ! Calculation of Phyto functional group index
    ! nr of functional group with the hihgest biomass which is larger tan 5.0
    !----------------------------------------------------
    cx_any=ZERO; h=ZERO ;PTi=ZERO
    do i=1,iiPhytoPlankton
      h= PhytoPlankton(i,iiC)
      where (h> 5.0 .and. h>cx_any)
         cx_any=h; PTi=real(i)
      endwhere
    enddo
  end select
  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
