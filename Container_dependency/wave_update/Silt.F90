#include "DEBUG.h"
#include "INCLUDE.h"
#include"cppdefs.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: Silt
!
! DESCRIPTION
!
!

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine SiltDynamics
!
! !USES:


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN,ZERO
#ifdef nopointers
  use mem,  ONLY: D3STATE
#else
  use mem,  ONLY: Q9x,Qp9x,R9x, ppR9x,ERHO,Wind,ETW,ETAUB,EUWIND,EVWIND, &
                  EUCURR_LEVEL,EVCURR_LEVEL,R2c,sediR9,NO_BOXES_XY, &
                  Hs_out,Tz_out,u_orb_out,TauW,TauC,TauBed,ws_out,eta_out
  use mem,  ONLY: iiPelSinkREF
#endif
  use mem, ONLY: Depth,ETW, NO_BOXES, iiBen, iiPel, flux_vector,InitializeModel
  use mem, ONLY: jrESS,dry_z,Depth

  use mem_Param,  ONLY: p_small, p_poro,p_dry_ben
  use global_interface,   ONLY: eTq
  use mem_Silt
  use LimitRates, ONLY:LimitChange

  use wave, ONLY: do_wave
  use spm_util, ONLY: wave_parameters, reference_concentration
  use bio_var, ONLY: nuh_l,ws
  use LimitRates,ONLY:LimitChange
  use constants,only:ANY,SEC_PER_DAY,NEGATIVE


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following vector functions are used:MM_vector, eTq_vector, insw_vector
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

!  
!
! !AUTHORS
!   Original version by P. Ruardij and  J. van der Molen
!
!
!
! !REVISION_HISTORY
!   ! 2007: JM introduced simple wave-driven resuspension model
!   ! 2009: JM introduced concentration-dependent settling rate
!   ! April 2010: JM: changed to work without using rate calculation
!                     changed to produce linear depth dependence
!   ! Nov 2013: JM: cleaned-up, new version calculating bottom concentration
!   !               and relying on getm/gotm for horizontal and vertical transport
!
! COPYING
!   
!   Copyright (C) 2006 P. Ruardij, the mfstep group, the ERSEM team 
!   (rua@nioz.nl, vichi@bo.ingv.it)
!
!   This program is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation;
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



  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
     REAL(RLEN),dimension(NO_BOXES_XY) :: psilt !fraction silt in sediment
     REAL(RLEN)                        :: tdepth,delta
     REAL(RLEN)                        :: U, Hs, Tz, tau, uw, R9x_new, R9x_old
     REAL(RLEN)                        :: tau_w, tau_c, tau_bed, eta
     REAL(RLEN)                        :: R9x_surf, R9x_bot, R9x_grad
     REAL(RLEN), PARAMETER             :: g=9.81, fw=0.1
     REAL(RLEN), PARAMETER             :: grainsize=250*1E-6, ws_max=0.002, rho_s=2650, nu_kin_ref=1.0E-6, t_ref=20.0
     REAL(RLEN), PARAMETER             :: d_eromin=0.01, dlim1=10.0D0, dlim2=5.0D0, dlimi=dlim2+1.0D+00
     REAL(RLEN)                        :: silt_dry,nu_kin
     REAL(RLEN), PARAMETER             :: kg_to_mg=1E6    !unit conversion factor kg to mg
     REAL(RLEN), PARAMETER             :: a_settle=0.15, b_settle=1.1, kappa=0.4, k_timelag=1000000.0 !a_settle=40.0E-5
! now read from Silt.nml     REAL(RLEN), PARAMETER             :: stick_fact=0.25 !0.15 !0.25! assume that part of the silt sticks to sand
                                                        ! and is not available for resuspension
     INTEGER                           :: n
!    INTEGER, PARAMETER                :: suspmode=3 ! 1: waves, 2: currents, 3: both
     INTEGER, PARAMETER                :: refconc_mode=1
     INTEGER, PARAMETER                :: tau_mode=2
     INTEGER, PARAMETER                :: tauc_mode=1
     INTEGER, PARAMETER                :: tauw_mode=1
     INTEGER, PARAMETER                :: taucrit_mode=0
     REAL(RLEN)                        :: phi_c,phi_w,phi_wc
     REAL(RLEN)                        :: b,z,cref,z_a,u_orb,u_star,cbottmax
     REAL(RLEN)                        :: lim_1,mass_t_m2,lim_2,alpha,r
     REAl(RLEN)                        :: lim_resus,av_ws_uncorr
 
     real(RLEN), external  :: GetDelta
     character :: ch

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
     jrESS=ZERO

     call CouplingBioSilt(0,psilt)
     delta=GetDelta()
!    psilt=max(0.0,(p_poro(1) - 0.38662 )/ 0.415)
     tdepth=sum(Depth)
     nu_kin=nu_kin_ref*2.0/(1+ETW(1)/t_ref)   !van Rijn, 1994
     !allow only "dry" tidal flats if p_dry_ben==.true. (p_dry_ben defined in param.nml 
     ! (default=.false.))
     silt_dry=DONE; if ( p_dry_ben.and.method.ge.2) silt_dry=dry_z(1)
!    method =2 : Johan's official implementation
!    method =3 : Modification of method 2 with sedimentation in shallow areas with 
!                a depth smaller than dlim_2
!    method =4 : netto change in Silt is limited by LimitChange function with a condition that
!                the total silt in the water column may change between  -[silt]/delta and +[silt]/delta
     select case (method)
       case (1)          ! Piet's old method; removed
       case (2,3,4)          ! operational method
         start_R9x=R9x
!LEVEL1 'Silt start: R9x(1)',R9x(1)
         ! wave calculation; wave_method from Silt.nml (default 1 = jonswap)
         call do_wave(wave_method,Hs,Tz,phi_w,Wind,EUWIND(1),EVWIND(1),tdepth)
         call wave_parameters(Hs,Tz,EUCURR_LEVEL(1),EVCURR_LEVEL(1), &
                                                   phi_w,phi_wc,u_orb,tdepth)

         !JM EUCURR_LEVEL and EVCURR_LEVEL are actually bottom currents
         ! reference concentration for sand fraction
         call reference_concentration(cref,z_a,refconc_mode,tau_mode, &
   &           tauc_mode,tauw_mode,taucrit_mode,tau_w,tau_c,tau_bed,nu_kin, &
   &           phi_wc,u_orb,u_star,ERHO(1),rho_s,grainsize, &
   &           EUCURR_LEVEL(1),EVCURR_LEVEL(1),Depth(1),Tz,eta)

         ! concentration-dependent settling velocity; assume power law

!JM          ws(ppR9x,1:NO_BOXES)=-a_settle*((R9x(1:NO_BOXES)/kg_to_mg)**b_settle) 
          ws(1:NO_BOXES,ppR9x)=-a_settle*((R9x(1:NO_BOXES)/kg_to_mg)**b_settle) 
!JM          where (-ws(ppR9x,:)>ws_max) ws(ppR9x,:)=-ws_max
          where (-ws(:,ppR9x)>ws_max) ws(:,ppR9x)=-ws_max
         ! limit settling velocity in shallow situations to prevent instability 
!        of the model
         lim_resus=1.0
         if ( method.eq.2) then
           ! calculate  limiting of  sedimentation for depths dlim1->dlim2
           if (tdepth<dlim1 ) then
!JM             av_ws_uncorr=NZERO+sum(ws(ppR9x,1:NO_BOXES)*Depth(:))
             av_ws_uncorr=NZERO+sum(ws(1:NO_BOXES,ppR9x)*Depth(:))
             if (tdepth<dlim1 .and. tdepth>dlim2) &
!JM               where (ws(ppR9x,1:NO_BOXES)<-ws_max*(1.0D0-((dlim1-tdepth)/(dlim1-dlim2))))  &
!JM                      ws(ppR9x,1:NO_BOXES)=-ws_max*(1.0D0-((dlim1-tdepth)/(dlim1-dlim2)))
               where (ws(1:NO_BOXES,ppR9x)<-ws_max*(1.0D0-((dlim1-tdepth)/(dlim1-dlim2))))  &
                      ws(1:NO_BOXES,ppR9x)=-ws_max*(1.0D0-((dlim1-tdepth)/(dlim1-dlim2)))
           !assume that below depth dlim2 sedimentation rate is 0.
!JM             if (tdepth<=dlim2) ws(ppR9x,:)=0.0D0
!JM             lim_resus=sum(ws(ppR9x,1:NO_BOXES)*Depth(:))/av_ws_uncorr
             if (tdepth<=dlim2) ws(:,ppR9x)=0.0D0
             lim_resus=sum(ws(1:NO_BOXES,ppR9x)*Depth(:))/av_ws_uncorr
           endif
         elseif( method.eq.3) then
           ! dlimi is value between dlim1 and dlim2
           ! calculate alpha on basis of the fact that the value for 
           ! dlimi is equal to the one in the original calculation. With 
           ! this alpha a exponential decrease of the limitation
           ! with respect to depth (dlimi--->0 )  is described.
           r= (1.0D0-((dlim1-dlimi)/(dlim1-dlim2)))
           alpha=log(r)/dlimi ! alpha will be a negative value
           ! use orginal calculation for limiting sedimentation for depths 
           ! dlim1->dlimi
           if (tdepth<dlim1 ) then
!JM             av_ws_uncorr=NZERO+sum(ws(ppR9x,1:NO_BOXES)*Depth(:))
             av_ws_uncorr=NZERO+sum(ws(1:NO_BOXES,ppR9x)*Depth(:))
!JM             if (tdepth<dlim1 .and. tdepth>=dlimi) ws(ppR9x,1:NO_BOXES)= &
!JM               max(ws(ppR9x,:),-ws_max*(1.0D0-(dlim1-tdepth)/(dlim1-dlim2))) 
             if (tdepth<dlim1 .and. tdepth>=dlimi) ws(1:NO_BOXES,ppR9x)= &
               max(ws(:,ppR9x),-ws_max*(1.0D0-(dlim1-tdepth)/(dlim1-dlim2))) 
           ! limit sedimentation for depths dlimi->0
!JM             if (tdepth<dlimi)ws(ppR9x,1:NO_BOXES)=  &
!JM               max(ws(ppR9x,1:NO_BOXES),-ws_max*(exp(alpha*(dlim1-tdepth)))) 
!JM             lim_resus=sum(ws(ppR9x,1:NO_BOXES)*Depth(:))/av_ws_uncorr
             if (tdepth<dlimi)ws(1:NO_BOXES,ppR9x)=  &
               max(ws(1:NO_BOXES,ppR9x),-ws_max*(exp(alpha*(dlim1-tdepth)))) 
             lim_resus=sum(ws(1:NO_BOXES,ppR9x)*Depth(:))/av_ws_uncorr
           endif
         endif
!LEVEL1 'Silt: ws(1)',ws(ppR9x,1)
!LEVEL1 'Silt: cref',cref

         ! update bottom concentration; use rouse profile
         z=Depth(1)/2
!JM         b=-ws(ppR9x,1)/(kappa*sqrt(max(tau_bed,1.0D-10)/ERHO(1)))
         b=-ws(1,ppR9x)/(kappa*sqrt(max(tau_bed,1.0D-10)/ERHO(1)))
!LEVEL1 'Silt: b',b
         !limit uptake considering active bottom layer
         cbottmax=lim_resus*psilt(1)*(max(0.5D+00*eta,d_eromin)/tdepth)*rho_s
!LEVEL1 'Silt: cbottmax',cbottmax
!LEVEL1 'Silt: z',z
!LEVEL1 'Silt: z_a',z_a

         !if grid-point is nearly dyr no resuspension fo silt: R9x keep it old value.
         if (silt_dry > 0.05 ) then
         ! Rouse profile to calculate concentration in bottom grid cell
            R9x(1)=kg_to_mg * stick_fact * min(cbottmax,psilt(1)*cref)  &
                               *  ((z/z_a)*(tdepth-z_a)/(tdepth-z))**(-b)
         endif
         if (method.eq.4) then
           !prevent negative values on settling
           R9x(1)=max(R9x(1),0.4*start_R9x(1))  
           !prevent too positivee values on settling
           R9x(1)=R9x(1)*1.0D6/(R9x(1)+1.0D6)
         endif

!JM           ws(ppR9x,0)=ws(ppR9x,1)
           ws(0,ppR9x)=ws(1,ppR9x)
! interpolate ws to (top) cell boundaries
!JM           sediR9(1:NO_BOXES)=-ws(ppR9x,1:NO_BOXES)*SEC_PER_DAY
!JM           ws(ppR9x,1:NO_BOXES-1)=ws(ppR9x,1:NO_BOXES-1) + &
!JM             (ws(ppR9x,2:NO_BOXES)-ws(ppR9x,1:NO_BOXES-1)) &
!JM                *(0.5*Depth(1:NO_BOXES-1)/ &
!JM                (0.5*Depth(1:NO_BOXES-1)+0.5*Depth(2:NO_BOXES))) 
           sediR9(1:NO_BOXES)=-ws(1:NO_BOXES,ppR9x)*SEC_PER_DAY
           ws(1:NO_BOXES-1,ppR9x)=ws(1:NO_BOXES-1,ppR9x) + &
             (ws(2:NO_BOXES,ppR9x)-ws(1:NO_BOXES-1,ppR9x)) &
                *(0.5*Depth(1:NO_BOXES-1)/ &
                (0.5*Depth(1:NO_BOXES-1)+0.5*Depth(2:NO_BOXES))) 

         if (R9x(1)<0.0) then
           write(LOGUNIT,*)'Silt.F90: R9x(1)<0: ',R9x(1)
           stop
         endif

         ! sediment uptake
         if (silt_dry > 0.05 ) then
           !if Qp9x is still zero answer will be zero
           !detritus resuspension run 1 step behind............
           ! jrESS >0 resuspension
           ! jrESS <0 sedimentation
           jrESS=-(Qp9x-Q9x)/delta
!          jrESS=-min(ZERO,Qp9x-Q9x)/delta
         endif

         ! output
        Hs_out=Hs
        Tz_out=Tz
        u_orb_out=u_orb
        TauW=tau_w
        TauC=tau_c
        TauBed=tau_bed
!JM        ws_out=ws(ppR9x,:)
        ws_out=ws(:,ppR9x)
        eta_out=eta
 
        select case (InitializeModel)         
          case (0)
          case (1)
            R9x(:)=1000.0
        end select

     end select

  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
