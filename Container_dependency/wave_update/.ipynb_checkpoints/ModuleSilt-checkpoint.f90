!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: Silt
!
! DESCRIPTION

! !INTERFACE
  module mem_Silt
!
! !USES:

  use global_mem
  use mem,only: iiPELSINKREF,ppR9x,NO_BOXES

!  
!
! !AUTHORS
!   Original version by P. Ruardij and M. Vichi
!
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
  ! Default all is public
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  public
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Silt PARAMETERS (read from nml)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! method: 1:  Piet's old method; removed
  !         2:  limit sedimentation for shallow water
  !         3:  limit sedimentation for shallow water (updated method as used in JSR paper)
  !         4:  much simpler and more universal limitation of sedimentation rates (preferred)
  ! p_SampleDepth: obsolete
  ! stick_fact:    factor to limit resuspension assuming a proportion of (1-stick_fact)
  !                remains tied in with the coarse fraction. Default value is for 
  !                3D runs. For 1D runs use 0.25
  integer     :: method
  real(RLEN)  :: p_SampleDepth=0.01, stick_fact=0.15
  real(RLEN),dimension(:),allocatable    :: start_R9x

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Wave formulation used by the Silt resuspension (see Silt/wave.F90)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! wave_method: 1:  jonswap() as originally compiled (default)
  !              2:  jonswap() with the depth-duration law and breaker ratio
  !                  from this namelist (wave_p, wave_Tdmax, wave_gamma)
  !              3:  effective-fetch scheme: Breugem & Holthuijsen local sea on
  !                  the directional effective fetch, plus propagated sea on the
  !                  exposure fetch; SPM peak period on the period fetch
  !              4:  jonswap() height with an effective-fetch cap; period as 1
  !  Methods 3 and 4 need the fetch fields: wave_fetch_file in getm_bio.inp
  !  (&getm_bfm_wave_nml), created by DWS_CLI/tools/wave_fetch_grid.py.
  ! wave_p, wave_Tdmax: depth exponent and maximum duration (s) of
  !                  Tdm = (depth/10)^p * Tdmax (methods 2 and 4)
  ! wave_gamma:      breaker ratio, Hs <= wave_gamma*depth; <0: 0.4 for
  !                  methods 2 and 4, 0.25 for method 3
  ! wave_alpha_exp:  weight of the propagated sea (method 3)
  ! wave_Tdmax_exp:  duration limit of the propagated sea (s, method 3)
  ! wave_Tdmax_tp:   duration limit of the peak period (s, method 3)
  ! wave_Tz_min:     lower limit of the returned period (s, methods 2-4)
  ! Calibrated values (DWS 500 m, 2015): method 3 with the defaults below;
  ! method 2 with wave_p=1.7, wave_Tdmax=3600.
  integer, parameter :: NBEAR_FETCH=36
  integer     :: wave_method=1
  real(RLEN)  :: wave_p=2.0D0, wave_Tdmax=10800.0D0, wave_gamma=-1.0D0
  real(RLEN)  :: wave_alpha_exp=0.8D0, wave_Tdmax_exp=18000.0D0, wave_Tdmax_tp=18000.0D0
  real(RLEN)  :: wave_Tz_min=1.0D0

  ! Column values for the fetch-based methods, set by the 3D driver for every
  ! water column through set_wave_column() before the biology is called.
  !   wave_convc:      grid convergence (deg), compass bearing of grid north
  !   wave_fetch_dir:  directional effective fetch (m) for bearings the wind
  !                    blows FROM, 0,10,...,350 deg geographic
  !   wave_fetch_cap:  as wave_fetch_dir with the method-4 blocking depth
  !   wave_fetch_exp:  direction-mean exposure fetch (m)
  !   wave_fetch_tp:   direction-mean period fetch (m)
  real(RLEN)  :: wave_convc=0.0D0
  real(RLEN)  :: wave_fetch_dir(NBEAR_FETCH)=-1.0D0, wave_fetch_cap(NBEAR_FETCH)=-1.0D0
  real(RLEN)  :: wave_fetch_exp=-1.0D0, wave_fetch_tp=-1.0D0
  logical     :: wave_fetch_set=.false., wave_cap_set=.false.


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! SHARED PUBLIC FUNCTIONS (must be explicited below "contains")

  public InitSilt, set_wave_column
  contains

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  subroutine InitSilt()

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer              :: status
  namelist /Silt_parameters/ method,p_SampleDepth,stick_fact, &
                             wave_method,wave_p,wave_Tdmax,wave_gamma, &
                             wave_alpha_exp,wave_Tdmax_exp,wave_Tdmax_tp, &
                             wave_Tz_min
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !BEGIN compute
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !  Open the namelist file(s)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

   allocate(start_R9x(1:NO_BOXES),stat=status)
   if (status /= 0) call error_msg_prn(ALLOC,"Silt","Start_R9x")
   start_R9x=ZERO
   iiPELSINKREF(ppR9x)=-ppR9x  
   write(LOGUNIT,*) "#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-"
   write(LOGUNIT,*) "#  Reading Silt parameters.."
   open(NMLUNIT,file='Silt.nml',status='old',action='read',err=100)
    read(NMLUNIT,nml=Silt_parameters,err=101)
    close(NMLUNIT)
    write(LOGUNIT,*) "#  Namelist is:"
    write(LOGUNIT,nml=Silt_parameters)
    if (wave_method.lt.1 .or. wave_method.gt.4) then
      write(LOGUNIT,*) "#  Silt: wave_method must be 1, 2, 3 or 4; found ",wave_method
      call error_msg_prn(NML_READ,"InitSilt.f90","Silt_parameters: wave_method")
    endif
    if (wave_gamma.lt.ZERO) then
      if (wave_method.eq.3) then
        wave_gamma=0.25D0
      else
        wave_gamma=0.4D0
      endif
    endif
    write(LOGUNIT,*) "#  Silt: wave_method=",wave_method," wave_gamma=",wave_gamma
!   if (sw_adv==0 ) then
!     iiPELSINKREF(ppR9x)=-1000000 
!     write(LOGUNIT,*) "#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-"
!     write(LOGUNIT,*) "Vertical advective transport is excluded for R9x (Silt)"
!     write(LOGUNIT,*) "#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-"
!   endif
    
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !END compute
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  return
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Error Messages
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
100 call error_msg_prn(NML_OPEN,"InitSilt.f90","Silt.nml")
101 call error_msg_prn(NML_READ,"InitSilt.f90","Silt_parameters")
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  end  subroutine InitSilt

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Store the fetch fields and grid convergence of the current water column.
  ! Called by the 3D driver (GETM getm_bio.F90) before the biology of each
  ! column when a wave fetch file is used.
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  subroutine set_wave_column(convc,fetch_dir,fetch_exp,fetch_tp,fetch_cap)
  real(RLEN),intent(in)          :: convc,fetch_exp,fetch_tp
  real(RLEN),intent(in)          :: fetch_dir(NBEAR_FETCH)
  real(RLEN),intent(in),optional :: fetch_cap(NBEAR_FETCH)

  wave_convc=convc
  wave_fetch_dir=fetch_dir
  wave_fetch_exp=fetch_exp
  wave_fetch_tp=fetch_tp
  wave_fetch_set=.true.
  if (present(fetch_cap)) then
    wave_fetch_cap=fetch_cap
    wave_cap_set=.true.
  endif
  end subroutine set_wave_column

  end module mem_Silt
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
