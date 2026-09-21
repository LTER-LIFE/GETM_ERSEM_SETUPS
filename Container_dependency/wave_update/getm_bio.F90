#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: getm_bio()
!
! !INTERFACE:
   module getm_bio
!
! !DESCRIPTION:
!
! !USES:
   use parameters, only: g,rho_0
   use domain, only: imin,imax,jmin,jmax,kmax
   use domain, only: az
   use time, only: secondsofday
   use variables_2d, only: D
   use variables_3d, only: taub
   use variables_3d, only: uu,vv,ww,hun,hvn,ho,hn
   use variables_3d, only: nuh,T,S,rad,rho,light
#ifdef BFM_GOTM
#if defined(SPHERICAL) || defined(CURVILINEAR)
   use domain, only: dxu,dxv,dyu,dyv,arcd1
#else
   use domain, only: dx,dy,ard1
#endif
   use domain, only: ioff,joff
   use domain, only: az,au,av
   use variables_bio_3d, only: cc3d,ffp,ffb,bio_missing,d3_pelvar_type
   use variables_bio_3d, only: ccb3d,sw_CalcPhyto_2d, &
                           cc3d_out,ccb3d_out, counter_ave,flag_out, &
                           n_cc3d_out,n_ccb3d_out,adv3d_courant,adv3d_number
#ifdef INCLUDE_DIAGNOS_PRF
   use variables_bio_3d, only: ccb3d_prf,n_ccb3d_prf
#endif
   use variables_3d, only: taub,dt
!JM uncomment this for calculation horizontal fluxes, currently only works for tracking
  use variables_bio_3d, only: cut,cvt      !JM added
#else
   use variables_3d, only: cc3d
#endif
   use meteo, only: swr,u10,v10
   use halo_zones, only: update_3d_halo,wait_halo,D_TAG,H_TAG
#ifdef BFM_GOTM
   use bio, only: init_bio, do_bio
   use bio, only: bio_calc
   use bio_var, only: numc
   use bio_bfm,only:set_env_bio_bfm !JM in bfm2016 this is done in set_env_bio_bfm: ,copy_from_gotm_to_bfm
   use time,only:julianday,start
   use bio_var, only: numbc,ccb,c1dimz, &
          numc_diag,diag, numbc_diag,diagb, bio_setup, &
          adv1d_courant,adv1d_number,cc
   use bfm_output, only: var_ave, var_ids,var_names,&
          stPelStateS,stPelDiagS,stPelFluxS,stBenStateS,stBenDiagS,stBenFluxS, &
          stPelStateE,stPelDiagE,stPelFluxE,stBenStateE,stBenDiagE,stBenFluxE
#ifdef INCLUDE_DIAGNOS_PRF
   use bio_var, only: numbc_prf,nprf
   use bfm_output, only: stPRFDiagS,stPRFDiagE
#endif
   use gotm_error_msg, only: set_d3_model_flag_for_gotm,  &
                                     output_gotm_error,get_warning_for_getm
   use exceptions, only:getm_error
   use coupling_getm_bfm,only:fill_diagn_bfm_vars,DIAG_ADD,DIAG_INFO, &
      getm_bfm_bennut_calc_initial,set_2d_grid_parameters,init_2d_grid, &
      make_uv_flux_output,read_poro,diag_end_sections
   use coupling_getm_bfm_rivers,only:make_river_flux_output
   use mem,only:iiPhytoPlankton,sw_CalcPhyto,R6c
#else
   use bio, only: init_bio, init_var_bio, set_env_bio, do_bio
   use bio, only: bio_calc
   use bio_var, only: numc
   use bio_var, only: cc,ws,var_names,var_units,var_long
#endif
   IMPLICIT NONE
!
! !PUBLIC DATA MEMBERS:
   public init_getm_bio, do_getm_bio
   integer, public           :: bio_init_method=0
#ifdef BFM_GOTM
!  only for BFM:
   integer                   :: calc_init_bennut_states=0
   integer                   :: calc_uv_fluxes
   logical, public           :: hotstart_bio=.true.
   logical,public            :: no_hortr_silt=.false.
!  effective-fetch fields for the BFM Silt wave scheme (Silt.nml wave_method
!  3 and 4), read from wave_fetch_file (getm_bio.inp, &getm_bfm_wave_nml)
   character(len=PATH_MAX)   :: wave_fetch_file=''
   logical                   :: have_wave_fetch=.false.
   logical                   :: have_wave_cap=.false.
   REALTYPE, dimension(:,:,:), allocatable :: wave_fetch_dir,wave_fetch_cap
   REALTYPE, dimension(:,:),   allocatable :: wave_fetch_exp,wave_fetch_tp
#endif
!
! !PRIVATE DATA MEMBERS:
   integer  :: bio_adv_split=0
   integer  :: bio_adv_hor=1
   integer  :: bio_adv_ver=1
   REALTYPE :: bio_AH=-_ONE_
#ifdef BFM_GOTM
   integer         :: bio_hor_adv=1
   integer         :: bio_ver_adv=1
#ifdef STATIC
   REALTYPE        :: delxu(I2DFIELD),delxv(I2DFIELD)
   REALTYPE        :: delyu(I2DFIELD),delyv(I2DFIELD)
   REALTYPE        :: area_inv(I2DFIELD)
#else
   REALTYPE, dimension(:,:), allocatable :: delxu,delxv
   REALTYPE, dimension(:,:), allocatable :: delyu,delyv
   REALTYPE, dimension(:,:), allocatable :: area_inv
#endif
#endif
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_getm_bio
!
! !INTERFACE:
#ifdef BFM_GOTM
!JM   subroutine init_getm_bio(namlst,input_file)
   subroutine init_getm_bio(namlst,input_file,steps_in_outdelt)
#else
   subroutine init_getm_bio(nml_file)
#endif
!
! !DESCRIPTION:
!  Reads the namelist and makes calls to the init functions of the
!  various model components.
!
! !USES:
   use advection, only: J7
   use advection_3d, only: print_adv_settings_3d
#ifdef BFM_GOTM
   use bio_var,only: pelvar_type
   use mem_Silt,only: wave_method
#ifdef MACROPHYT
   use mem_2DMacroPhyto,ONLY:Init2dMacroPhyto
#endif
#endif
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
#ifdef BFM_GOTM
   integer,intent(IN)                             :: namlst
   character(len=*),intent(IN)                    :: input_file
   integer, intent(in)                 :: steps_in_outdelt         !JM added
#else
   character(len=*), intent(in)   :: nml_file
#endif
!
! !REVISION HISTORY:
!  See the log for the module
!
!  !LOCAL VARIABLES
   integer, parameter        :: unit_bio=63
   integer                   :: rc
   integer                   :: i,j,n
   character(len=PATH_MAX)   :: bio_init_file
   integer                   :: bio_init_format, bio_field_no
#ifdef BFM_GOTM
   integer, parameter                  :: unit_bfm=64
   REALTYPE                            :: h(0:kmax)
   REALTYPE                            :: s(0:kmax)
   REALTYPE                            :: t(0:kmax)
   REALTYPE                            :: rho(0:kmax)
   REALTYPE                            :: nuh(0:kmax)
   REALTYPE                            :: dzero
   character(len=80)                   :: brrbrr,msg
   character(len=80)     :: ben_init_file
   character(len=80)     :: ben_param_file
#endif
#ifdef BFM_GOTM
   namelist /getm_bio_nml/ hotstart_bio,bio_hor_adv,bio_ver_adv, &
                           bio_adv_split,bio_AH
   namelist /getm_bfm_nml/ben_param_file,ben_init_file, &
                          calc_init_bennut_states,no_hortr_silt
   namelist /getm_bfm_wave_nml/wave_fetch_file
#else
   namelist /getm_bio_nml/ bio_init_method, &
                           bio_init_file,bio_init_format,bio_field_no, &
                           bio_adv_split,bio_adv_hor,bio_adv_ver,bio_AH
#endif
!EOP
!-------------------------------------------------------------------------
!BOC
   LEVEL2 'init_getm_bio()'

   call init_bio(NAMLST2,'bio.nml',unit_bio,kmax)

#ifdef BFM_GOTM
   h=30.0D+00
   s=35.0
   t=4.0
   rho=1000
   nuh=0.0D+00
   dzero=_ZERO_
#endif

#ifdef BFM_GOTM
   if (bio_calc) then
     if ( bio_setup /=2 ) then
     call set_env_bio_bfm(kmax,dt,dzero,dzero,h,t,s,rho,nuh,12.0D+00,&
          dzero,dzero,dzero,dzero, dzero,julianday)
!JM bfm2024     call set_env_bio_bfm(kmax,dzero,dzero,dzero,dzero,&
!JM bfm2024          dzero,dzero,dzero,dzero,julianday,1.0D0) !AN
!JM     call init_var_bio(namlst,unit_bfm)
!JM not in bfm2016     call copy_from_gotm_to_bfm( &
!JM                 kmax,dt,h,t,s,rho,nuh)
!JM     call init_var_bfm(namlst,'bio_bfm.nml',unit_bfm,bio_setup)
!JM 2024?    call init_var_bfm(namlst,'bio_bfm.nml',unit_bfm,bio_setup,steps_in_outdelt)
     call init_var_bfm(namlst,'bio_bfm.nml',unit_bfm,bio_setup)
     n_cc3d_out=count(var_ave(stPelStateS:stPelStateE) ) + &
                count( var_ids(stPelDiagS:stPelFluxE)/= 0 )
     if ( n_cc3d_out > 0 )  then
       ! pel.biological fields of diagnos.
       allocate(cc3d_out(I3DFIELD,n_cc3d_out),stat=rc)
       if (rc /= 0) stop 'init_getm_bio: Error allocating memory (cc3d_out)'
       cc3d_out=bio_missing
     endif
     ! pel.biological fields of diagnos.
     allocate(counter_ave(I2DFIELD),stat=rc)
     if (rc /= 0) stop 'init_getm_bio: Error allocating memory (counter_ave)'
       counter_ave=0
       if ( numc > 0 ) then
         ! pel.biological fields of states
         allocate(cc3d(I3DFIELD,numc),stat=rc) 
         if (rc /= 0) stop 'init_getm_bio: Error allocating memory (cc3d)'
         ! pel.biological fields of diagnos.
         allocate(adv3d_courant(numc,I2DFIELD),stat=rc)
         if (rc /= 0) stop 'init_getm_bio:Error allocating memory (counter_ave)'
         adv3d_courant=0
         ! pel.biological fields of diagnos.
         allocate(adv3d_number(numc,I2DFIELD),stat=rc)
         if (rc /= 0) stop 'init_getm_bio: Error allocating memory (counter_ave)'
         adv3d_number=0
       endif
     endif

     if ( bio_setup >=2) then
       if ( numbc > 0 )  then
         ! bent, biological fields of states
         allocate(ccb3d(I2DFIELD,0:1,numbc),stat=rc)
         if (rc /= 0) stop 'init_getm_bio: Error allocating memory (ccb3d)'
         ccb3d=bio_missing
       endif
       n_ccb3d_out=count( var_ave(stBenStateS:stBenStateE)) + &
              count( var_ids(stBenDiagS:stBenFluxE) /= 0)
       ! ben. biological fields of flux.
       if ( n_ccb3d_out > 0 )  then
         allocate(ccb3d_out(I2DFIELD,0:1,n_ccb3d_out),stat=rc)
         if (rc /= 0) stop 'init_getm_bio: Error allocating memory (ccb3d_out)'
         ccb3d_out=bio_missing
       endif
       !control to set temporary phytoplankton group on/off
       allocate(sw_CalcPhyto_2d(I2DFIELD,iiPhytoPlankton),stat=rc)
       if (rc /= 0)stop 'init_getm_bio:Error allocating memory(sw_CalcPhyto_2d)'
       sw_CalcPhyto_2d=1
       !control to set temporary phytoplankton group on/off
       allocate(d3_pelvar_type(numc),stat=rc)
       if (rc /= 0) stop 'init_getm_bio:Error allocating memory (d3_pelvar_type)'
       d3_pelvar_type=pelvar_type
#ifdef INCLUDE_DIAGNOS_PRF
       n_ccb3d_prf= count( var_ids(stPRFDiagS:stPRFDiagE) /= 0)
       ! ben. biological fields of flux.
       if ( n_ccb3d_prf > 0 )  then
         allocate(ccb3d_prf(n_ccb3d_prf,I2DFIELD,0:nprf),stat=rc)
         if (rc /= 0) stop 'init_getm_bio: Error allocating memory (ccb3d_prf)'
         ccb3d_prf=bio_missing
       endif
#endif
     endif

     ben_init_file=''
     LEVEL2 "Reading from "//trim(input_file)
     open(namlst,status='unknown',file=trim(input_file),err=90)
     brrbrr='getm_bio_nml'
     read(namlst,NML=getm_bio_nml,end=91,err=92)

     LEVEL2 "Settings related to 3D biological calculations"
     LEVEL3 'bio_hor_adv=   ',bio_hor_adv
     LEVEL3 'bio_ver_adv=   ',bio_ver_adv
     LEVEL3 'bio_adv_split= ',bio_adv_split
     LEVEL3 'bio_AH=        ',bio_AH

     call fill_diagn_bfm_vars( DIAG_INFO,.true., 0,0, h,1,_ZERO_ );n=0
     ! This "initialization' is only done to check if special arrays have to be allocated.
     LEVEL2 "Set all elements of cc3d_out which contain river outputs on 0 &
            & and NOT on 'bio_missing'..."
     LEVEL2 "In this way all active grid points in which no river ents  &
             & are set on 0"
     c1dimz(1:kmax)=_ZERO_
     do j=jmin,jmax
       do i=imin,imax
         if (az(i,j) .ge. 1 ) &
         call make_river_flux_output(_ZERO_,_ZERO_,_ZERO_,c1dimz,numc,i,j)
       enddo
     end do
     n=diag_end_sections(1,3)
     call make_uv_flux_output(0,0,stPelFluxS,stPelFluxE,.TRUE.,n,i)
!JM for calculation horizontal fluxes, code below needs to be uncommented
    if ( i>=1)  then
      calc_uv_fluxes=i
      LEVEL2 "In for the output flow fields will be generated"
      LEVEL3 "therefor memory for variabels 'cut' and 'cvt' will be allocated"
      if (calc_uv_fluxes ==2) LEVEL3 "In ths run waterflux fields will be calculated"
!      allocate(cut(I3DFIELD),stat=rc)    ! work array
      if (rc /= 0) stop 'init_advection_3d: Error allocating memory (cut)'
!      allocate(cvt(I3DFIELD),stat=rc)    ! work array
      if (rc /= 0) stop 'init_advection_3d: Error allocating memory (cvt)'
    endif

     if ( bio_setup /=2 ) then
       if (hotstart_bio) then
         LEVEL2 "Reading biological fields from hotstart file"
       else
         LEVEL2 "Initialise biological fields from namelist"
         do j=jmin,jmax
           do i=imin,imax
!AN             if (az(i,j) .ge. 1 ) cc3d(1:numc,i,j,:)=cc(1:numc,:)
!JM            if (az(i,j) .ge. 1 ) cc3d(:,i,j,1:numc)=cc(:,1:numc) 
            if (az(i,j) .ge. 1 ) cc3d(i,j,:,1:numc)=cc(:,1:numc) 
           enddo
         end do
       endif
     end if

!    optional fetch fields for the Silt wave scheme; the namelist group may be
!    absent, in which case wave_fetch_file stays empty and nothing is read
     wave_fetch_file=''
     rewind(namlst)
     read(namlst,NML=getm_bfm_wave_nml,iostat=rc)
     if (rc .gt. 0) then
       call getm_error('init_getm_bio', &
            'error reading getm_bfm_wave_nml in '//trim(input_file))
     end if
     rewind(namlst)
     if (len_trim(wave_fetch_file) .gt. 0) then
       call init_wave_fetch()
     else
       LEVEL2 'no wave_fetch_file: Silt wave_method 3 and 4 not available'
     end if
     LEVEL2 'Silt wave_method= ',wave_method
     if (wave_method .ge. 3 .and. .not. have_wave_fetch) then
       call getm_error('init_getm_bio', &
            'Silt.nml wave_method 3/4 needs wave_fetch_file in &getm_bfm_wave_nml of ' &
            //trim(input_file))
     end if
     if (wave_method .eq. 4 .and. .not. have_wave_cap) then
       call getm_error('init_getm_bio', &
            'Silt.nml wave_method 4 needs fetch_cap_* fields in '//trim(wave_fetch_file))
     end if

     if ( bio_setup >=2 ) then
       ben_param_file=''
       ben_init_file=''
       brrbrr='getm_bfm_nml'
       read(namlst,NML=getm_bfm_nml)
       close(namlst)
       write(stderr,nml=getm_bfm_nml)
       ! if output through check_presence is .true. a data file is present!
       call set_2d_grid_parameters(read_poro,file=ben_param_file, &
                          bio_missing=bio_missing,test_presence=read_poro)
       if (hotstart_bio.and.(ben_init_file.ne.'')) then
         LEVEL2 "Reading benthic biological fields from hotstart file"
       else
         LEVEL2 "Initialise benthic biological fields from namelist"
         do j=jmin,jmax
           do i=imin,imax
            ! if (az(i,j) .eq. 1 ) ccb3d(1:numbc,i,j,:)=ccb(1:numbc,:)
            !JM above compiles but crashes            
             ! if (az(i,j) .eq. 1 ) ccb3d(i,j,:,1:numbc)=reshape(ccb(1:numbc,:), (/numbc,2/), order=(/2,1/))
             ! AN I think the reshape is not needed anynore
             if (az(i,j) .eq. 1 ) ccb3d(i,j,:,1:numbc)=ccb(:,1:numbc)
           enddo
         end do
         call init_2d_grid(ben_init_file)
         do j=jmin,jmax
           do i=imin,imax
             n=count(ccb3d(i,j,1,1:numbc)<  _ZERO_) 
             if (az(i,j) .ge. 1.and. n> 0 ) then
               LEVEL3 'missing intial benthos data  gridpoint i,j=',i+ioff,j+joff
!AN            ccb3d(1:numbc,i,j,:)=ccb(1:numbc,:)
               ccb3d(i,j,:,1:numbc)=ccb(:,1:numbc) !AN
               LEVEL3 'Initial values defined in bio_bfm.inp are used'
             endif
           enddo
         end do
       endif
     end if

#ifndef STATIC
     allocate(delxu(I2DFIELD),stat=rc)
     if (rc /= 0) stop 'init_getm_bio: Error allocating memory (delxu)'

     allocate(delxv(I2DFIELD),stat=rc)
     if (rc /= 0) stop 'init_getm_bio: Error allocating memory (delxv)'

     allocate(delyu(I2DFIELD),stat=rc)
     if (rc /= 0) stop 'init_getm_bio: Error allocating memory (delyu)'

     allocate(delyv(I2DFIELD),stat=rc)
     if (rc /= 0) stop 'init_getm_bio: Error allocating memory (delyv)'

     allocate(area_inv(I2DFIELD),stat=rc)
     if (rc /= 0) stop 'init_getm_bio: Error allocating memory (area_inv)'

#endif
#if defined(SPHERICAL) || defined(CURVILINEAR)
     delxu=dxu
     delxv=dxv
     delyu=dyu
     delyv=dyv
     area_inv=arcd1
#else
     delxu=dx
     delxv=dx
     delyu=dy
     delyv=dy
     area_inv=ard1
#endif

     do n=1,numc
       call update_3d_halo(cc3d(:,:,:,n),cc3d(:,:,:,n),az, & 
                        imin,jmin,imax,jmax,kmax,D_TAG)
       call wait_halo(D_TAG)
     enddo

   end if
   call set_d3_model_flag_for_gotm(.TRUE.)
#ifdef MACROPHYT
   call Init2dMacroPhyto
#endif
#else
   if (bio_calc) then

!AN      call init_var_bio

      allocate(cc3d(I3DFIELD,numc),stat=rc) 
      if (rc /= 0) stop 'init_getm_bio: Error allocating memory (cc3d)'
      cc3d = _ZERO_
!      cc3d(10,:,:,:) = 0.0001
      open(NAMLST2,status='unknown',file=trim(nml_file))
      read(NAMLST2,NML=getm_bio_nml)
      close(NAMLST2)

      LEVEL2 'Advection of biological fields'
      if (bio_adv_hor .eq. J7) stop 'init_bio: J7 not implemented yet'
      call print_adv_settings_3d(bio_adv_split,bio_adv_hor,bio_adv_ver,bio_AH)

      select case (bio_init_method)
         case(0)
            LEVEL3 'getting initial bio fields from hotstart file'
         case(1)
            LEVEL3 "initial biological fields from namelist - bio_<model>.nml"
            do j=jmin,jmax
               do i=imin,imax
                  if (az(i,j) .ge. 1 ) then
!JM                     cc3d(:,i,j,:)=cc
                     cc3d(i,j,:,:)=cc
                  end if
               end do
            end do
         case(2)
            LEVEL3 'reading initial bio-fields from ',trim(bio_init_file)
            do n=1,numc
               LEVEL4 'inquiring ',trim(var_names(n))
               call get_field(bio_init_file,trim(var_names(n)),bio_field_no, &
                              cc3d(:,:,:,n))
            end do
         case default
            FATAL 'Not valid bio_init_method specified'
            stop 'init_getm_bio()'
      end select

      do n=1,numc
         call update_3d_halo(cc3d(:,:,:,n),cc3d(:,:,:,n),az, &
                             imin,jmin,imax,jmax,kmax,D_TAG)
         call wait_halo(D_TAG)
      end do

   end if
#endif
   return
#ifdef BFM_GOTM
90 msg='can not open ' ;goto 100
91 msg= 'EOF';          goto 100
92 msg='error reading'
100 brrbrr=TRIM(input_file)// ' ('//trim(brrbrr)//' )'
    call getm_error( 'init_getm_bio',trim(msg)//TRIM(brrbrr))
#endif
   end subroutine init_getm_bio
!EOC

#ifdef BFM_GOTM
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_wave_fetch
!
! !INTERFACE:
   subroutine init_wave_fetch()
!
! !DESCRIPTION:
!  Reads the effective-fetch fields used by the BFM Silt wave scheme
!  (Silt.nml: wave_method 3 or 4) from wave_fetch_file. The file has the
!  dimensions of the bathymetry file and holds
!    fetch_dir_000 ... fetch_dir_350   directional fetch, bearing wind comes from
!    fetch_exp, fetch_tp               direction-mean exposure and period fetch
!    fetch_cap_000 ... fetch_cap_350   optional, for wave_method 4
!  Bearings are geographic; do_getm_bio passes convc so that BFM can convert
!  the grid-relative wind direction. Create the file with
!  DWS_CLI/tools/wave_fetch_grid.py.
!
! !USES:
   use domain, only: ilg,ihg,jlg,jhg,ill,ihl,jll,jhl
   use mem_Silt, only: NBEAR_FETCH
   IMPLICIT NONE
!
   interface
      subroutine get_2d_field(fn,varname,il,ih,jl,jh,break_on_missing,f)
         character(len=*),intent(in)   :: fn,varname
         integer, intent(in)           :: il,ih,jl,jh
         logical, intent(in)           :: break_on_missing
         REALTYPE, intent(out)         :: f(:,:)
      end subroutine get_2d_field
   end interface
!
! !LOCAL VARIABLES:
   integer                   :: rc,k,i,j,nbad
   character(len=16)         :: vname
   REALTYPE                  :: fmin,fmax
!EOP
!-------------------------------------------------------------------------
!BOC
   LEVEL2 'init_wave_fetch(): reading ',trim(wave_fetch_file)

   allocate(wave_fetch_dir(NBEAR_FETCH,E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_wave_fetch: Error allocating memory (wave_fetch_dir)'
   allocate(wave_fetch_cap(NBEAR_FETCH,E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_wave_fetch: Error allocating memory (wave_fetch_cap)'
   allocate(wave_fetch_exp(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_wave_fetch: Error allocating memory (wave_fetch_exp)'
   allocate(wave_fetch_tp(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_wave_fetch: Error allocating memory (wave_fetch_tp)'
   wave_fetch_dir=-_ONE_
   wave_fetch_cap=-_ONE_
   wave_fetch_exp=-_ONE_
   wave_fetch_tp=-_ONE_

   do k=1,NBEAR_FETCH
      write(vname,'(A,I3.3)') 'fetch_dir_',(k-1)*(360/NBEAR_FETCH)
      call get_2d_field(trim(wave_fetch_file),trim(vname),ilg,ihg,jlg,jhg,.true., &
                        wave_fetch_dir(k,ill:ihl,jll:jhl))
   end do
   call get_2d_field(trim(wave_fetch_file),'fetch_exp',ilg,ihg,jlg,jhg,.true., &
                     wave_fetch_exp(ill:ihl,jll:jhl))
   call get_2d_field(trim(wave_fetch_file),'fetch_tp',ilg,ihg,jlg,jhg,.true., &
                     wave_fetch_tp(ill:ihl,jll:jhl))
   do k=1,NBEAR_FETCH
      write(vname,'(A,I3.3)') 'fetch_cap_',(k-1)*(360/NBEAR_FETCH)
      call get_2d_field(trim(wave_fetch_file),trim(vname),ilg,ihg,jlg,jhg,.false., &
                        wave_fetch_cap(k,ill:ihl,jll:jhl))
   end do

!  every wet column needs valid values
   nbad=0
   do j=jmin,jmax
      do i=imin,imax
         if (az(i,j) .ge. 1) then
            if (any(.not.(wave_fetch_dir(:,i,j) .ge. _ZERO_)) .or. &
                .not.(wave_fetch_exp(i,j) .ge. _ZERO_) .or.       &
                .not.(wave_fetch_tp(i,j) .ge. _ZERO_)) nbad=nbad+1
         end if
      end do
   end do
   if (nbad .gt. 0) then
      LEVEL3 'wet points without valid fetch values: ',nbad
      call getm_error('init_wave_fetch', &
           'missing or negative fetch in '//trim(wave_fetch_file))
   end if

   have_wave_cap=.true.
   do j=jmin,jmax
      do i=imin,imax
         if (az(i,j) .ge. 1) then
            if (any(.not.(wave_fetch_cap(:,i,j) .ge. _ZERO_))) have_wave_cap=.false.
         end if
      end do
   end do
   have_wave_fetch=.true.

   fmin=minval(wave_fetch_exp(imin:imax,jmin:jmax),mask=az(imin:imax,jmin:jmax).ge.1)
   fmax=maxval(wave_fetch_exp(imin:imax,jmin:jmax),mask=az(imin:imax,jmin:jmax).ge.1)
   LEVEL3 'fetch_exp range (m): ',fmin,fmax
   if (have_wave_cap) then
      LEVEL3 'fetch_cap fields found: wave_method 4 available'
   else
      LEVEL3 'no fetch_cap fields: wave_method 4 not available'
   end if

   return
   end subroutine init_wave_fetch
!EOC
#endif

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  do_getm_bio()
!
! !INTERFACE:
#ifndef BFM_GOTM
   subroutine do_getm_bio(dt,write_3d)
#else
  subroutine do_getm_bio(dt,write_3d,llinit)
#endif
!
! !DESCRIPTION:
!
! !USES:
   use advection_3d, only: do_advection_3d
   use getm_timers, only: tic, toc, TIM_GETM_BIO, TIM_ADVECTBIO
#ifdef BFM_GOTM
   use bio,only:getm_stderr_control,warning_level,get_bio_updates !JM, &
!JM                struct_test_3d,struct_test_2d,struct_test_part
!JM   use bio_bfm,only:test_structure,test_model_states
   use global_mem,only:SILTTRANSPORT
   use variables_bio_3d,only:counter_reset
   use controlled_messages,only:controlled_output
   use meteo, only:swr,u10,v10,dl                                     !BFM
   use domain,only:dry_z,H
   use domain,only:convc
   use mem_Silt,only:set_wave_column
#ifdef MACROPHYT
   use mem_2DMacroPhyto,ONLY:Do2dMacroPhyto, Fin2dMacroPhyto
#endif
#endif

   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: dt
#ifdef BFM_GOTM
   logical, intent(in)                 :: write_3d
   logical, intent(in),optional        :: llinit
#endif
!
! !REVISION HISTORY:
!  See the log for the module
!
! !LOCAL VARIABLES:
   integer         :: n
   integer         :: i,j,k
   REALTYPE,dimension(I3DFIELD) :: fadv3d
   REALTYPE        :: h1d(0:kmax),T1d(0:kmax),S1d(0:kmax),rho1d(0:kmax)
   REALTYPE        :: nuh1d(0:kmax),rad1d(0:kmax),light1d(0:kmax)
   REALTYPE        :: bioshade1d(0:kmax)
   REALTYPE        :: wind_speed,I_0
   integer         :: ltimes !AN
#ifdef BFM_GOTM
   integer         :: iout,shiftcounter
   REALTYPE        :: wind,uv_b,bath_dep
   REALTYPE        :: r,remember
   REALTYPE,dimension(I3DFIELD) :: cc3d_halo
   integer,parameter   :: notransport=0
   CHARACTER(LEN=100) :: sub,outp
   CHARACTER(LEN=100) :: msg=""
   LOGICAL         :: error_flag,warning_flag
   REALTYPE,dimension(:,:,:),pointer,contiguous :: p_cut,p_cvt
#endif
!EOP
!-----------------------------------------------------------------------
!BOC

   call tic(TIM_GETM_BIO)

#ifdef BFM_GOTM
!   LEVEL2 'getm_bio: do advection '
!LEVEL1 'R2c',cc3d(2,9,:,42)
!LEVEL1 'N1p',cc3d(4,6,20,3)
!LEVEL1 'N3n',cc3d(4,6,20,4)
!LEVEL1 'uu',uu(4,6,20)
    p_cut => cut
    p_cvt => cvt

!  then we do the advection of the biological variables
   if ( .not.bio_calc) return
   if ( bio_setup /= 2 .and.(.not.present(llinit)) )  then
     j=1
     ! Check and if needed reinitialize output variables in which uv-fluxes are stored...
!    k=diag_end_sections(2,3)
!    call make_nettransport_flux_output(0,dt,0,stBenFluxS,stBenFluxE,.TRUE.,k,j)
     k=diag_end_sections(1,3)
     call make_uv_flux_output(0,0,stPelFluxS,stPelFluxE,.TRUE.,k,j)
     k=diag_end_sections(1,5)
     j=1;if ( calc_uv_fluxes==2) j=0
     call tic(TIM_ADVECTBIO)
     do n=j,numc
       i=-1
       if (n==0 ) then
         i=n;ffp=1.0D+00
       elseif ( d3_pelvar_type(n)==SILTTRANSPORT) then
         i=n;ffp(:,:,:) = cc3d(:,:,:,n);if (no_hortr_silt)i=0
       elseif ( d3_pelvar_type(n)/=notransport) then
         i=n;ffp(:,:,:) = cc3d(:,:,:,n)
       endif
!if ( n.eq.4 .or. n.eq.3) then
!LEVEL1 'calling do_advection, n,i',n,i
!LEVEL1 'dt,bio_adv_split,bio_hor_adv,bio_ver_adv,bio_AH,H_TAG',dt,bio_adv_split,bio_hor_adv,bio_ver_adv,bio_AH,H_TAG
!LEVEL1 'ffp',ffp(4,6,:)
!LEVEL1 'cc3d',cc3d(4,6,:,n)
!LEVEL1 'uu',uu(4,6,:)
!LEVEL1 'vv',vv(4,6,:)
!LEVEL1 'ww',ww(4,6,:)
!LEVEL1 'hun',hun(4,6,:)
!LEVEL1 'hvn',hvn(4,6,:)
!LEVEL1 'ho',ho(4,6,:)
!LEVEL1 'hn',hn(4,6,:)
!if ( n.eq.4 ) then
!stop
!endif
!if ( n.eq.42 ) then
!LEVEL1 'calling do_advection, n,i',n,i
!LEVEL1 'dt,bio_adv_split,bio_hor_adv,bio_ver_adv,bio_AH,H_TAG',dt,bio_adv_split,bio_hor_adv,bio_ver_adv,bio_AH,H_TAG
!LEVEL1 'ffp',ffp(2,9,:)
!LEVEL1 'cc3d',cc3d(2,9,:,n)
!LEVEL1 'uu',uu(2,9,:)
!LEVEL1 'vv',vv(2,9,:)
!LEVEL1 'ww',ww(2,9,:)
!LEVEL1 'hun',hun(2,9,:)
!LEVEL1 'hvn',hvn(2,9,:)
!LEVEL1 'ho',ho(2,9,:)
!LEVEL1 'hn',hn(2,9,:)
!LEVEL1 'ffphor'
!LEVEL0 ffp(1:12,1:10,25)
!LEVEL1 'uuhor'
!LEVEL1 uu(1:12,1:10,25)
!LEVEL1 'vvhor'
!LEVEL1 vv(1:12,1:10,25)
!LEVEL1 'wwhor'
!LEVEL1 ww(1:12,1:10,25)
!LEVEL1 'hunhor'
!LEVEL1 hun(1:12,1:10,25)
!LEVEL1 'hvnhor'
!LEVEL1 hvn(1:12,1:10,25)
!LEVEL1 'hohor'
!LEVEL1 ho(1:12,1:10,25)
!LEVEL1 'hnhor'
!LEVEL1 hn(1:12,1:10,25)
!LEVEL1 'faked dt'
!        call do_advection_3d(0.1D0,ffp,uu,vv,ww,hun,hvn,ho,hn, &
!              bio_adv_split,bio_hor_adv,bio_ver_adv,_ZERO_,H_TAG &
!              )
!LEVEL1 'after advection, n',n
!LEVEL1 'ffp',ffp(1:12,1:10,25)
!stop
!endif

       if ( i>=0) &
!JM        call do_advection_3d(dt,ffp,uu,vv,ww,hun,hvn,ho,hn, &
!JM              delxu,delxv,delyu,delyv,area_inv,az,au,av, &
!JM              bio_hor_adv,bio_ver_adv,bio_adv_split,bio_AH)
!        call do_advection_3d(dt,ffp,uu,vv,ww,hun,hvn,ho,hn, &
!              bio_adv_split,bio_hor_adv,bio_ver_adv,bio_AH,H_TAG &
!              )
!JM for horizontal flux calculations: add cut,cvt as output variables to call below
! currently only works for tracking
!        call do_advection_3d(dt,ffp,uu,vv,ww,hun,hvn,ho,hn, &
!              bio_adv_split,bio_hor_adv,bio_ver_adv,_ZERO_,H_TAG &
!              )
        call do_advection_3d(dt,ffp,uu,vv,ww,hun,hvn,ho,hn, &
              bio_adv_split,bio_hor_adv,bio_ver_adv,_ZERO_,H_TAG, &
              ffluxu=p_cut,ffluxv=p_cvt)
!if ( n.eq.4 .or. n.eq.3) then
!LEVEL1 'after advection, n',n
!LEVEL1 'ffp',ffp(4,6,:)
!endif
!if ( n.eq.42 ) then
!LEVEL1 'after advection, n',n
!LEVEL1 'ffp',ffp(2,9,:)
!endif
       if ( i>=1 ) then
         cc3d(:,:,:,n) =ffp(:,:,:)
       endif
       k=diag_end_sections(1,3)
       call make_uv_flux_output(1,n,stPelFluxS,stPelFluxE,.TRUE.,k,j)
     end do
   endif
   call toc(TIM_ADVECTBIO)
   call tic(TIM_GETM_BIO)

!LEVEL1 'after advection'
!LEVEL1 'R2c',cc3d(2,9,:,42)
!stop
!LEVEL1 'N1p',cc3d(4,6,20,3)
!LEVEL1 'N3n',cc3d(4,6,20,4)

!  First we do all the vertical processes
!   LEVEL3 'getm_bio: do column '
   shiftcounter=0
   do j=jmin,jmax
      do i=imin,imax
         if (az(i,j) .eq. 1 ) then
           I_0=swr(i,j)
           h1d=hn(i,j,:)
           T1d=T(i,j,:)
           S1d=S(i,j,:)
           Rho1D=rho(i,j,:)
           nuh1d=nuh(i,j,:)
           light1d=light(i,j,:)
           bath_dep=H(i,j)
!JM           call test_structure(struct_test_3d,struct_test_2d,struct_test_part)
           getm_stderr_control=i.eq.1.and.j.eq.1
           if ( bio_setup /= 2 ) then
!JM             cc(:,:)=cc3d(:,i,j,:)
             cc(:,:)=cc3d(i,j,:,:)
             if ( warning_level >2 ) then
!JM               adv1d_courant=adv3d_courant(:,i,j)
!JM               adv1d_number=adv3d_number(:,i,j)
               adv1d_courant=adv3d_courant(i,j,:)
               adv1d_number=adv3d_number(i,j,:)
             endif
           endif
           if ( bio_setup >= 2 ) then
                call set_2d_grid_parameters(read_poro,igrid=i,jgrid=j)
!JM                ccb(:,:)=ccb3d(:,i,j,:)
                ccb(:,:)=ccb3d(i,j,:,:)
           endif
           call set_env_bio_bfm(kmax,dt,bath_dep,sqrt(taub(i,j)), &
             h1d,T1d,S1d,Rho1d,nuh1d,dl(i,j),&
             u10(i,j),v10(i,j),uu(i,j,1)/hun(i,j,1),vv(i,j,1)/hvn(i,j,1), &
                                  I_0,julianday,dry_z(i,j))
!          fetch fields and grid convergence for the Silt wave scheme
           if (have_wave_fetch) then
             if (have_wave_cap) then
               call set_wave_column(convc(i,j),wave_fetch_dir(:,i,j), &
                      wave_fetch_exp(i,j),wave_fetch_tp(i,j),         &
                      fetch_cap=wave_fetch_cap(:,i,j))
             else
               call set_wave_column(convc(i,j),wave_fetch_dir(:,i,j), &
                      wave_fetch_exp(i,j),wave_fetch_tp(i,j))
             end if
           end if
!AN           call set_env_bio_bfm(bath_dep,sqrt(taub(i,j)),dl(i,j),u10(i,j),v10(i,j),&
!AN                                I_0,julianday,dry_z(i,j) &
!AN                                )
!JM bfm2024           call set_env_bio_bfm(kmax,bath_dep,sqrt(taub(i,j)),dl(i,j),u10(i,j),v10(i,j),&
!JM bfm2024                                uu(i,j,1)/hun(i,j,1),vv(i,j,1)/hvn(i,j,1),I_0,julianday,dry_z(i,j) &
!JM bfm2024                                )
           sw_CalcPhyto(:,1)=sw_CalcPhyto_2d(i,j,:)
#ifdef MACROPHYT
           call Do2dMacroPhyto(julianday,i,j,1)
#endif
!LEVEL1 'i,j,az ',i,j,az(i,j)
!JM old call style          call do_bio(dt,kmax,write_3d,shiftcounter, &
!JM                                             counter_reset, &
!JM                                             h1d,nuh1d)
!LEVEL1 'getm_bio, calling do_bio, i,j',i,j
!LEVEL1 'az',az(i,j)
!LEVEL1 'cc',cc(:,3)
!LEVEL1 'ccb',ccb(:,2)
!LEVEL1 'Dcm',ccb(:,84)
!LEVEL1 'N1p',cc3d(i,j,20,3)
!LEVEL1 'N3n',cc3d(i,j,20,4)
!stop
           call do_bio(write_3d)
!LEVEL1 'N1p',cc3d(i,j,20,3)
!LEVEL1 'N3n',cc3d(i,j,20,4)
!stop
           if ( bio_setup /= 2 ) then
!JM             cc3d(:,i,j,:)=cc(:,:)
             cc3d(i,j,:,:)=cc(:,:)
             sw_CalcPhyto_2d(i,j,:)= sw_CalcPhyto(:,1)
             call get_warning_for_getm(warning_flag )
             if (write_3d.and.warning_level>2) then
               k=0;r=counter_ave(i,j)
               do n=1,numc
                 if (adv1d_number(n)/r > _ONE_) then
                   write(outp, &
                   '(A,'':adv_center:MaxCourantN.='',G12.4,'' iter.='',F5.1)')&
                   trim(var_names(n)),adv1d_courant(n)/r,adv1d_number(n)/r
                     k=len_trim(outp);LEVEL1 outp(1:k)
                 endif
               enddo
               if ( k>0) then
                 warning_flag=.false.
#ifdef GETM_PARALLEL
                 write(outp,'(''i,j,='',I4,''('',I2,'') '',I4,''('',I2,'') &
                   &Depth='',F10.3)') i+ioff,i,j+joff,j,sum(h1d)
#else
                 write(outp,'(''i,j,='',I4,''('',I3,'') '',I4,''('',I3,'') &
                   &Depth='',F10.3)') i+ioff,i,j+joff,j,sum(h1d)
#endif
                  k=len_trim(outp);LEVEL3 outp(1:k)
               endif
               adv1d_number=_ZERO_
               adv1d_courant=_ZERO_
             endif
             if (warning_level >2) then
!JM               adv3d_courant(:,i,j)=adv1d_courant
!JM               adv3d_number(:,i,j)=adv1d_number
               adv3d_courant(i,j,:)=adv1d_courant
               adv3d_number(i,j,:)=adv1d_number
             endif
           endif
!           if (warning_flag ) then
!             r=sum(h1d)
!             k=len_trim(msg);LEVEL3 msg(1:k)
!LEVEL1 'getm_bio: Warning flag set for grid cell i,j:',i,j
!LEVEL1 'Stopping run'
!stop
!           endif

!JM           if ( bio_setup >= 2 ) ccb3d(:,i,j,:)=ccb(:,:)
           if ( bio_setup >= 2 ) ccb3d(i,j,:,:)=ccb(:,:)
           call fill_diagn_bfm_vars( DIAG_ADD,write_3d, i,j, h1d,kmax,dt)
           call output_gotm_error( error_flag, sub, msg)
           if ( error_flag ) then
             k=len_trim(outp)
             call getm_error( sub, msg(1:k))
           endif
           call get_bio_updates(kmax,bioshade1d)
           light(i,j,:)=bioshade1d
#ifdef MACROPHYT
           call Fin2dMacroPhyto(i,j,1)
#endif
!          call test_model_states(11,0,numbc,1,ccb)
         end if
      end do
   end do
!LEVEL1 'getm-bio after loops'
!stop
   ltimes=0;if (write_3d) ltimes=1 !AN
!   call controlled_output(ltimes,'') !AN
   if ( present(llinit)) return
   if (shiftcounter > 0) then
     if (warning_level >2) then
        STDERR 'number of shiftings of masses between vertical ', &
         'layers < 0.001%=',shiftcounter
     else
        STDERR 'number of shiftings of masses between vertical ', &
         'layers=',shiftcounter
     endif
   endif
!LEVEL1 'bio_setup,D_TAG,numc: ',bio_setup,D_TAG,numc
!stop
   if ( bio_setup /= 2 )  then
     do n=1,numc
!LEVEL1 'numc, n:',numc,n
        cc3d_halo=cc3d(:,:,:,n)
        call update_3d_halo(cc3d_halo,cc3d_halo,az, &
                         imin,jmin,imax,jmax,kmax,D_TAG)
        call wait_halo(D_TAG)
        cc3d(:,:,:,n)=cc3d_halo
     end do
   endif
!LEVEL1 'after update halo'
!stop
!  if ( bio_setup >= 2 )  then
!    do n=1,numbc
!      ffb(:,:)=ccb3d(n,:,:,1)
!      call update_2d_halo(ffb,ffb,az,imin,jmin,imax,jmax,D_TAG)
!      call wait_halo(D_TAG)
!      ccb3d(n,:,:,1)=ffb(:,:)
!    enddo
!  endif
   call toc(TIM_GETM_BIO)
!LEVEL1 'after bio_setup'
!stop
#else
!not BFM
!  First we do all the vertical processes
   do j=jmin,jmax
      do i=imin,imax
         if (az(i,j) .ge. 1 ) then
#ifdef GOTM_V3
            I_0=swr(i,j)
            h1d=hn(i,j,:)
            T1d=T(i,j,:)
            S1d=S(i,j,:)
            nuh1d=nuh(i,j,:)
            light1d=light(i,j,:)
!JM            cc=cc3d(:,i,j,:)
            cc=cc3d(i,j,:,:)
            call do_bio(kmax,I_0,dt,h1d,T1d,S1d,nuh1d,light1d,bioshade1d)
!JM            cc3d(:,i,j,:)=cc
            cc3d(i,j,:,:)=cc
            light(i,j,:)=bioshade1d
#else
            h1d=hn(i,j,:)
            T1d=T(i,j,:)
            S1d=S(i,j,:)
            rho1d=rho(i,j,:)
            nuh1d=nuh(i,j,:)
            rad1d=rad(i,j,:)
            if (allocated(u10) .and. allocated(v10)) then
               wind_speed=sqrt(u10(i,j)*u10(i,j)+v10(i,j)*v10(i,j))
            else
               wind_speed=_ZERO_
            end if
            if (allocated(swr)) then
               I_0=swr(i,j)
            else
               I_0=_ZERO_
            end if
            light1d=light(i,j,:)
!JM            cc=cc3d(:,i,j,:)
            cc=cc3d(i,j,:,:)
            call set_env_bio(kmax,dt,-D(i,j),sqrt(taub(i,j)), &
                             h1d,T1d,S1d,rho1d,nuh1d,rad1d,   &
                             wind_speed,I_0,secondsofday,w_adv_ctr_=0)
            call do_bio()
!JM            cc3d(:,i,j,:)=cc
            cc3d(i,j,:,:)=cc
!            light(i,j,:)=bioshade1d
#endif
         end if
      end do
   end do

!  then we do the advection of the biological variables
   call tic(TIM_ADVECTBIO)
   do n=1,numc

#if 1
      fadv3d = cc3d(:,:,:,n)
      call update_3d_halo(fadv3d,fadv3d,az, &
                          imin,jmin,imax,jmax,kmax,D_TAG)
      call wait_halo(D_TAG)

!     KK-TODO: bio_AH_method + include bio_AH_method=1 into advection

      call do_advection_3d(dt,fadv3d,uu,vv,ww,hun,hvn,ho,hn,                   &
                           bio_adv_split,bio_adv_hor,bio_adv_ver,bio_AH,H_TAG)

!      if (bio_AH_method .gt. 1) then
!         call update_3d_halo(fadv3d,fadv3d,az,imin,jmin,imax,jmax,kmax,D_TAG)
!         call wait_halo(D_TAG)
!         call tracer_diffusion(fadv3d,hn,bio_AH_method,bio_AH_const,bio_AH_Prt,bio_AH_stirr_const)
!      end if

      cc3d(:,:,:,n) = fadv3d
#else
      call update_3d_halo(cc3d(:,:,:,n),cc3d(:,:,:,n),az, &
                          imin,jmin,imax,jmax,kmax,D_TAG)
      call wait_halo(D_TAG)

!     KK-TODO: bio_AH_method + include bio_AH_method=1 into advection

      call do_advection_3d(dt,cc3d(:,:,:,n),uu,vv,ww,hun,hvn,ho,hn,            &
                           bio_adv_split,bio_adv_hor,bio_adv_ver,bio_AH,H_TAG)

!      if (bio_AH_method .gt. 1) then
!         call update_3d_halo(cc3d(n,:,:,:),cc3d(n,:,:,:),az,imin,jmin,imax,jmax,kmax,D_TAG)
!         call wait_halo(D_TAG)
!         call tracer_diffusion(cc3d(n,:,:,:),hn,bio_AH_method,bio_AH_const,bio_AH_Prt,bio_AH_stirr_const)
!      end if
#endif
   end do

   call toc(TIM_ADVECTBIO)
   call toc(TIM_GETM_BIO)
#endif
!LEVEL1 'end do_getm_bio'
!stop
   return
   end subroutine do_getm_bio
!EOC
#ifdef BFM_GOTM
!--------------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: init_calc_getm_bio
!
! !INTERFACE:
     subroutine init_calc_getm_bio
!
! DESCRIPTION
!    Routine which call routines specific valid for biomodels
!    which calc initial conditions from other initial conditions.
!
! !USES:
     use bio_var,   only: bio_model
!
!
!
!EOP
!-------------------------------------------------------------------------
!BOC
      if ( .not.bio_calc) return
      select case (bio_model)
      case (6)  ! The BFM (ERSEM) model
          LEVEL3 "init_calc_getm_bio: bio_model=",bio_model
          call getm_bfm_bennut_calc_initial(calc_init_bennut_states,read_poro,hotstart_bio)
      case  default
         continue
      end select

     end subroutine init_calc_getm_bio
!EOC
#endif

!-----------------------------------------------------------------------

   end module getm_bio

!-----------------------------------------------------------------------
! Copyright (C) 2007 - Karsten Bolding and Hans Burchard               !
!-----------------------------------------------------------------------
