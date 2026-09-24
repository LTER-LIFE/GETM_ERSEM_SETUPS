!$Id$
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: wave --- wave calculation \label{sec:wave}
!
! !INTERFACE:
   module wave
!
! !DESCRIPTION:
!  Wave calculation
!
!  Contains main switch box
!  to various wave calculation methods
!
!  wave_mode 1: simple JONSWAP equilibrium wave method, as originally compiled
!  wave_mode 2: effective-fetch scheme (Breugem & Holthuijsen local sea,
!               propagated sea on the exposure fetch, SPM peak period)
!  Parameters and column fetch values of mode 2 are in mem_Silt
!  (Silt.nml, set_wave_column).
!
!
!
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! !USES:

!
!  default: all is private.
!   private
!
! !PUBLIC MEMBER FUNCTIONS:
!   public do_wave
   public
!
! !PRIVATE MEMBER FUNCTIONS:

!
! !REVISION HISTORY:
!  Original author(s): Johan Van Der Molen 
!
! Nov 2013: taken from  gotm-3.3.2_withspm_may2011 without modifications
!
!  $Log$
!
! !LOCAL VARIABLES:

!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Manages SPM equations
!
! !INTERFACE:
   subroutine do_wave(wave_mode,Hhs,Ttz,Phiw,wind,windx,windy,depth)
!
! !DESCRIPTION:
!  This routine is the main switch box that selects the method
!  to determine spm concentrations
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: wave_mode
   REALTYPE, intent(in)                 :: wind,windx,windy,depth
   REALTYPE, intent(out)                :: Hhs,Ttz,Phiw
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard 
!
! !LOCAL VARIABLES

!EOP
!-----------------------------------------------------------------------
!BOC

! establish wave conditions

   select case (wave_mode)
   case (1)   ! simple JONSWAP equilibrium model
     call jonswap(Hhs,Ttz,Phiw,wind,windx,windy,depth)
   case (2)   ! effective-fetch scheme
     call effective_fetch_waves(Hhs,Ttz,Phiw,wind,windx,windy,depth)
   case default
     write(*,*) 'do_wave: unknown wave_mode ',wave_mode
     stop 'do_wave'
   end select

   return
   end subroutine do_wave
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Calculate wave conditions using JONWAP method
!
! !INTERFACE:
   subroutine jonswap(Hhs,Ttz,Phiw,wind,windx,windy,depth)
!
! !DESCRIPTION:
! Calculates wave characteristics assuming equilibrium with wind and JONSWAP spectrum
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                 :: wind,windx,windy,depth
   REALTYPE, intent(out)                :: Hhs,Ttz,Phiw
!

! !LOCAL VARIABLES
   REALTYPE, parameter  :: g=9.81 
   REALTYPE             :: U,F,Fstar,Td,Tdm
   REALTYPE, parameter :: pi=3.141592654, d360=360.0
!JM Tuning params: could come in through namelist in future
   REALTYPE, parameter  :: Hs_min=0.0, Tz_min=1.0, Tdmax=3*3600 
!EOP
!-----------------------------------------------------------------------
!BOC
   
   if (wind.gt.1E-3) then
     ! assume reduced duration in shallow water
     if (depth.lt.10.0) then 
       Tdm=((depth/10.0)**2)*Tdmax
     else
       Tdm=Tdmax
     endif
     ! simple wave generation: vRijn p 331
     F=wind*Tdm                        ! fetch
     U=0.7*(wind**1.2)
     Fstar=g*F/(U**2)
     Td=68.8*(Fstar**0.67)*U/g
     Hhs=max(Hs_min,min(0.243*(U**2)/g,0.0016*sqrt(Fstar)*(U**2)/g))
     Hhs=min(Hhs,0.4*depth)
     Ttz=max(Tz_min,min(8.14*U/g,0.286*(Fstar**0.33)*U/g))

     Phiw=atan2(windy,windx)
     if (Phiw .lt. _ZERO_) then       ! degrees wrt North
        Phiw=180.*(pi/2-Phiw)/pi
     else
        Phiw=180.*(2*pi-(Phiw-pi/2))/pi
     end if
     Phiw=mod(Phiw,d360)             ! reduce to <360

   else
     Hhs=0.0
     Ttz=Tz_min
     Phiw=0.0
   endif

   return
   end subroutine jonswap
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Effective-fetch wave scheme (wave_mode 2)
!
! !INTERFACE:
   subroutine effective_fetch_waves(Hhs,Ttz,Phiw,wind,windx,windy,depth)
!
! !DESCRIPTION:
! Height = min( sqrt(H_local^2 + (wave_alpha_exp*H_prop)^2), wave_gamma*depth )
!   H_local  Breugem & Holthuijsen (2007) with the directional effective fetch
!            at the wind direction and the local depth (coefficients as in
!            GETM src/waves/waves.F90, wind2waveHeight)
!   H_prop   SPM deep-water height on the exposure fetch, limited by
!            wind*wave_Tdmax_exp
! Period  = SPM peak-period curve on the period fetch, limited by
!            wind*wave_Tdmax_tp, no depth reduction. Returned in Ttz as
!            jonswap() does; it is a peak period.
! Calibrated against 24 RWS stations, 2015 (DWS 500 m).
!
! !USES:
   use mem_Silt, only: wave_gamma,wave_alpha_exp,wave_Tdmax_exp, &
                       wave_Tdmax_tp,wave_Tz_min,wave_convc, &
                       wave_fetch_dir,wave_fetch_exp,wave_fetch_tp,wave_fetch_set
   use mem, only: InitializeModel
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                 :: wind,windx,windy,depth
   REALTYPE, intent(out)                :: Hhs,Ttz,Phiw
!
! !LOCAL VARIABLES
   REALTYPE, parameter :: g=9.81D0
!  Breugem & Holthuijsen (2007)
   REALTYPE, parameter :: H8=0.24D0, k1=4.14D-4, m1=0.79D0, k3=0.343D0, m3=1.14D0, pw=0.572D0
   REALTYPE, parameter :: Hs_min=0.0D0
   REALTYPE :: Fd,dstar,fstar,tanhd,Hloc,U,Fe,Hprop,Ft
!EOP
!-----------------------------------------------------------------------
!BOC
   if (.not. wave_fetch_set) then
!    init_var_bfm calls SiltDynamics on a dummy column before GETM has read
!    the fetch file: return calm conditions during that call only
     if (InitializeModel .ne. 0) then
       Hhs=0.0D0
       Ttz=wave_Tz_min
       Phiw=0.0D0
       return
     endif
     write(*,*) 'effective_fetch_waves: wave_method 2 needs wave_fetch_file in getm_bio.inp'
     stop 'effective_fetch_waves'
   endif

   if (wind.gt.1.0D-3) then
     call wind_direction(windx,windy,Phiw)
     Fd=fetch_at_wind(Phiw,wave_convc,wave_fetch_dir)

!    local sea: Breugem & Holthuijsen
     dstar=g*max(depth,1.0D-6)/(wind*wind)
     fstar=g*Fd/(wind*wind)
     tanhd=max(tanh(k3*dstar**m3),1.0D-12)
     Hloc=H8*(tanhd*tanh(k1*fstar**m1/tanhd))**pw*wind*wind/g

!    propagated sea: SPM deep-water curve on the exposure fetch
     U=0.7D0*(wind**1.2D0)
     Fe=min(wind*wave_Tdmax_exp,wave_fetch_exp)
     Hprop=min(0.243D0*(U**2)/g,0.0016D0*sqrt(g*Fe/(U**2))*(U**2)/g)

     Hhs=sqrt(Hloc**2+(wave_alpha_exp*Hprop)**2)
     Hhs=max(Hs_min,min(Hhs,wave_gamma*depth))

!    period: SPM peak-period curve on the period fetch
     Ft=min(wind*wave_Tdmax_tp,wave_fetch_tp)
     Ttz=max(wave_Tz_min,min(8.14D0*U/g,0.286D0*((g*Ft/(U**2))**0.33D0)*U/g))
   else
     Hhs=0.0D0
     Ttz=wave_Tz_min
     Phiw=0.0D0
   endif

   return
   end subroutine effective_fetch_waves
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Wind direction as in jonswap()
!
! !INTERFACE:
   subroutine wind_direction(windx,windy,Phiw)
!
! !DESCRIPTION:
! Bearing (deg, 0-360) the wind blows TOWARD, relative to the frame of
! windx/windy. In GETM that is the model grid: grid north is at compass
! bearing convc.
!
! !USES:
   IMPLICIT NONE
   REALTYPE, intent(in)                 :: windx,windy
   REALTYPE, intent(out)                :: Phiw
   REALTYPE, parameter :: pi=3.141592654, d360=360.0
!EOP
!BOC
   Phiw=atan2(windy,windx)
   if (Phiw .lt. _ZERO_) then       ! degrees wrt North
      Phiw=180.*(pi/2-Phiw)/pi
   else
      Phiw=180.*(2*pi-(Phiw-pi/2))/pi
   end if
   Phiw=mod(Phiw,d360)             ! reduce to <360
   return
   end subroutine wind_direction
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Effective fetch at the wind direction
!
! !INTERFACE:
   REALTYPE function fetch_at_wind(Phiw,convc,table)
!
! !DESCRIPTION:
! Linear interpolation of a fetch table given for geographic bearings the
! wind blows FROM (360/n spacing, starting at 0) at the bearing
! Phiw + 180 + convc, where Phiw is the grid-relative direction the wind
! blows toward and convc the compass bearing of grid north.
!
! !USES:
   IMPLICIT NONE
   REALTYPE, intent(in)                 :: Phiw,convc
   REALTYPE, intent(in)                 :: table(:)
!
! !LOCAL VARIABLES
   REALTYPE, parameter :: d360=360.0
   REALTYPE            :: wfrom,x,f
   integer             :: n,i0,i1
!EOP
!BOC
   n=size(table)
   wfrom=modulo(Phiw+180.0D0+convc,d360)
   x=wfrom/(d360/n)
   i0=int(x)
   f=x-i0
   i0=mod(i0,n)+1
   i1=mod(i0,n)+1
   fetch_at_wind=(1.0D0-f)*table(i0)+f*table(i1)
   return
   end function fetch_at_wind
!EOC

   end module wave

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
