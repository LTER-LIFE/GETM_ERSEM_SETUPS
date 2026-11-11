#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: BenGlobal
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
  subroutine BenGlobalDynamics(mode)
!
! !USES:

  ! For the following Benthic-states fluxes are defined: D6m, D7m, D8m, D9m
  ! The following Benthic-states are used (NOT in fluxes): Y2c, Y5c, Y1c, Y4c
  ! The following Benthic 1-d global boxvars  are used: ETW_Ben
  ! The following global constants are used: RLEN

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:LOGUNIT,RLEN,ZERO
  use mem, ONLY: ruBPc, ruBPn, ruBPp, ruBPs, &
    jnetBPTc,jPIY3c,jZEY3c,jZIY3c,jY3N4n,jY3N1p,jY3RIc, jRIQIc,jRIQIs, &
    jO2Y2o,jmYIc,jP6Y3c,jPIQ6s,jBenFishInput,jCaCO3Y3c, &
    ResetSource_D2_vector,flux_vector,NO_BOXES_XY,  &
    ppBenthicCO2,iiBenthicCO2, iiBenPhyto,ppBenPhyto, &
    ppBenUrea,iiBenUrea, ppBenLabileDetritus,iiBenLabileDetritus, &
    iiBen,iiC,iiN,iiP,iiS,  ppQ6s

  ! The following vector functions are used:MM_vector
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!  
!
! !AUTHORS
!   P.Ruardij
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

   integer,intent(IN)            :: mode
   
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   integer                           :: i
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! In the betnhic processes all respirations and excretions
  ! are added to the rr???? and re??? variables.
  ! There rates are input to the Benthic Nutrient model
  ! first these variables are initialized:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  select case (mode)
     ! Pre processing: this done before calls to Ben*routines

    case(1)
      jnetBPTc(:)=   ZERO  ! mgC/m2 # Total Benthic primary production
      ruBPc(:)  =   ZERO  ! mgC/m2  # Total BenthicBac oxic N uptake
      ruBPn(:)  =   ZERO  ! mmol N/m2  # Total BenthicBac oxic N uptake
      ruBPp(:)  =   ZERO  ! mmol P/m2  # Total BenthicBac oxic P uptake
      ruBPs(:)  =   ZERO  ! mmol Si/m2  # Total BenthicBac oxic Si uptake

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! All nutrient flux rates are initalized here because nutrient are given
      ! back by Filterfeeders and by the nutrient regeration model 
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   

      jPIQ6s(:,:)=ZERO
     
!FilterFeeders:all vars below control the exchange between benthos and pelagic.
! Mortality of FilterFeeders is initialized in PelB/PelGLobal.F90 because
! the mortality of benthic larvae is added to Yy3c
       jmYIc(:,:)=ZERO
       jBenFishInput(:)=ZERO 
       jPIY3c(:,:)=ZERO
       jP6Y3c(:,:)=ZERO
       jZEY3c(:,:)=ZERO
       jZIY3c(:)=ZERO
       jY3N4n(:)=ZERO
       jY3N1p(:)=ZERO
       jY3RIc(:)=ZERO
       jRIQIc(:)=ZERO
       jRIQIs(:)=ZERO
       jO2Y2o(:)=ZERO
       jCaCO3Y3c=ZERO;

       do i=1,iiBenLabileDetritus
          call  ResetSource_D2_vector(ppBenLabileDetritus(i,iiC))
          call  ResetSource_D2_vector(ppBenLabileDetritus(i,iiN))
          call  ResetSource_D2_vector(ppBenLabileDetritus(i,iiP))
       enddo
       do i=1,iiBenUrea
          call  ResetSource_D2_vector(ppBenUrea(i,iiN))
       enddo
#ifdef INCLUDE_BENCO2
       do i=1,iiBenthicCO2
          call  ResetSource_D2_vector(ppBenthicCO2(i,iiC))
       enddo
#endif
     case(2)
     ! Post processing: this done after call to Ben routines
        do i=1,iiBenPhyto
!JM         call flux_vector(iiBen,ppBenPhyto(i,iiS),ppQ6s,jPIQ6s(i,:))
         call flux_vector(iiBen,ppBenPhyto(i,iiS),ppQ6s,jPIQ6s(:,i))
        enddo
        call CalcLossHac()
    end select
  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
