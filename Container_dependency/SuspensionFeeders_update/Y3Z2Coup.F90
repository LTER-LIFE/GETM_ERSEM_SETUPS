#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: BentoPelCoup
!
! DESCRIPTION
!   !

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine Y3Z2CoupDynamics
!
! !USES:


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN,DONE   ,LOGUNIT
  use mem, ONLY:    Yy3c,Z2c,sediMeZ,  D2STATE
  use mem, ONLY: max_change_per_step,NO_BOXES_XY, BoxNumberXY, BoxNumber, &
    ppYy3c,ppYy3n,ppYy3p,ppZ2c,ppZ2n,ppZ2p,iiPel,iiBen,iiZ2,iiYy3,&
    PelBoxAbove, Depth,jspaY3c,TauBed, flux

  use mem_Param,  ONLY: CalcMesoZooPlankton,CalcSuspensionFeeders
  use mem_FilterFeeder,ONLY:p_qnY3c=>p_qnc,p_qpY3c=>p_qpc,p_xtauC,p_sxresus, &
                                           p_xsteep
  use mem_MesoZoo,ONLY:p_qnMec=>p_qnc,p_qpMec=>p_qpc
  use mem_FilterFeeder,ONLY:p_heighty,p_pYy3Z2
  use LimitRates, ONLY:LimitChange
  use constants,only:ANY,POSITIVE
  use botflux,ONLY:addbotflux
  use global_interface,only:CorrectConcNearBed,CalcPelMassInM2

! use mem,  ONLY: Output2d_1,Output2d_2,Output2d_3, Output2d_4,Yy3c

!  
!
! !AUTHORS
!   Piet Ruardij  
!
!
! !REVISION_HISTORY
!   Created at Mon Apr 19 00:08:12 CEST 2004
!
!
! COPYING
!   
!   Copyright (C) 2006 P. Ruardij & M. Vichi
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
  real(RLEN),dimension(NO_BOXES_XY)   ::cx_any
  real(RLEN)  :: corr,px_any,rx_any
  real(RLEN)  :: ruYy3c,ruYy3n,ruYy3p
  real(RLEN)  :: toZ2c
  real(RLEN)  :: toZ2c_flux  !JM added

!write(LOGUNIT,*) 'subroutine Y3Z2coup'
!stop

  if (CalcMesoZooPlankton(iiZ2) .and. CalcSuspensionFeeders(iiYy3) ) then

    cx_any=CalcPelMassInM2(ppZ2c)
!write(LOGUNIT,*) 'before loop'
!stop
    DO BoxNumberXY=1,NO_BOXES_XY
!write(LOGUNIT,*) 'in loop loop',BoxNumberXY,NO_BOXES_XY
!stop
      BoxNumber=PelBoxAbove(BoxNumberXY)
!write(LOGUNIT,*) 'BoxNumber',BoxNumber
!stop

      ! flux Yy3 -> Z2
      toZ2c=jspaY3c(BoxNumberXY)* p_pYy3Z2 

!write(LOGUNIT,*) 'before resuspension'
!stop
      !resuspension ude of hight current.
      ! the positive flux ( flux to water column) is limtied
      ! by the bioasmass of Yy3. test on positive changes in this case.
      ! later toZ2c is subracted form Yy3c
      px_any=(DONE-TauBed(BoxNumberXY)/p_xtauC(iiYy3))
      rx_any=exp(-p_xsteep *px_any)*p_sxresus(iiYy3)*Yy3c(BoxNumberXY)
      call LimitChange(POSITIVE,rx_any,Yy3c(BoxNumberXY),max_change_per_step)
      toZ2c=toZ2c+rx_any
!write(LOGUNIT,*) 'before correctconc'
!stop

      call CorrectConcNearBed(Depth(BoxNumber),sediMeZ(iiZ2,BoxNumber) &
                                   , p_heighty,1.0D+40,corr)
      ! flux Z2 -> Yy3 
      ruYy3c  = min(Depth(BoxNumberXY),p_heighty)*corr* Z2c(BoxNumber)
      call LimitChange(POSITIVE,ruYy3c,Z2c(BoxNumber)*Depth(BoxNumber),&
                                                    max_change_per_step)
      call LimitChange(POSITIVE,ruYy3c,cx_any(BoxNumberXY), max_change_per_step)
!write(LOGUNIT,*) 'before use p_qnMEC'
!stop
      ruYy3n  =   ruYy3c*p_qnMec(iiZ2)
      ruYy3p  =   ruYy3c*p_qpMec(iiZ2)

      call addbotflux(ANY,BoxNumberXY,iiBen,ppYy3c,iiPel,ppZ2c,toZ2c-ruYy3c)
!JM      if (ppZ2n>0) then
!JM        call addbotflux(ANY,BoxNumberXY,iiBen,ppYy3n,iiPel, &
!JM                                             ppZ2n,toZ2c*p_qnY3c-ruYy3n)
!JM        call addbotflux(ANY,BoxNumberXY,iiBen,ppYy3p,iiPel, &
!JM                                             ppZ2p,toZ2c*p_qpY3c-ruYy3p)
!JM      else
!JM        call flux(BoxNumberXY,iiBen,ppYy3n,ppYy3n,-toZ2c*p_qnY3c+ruYy3n)
!JM        call flux(BoxNumberXY,iiBen,ppYy3p,ppYy3p,-toZ2c*p_qpY3c+ruYy3p)
!JM      endif
      if (ppZ2n>0) then
        toZ2c_flux=toZ2c*p_qnY3c-ruYy3n
        call addbotflux(ANY,BoxNumberXY,iiBen,ppYy3n,iiPel, &
                                             ppZ2n,toZ2c_flux)
        toZ2c_flux=toZ2c*p_qpY3c-ruYy3p
        call addbotflux(ANY,BoxNumberXY,iiBen,ppYy3p,iiPel, &
                                             ppZ2p,toZ2c_flux)
      else
        toZ2c_flux=-toZ2c*p_qnY3c+ruYy3n
        call flux(BoxNumberXY,iiBen,ppYy3n,ppYy3n,toZ2c_flux)
        toZ2c_flux=-toZ2c*p_qpY3c+ruYy3p
        call flux(BoxNumberXY,iiBen,ppYy3p,ppYy3p,toZ2c_flux)
      endif
     end DO
    
   endif
!write(LOGUNIT,*) 'end subroutine Y3Z2coup'
!stop

  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
