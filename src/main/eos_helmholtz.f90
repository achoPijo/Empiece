!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: eos_mesa
!
!  DESCRIPTION: None
!
!  REFERENCES: None
!
!  OWNER: 
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: helmholtz_microphysiss
!+
!-----------------------------------------------------------------
module eos_helmholtz

 use helmholtz_microphysics

 implicit none

contains

!----------------------------------------------------------------
!+
!  subroutine initialises the helmholtz eos tables
!+
!----------------------------------------------------------------
subroutine init_eos_helmholtz(ierr)
 
 integer, intent(out) :: ierr
 ierr = 0

 call read_eos_helmholtz(ierr)

end subroutine init_eos_helmholtz

!----------------------------------------------------------------
!+
!  subroutine returns pressure and sound speed as a function of
!  temperature and density
!+
!----------------------------------------------------------------
subroutine get_eos_pressure_soundspeed_helmholtz(temp,den,abar,zbar,pres,sound)
 use io,         only:fatal
 real, intent(in)    :: temp, den, abar, zbar
 real, intent(out)   :: pres, sound
 integer             :: ierr

 ierr = 0
 
 call helmeos(temp,den,abar,zbar,pres,sound,ierr)

 if (ierr /= 0 ) call fatal(label,'Temperature or density values outside helmholtz free energy tables range')

end subroutine get_eos_pressure_soundspeed_helmholtz



!----------------------------------------------------------------
!+
!  ADD ANY ADDITIONAL SUBROUTINE
!+
!----------------------------------------------------------------


end module eos_helmholtz
