!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: moddump
!
!  DESCRIPTION:
!   Take two relaxed stars from phanom_moddump_binary_system.f90
!   and place them in orbit at a given distance and with given 
!   initial conditions (irrotational - co-rotating)
!   Author: Jose Miguel Blanco (supervisor: Pablo Lorén-Aguilar)
!
!  REFERENCES: None
!
!  OWNER: Jose Miguel Blanco
!
!  $Id: -
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: centreofmass, extern_gwinspiral, externalforces, io,
!    options, part, physcon, prompting, timestep, units
!+
!--------------------------------------------------------------------------
module moddump
 implicit none

 logical, parameter :: use_defaults = .false.  ! if .true. will automatically use default values
                                               ! if .false., will ask user to prompt new values
 !--The default values
 real,    private   :: separation   = 50.
 logical, private   :: useirrinit        = .false.    !CHNGCODE

contains
!-----------------------------------------------------------------------

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu,nstar1,nstar2)
 use io,             only: iprint,fatal
 use prompting,      only: prompt
 use options,        only: iexternalforce,nfulldump,damp
 use part,           only: igas
 use units,          only: unit_velocity
 use physcon,        only: c,pi
 use timestep,       only: tmax, dtmax
 use centreofmass,   only: get_centreofmass,reset_centreofmass                 !CHNGCODE
 use externalforces, only: iext_gwinspiral
 use extern_gwinspiral, only: Nstar
                                                    !CHNGCODE Need to add  nstar1 and nstar2
 integer, intent(inout)    :: npart,nstar1,nstar2
 integer, intent(inout)    :: npartoftype(:)
 real,    intent(inout)    :: massoftype(:)
 real,    intent(inout)    :: xyzh(:,:),vxyzu(:,:)
 integer                   :: i
 real                      :: com(3),com_star1(3),com_star2(3),vcom(3)
 real                      :: rad1,rad2,mstar1,mstar2,mtotal, omega, omega2    !CHNGCODE added omega
 real                      :: xcm1,xcm2,ycm1,ycm2
 real                      :: c_code,tmax0 

 !
 !--Check particle numbers
 if (nstar1 <= 0 .or. nstar2 <= 0) call fatal('moddump','Require particle numbers in both stars')
 !
 !--Request parameters (unless hardcoded to use defaults)
 ! now determine the parameters of their new orbit
 if (.not.use_defaults) then
    call prompt('Enter desired separation:',separation,0.)
    call prompt('Use irrotational initial conditions?',useirrinit,.false.)              ! CHANGE to implement our desired initial conditions
 endif
 !
 !--Reset centre of mass location 
 call reset_centreofmass(nstar1,xyzh(:,1:nstar1),vxyzu(:,1:nstar1))
 call reset_centreofmass(nstar2,xyzh(:,nstar1+1:npart),vxyzu(:,nstar1+1:npart))

 mstar1 = nstar1 * massoftype(igas)
 mstar2 = nstar2 * massoftype(igas)
 mtotal = npart  * massoftype(igas)
 !
 !--Calcuate the new orbital parameters
 rad1   =  separation * mstar2/mtotal           ! distance of star 1 from the CoM
 rad2   =  separation * mstar1/mtotal           ! distance of star 2 from the CoM
 omega  =  sqrt(mtotal/separation**3)           ! rotational speed
 if (useirrinit) then
    omega2 = -1.0 * omega                       ! spin velocity of the stars, assuming sinchronaized spin speeds
 else
    omega2 = -0.0 * omega
 endif

 !
 !--Place stars on new orbits
 !  for simplicity, assume stars are on the x-axis                             !Take here into account both options co-rotating and irrotational
 vxyzu(1:3,:) = 0.0                             ! reset velocity
 xcm1 = rad1
 ycm1 = 0
 xcm2 = - rad2
 ycm2 = 0
 do i=1,nstar1
    xyzh(1,i)  =  xyzh(1,i) + rad1
    vxyzu(1,i) =  -xyzh(2,i)*(omega+omega2)+omega2*ycm1
    vxyzu(2,i) =  xyzh(1,i)*(omega+omega2)-omega2*xcm1
    vxyzu(3,i) =  0.0d0
 enddo
 do i=nstar1+1,npart
    xyzh(1,i)  = xyzh(1,i) - rad2
    vxyzu(1,i) =  -xyzh(2,i)*(omega+omega2)+omega2*ycm2
    vxyzu(2,i) =  xyzh(1,i)*(omega+omega2)-omega2*xcm2
    vxyzu(3,i) =  0.0d0 
 enddo
 !
 !--Set new runtime parameters
 tmax           = 1000.
 dtmax          =  100.
 damp           =    0.
 nfulldump      =   10
 iexternalforce =    0
 !
 !
 return
end subroutine modify_dump
!-----------------------------------------------------------------------
end module moddump
