!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: part
!
!  DESCRIPTION:
!  This module contains main particle data for the code and
!  functions used to query particle properties not stored.
!
!  Basically this module defines any quantity that is
!  stored on the particles and defines how that storage
!  is implemented. Thus any routine which requires knowledge
!  of the specifics of this storage should be placed here.
!  (for example, any routine that copies all of the variables
!   stored on a given particle).
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: dim, io, mpiutils
!+
!--------------------------------------------------------------------------
module part
 use dim, only:maxp,ndivcurlv,ndivcurlB,maxvxyzu, &
          maxalpha,maxptmass,maxstrain, &
          mhd,maxmhd,maxBevol,maxvecp,maxp_h2,periodic, &
          maxgrav,ngradh,maxtypes,h2chemistry,gravity, &
          switches_done_in_derivs,maxp_dustfrac,use_dustfrac, &
          lightcurve,maxlum,nalpha,maxmhdni
 implicit none
 character(len=80), parameter, public :: &  ! module version
    modid="$Id$"

!
!--basic storage needed for read/write of particle data
!
 real :: xyzh(4,maxp)
 real :: vxyzu(maxvxyzu,maxp)
 real(kind=4) :: alphaind(nalpha,maxalpha)
 real(kind=4) :: divcurlv(ndivcurlv,maxp)
 real(kind=4) :: divcurlB(ndivcurlB,maxp)
 real(kind=4) :: Bevol(maxBevol,maxmhd)
 real(kind=4) :: Bxyz(3,maxvecp)
 character(len=*), parameter :: xyzh_label(4) = (/'x','y','z','h'/)
 character(len=*), parameter :: vxyzu_label(4) = (/'vx','vy','vz','u '/)
 character(len=*), parameter :: Bevol_label(4) = (/'Bx ','By ','Bz ','psi'/)
!
!--storage in divcurlv
!
 integer, parameter :: idivv = 1
 integer, parameter :: icurlvx = 2
 integer, parameter :: icurlvy = 3
 integer, parameter :: icurlvz = 4
 character(len=*), parameter :: divcurlv_label(4) = &
   (/'divv  ','curlvx','curlvy','curlvz'/)
!
!--storage in divcurlB
!
 integer, parameter :: idivB = 1
 integer, parameter :: icurlBx = 2
 integer, parameter :: icurlBy = 3
 integer, parameter :: icurlBz = 4
 character(len=*), parameter :: divcurlB_label(4) = &
   (/'divB  ','curlBx','curlBy','curlBz'/)
!
!--physical viscosity
!
 real(kind=4) :: straintensor(6,maxstrain)
!
!--H2 chemistry
!
 integer, parameter :: nabundances = 5
 integer, parameter :: ih2ratio  = 1 ! ratio of H2 to H
 integer, parameter :: iHI       = 2 ! HI abundance
 integer, parameter :: iproton   = 3 ! proton abundance
 integer, parameter :: ielectron = 4 ! electron abundance
 integer, parameter :: iCO       = 5 ! CO abundance
 real :: abundance(nabundances,maxp_h2)
 character(len=*), parameter :: abundance_label(5) = &
   (/'h2ratio','abHIq  ','abhpq  ','abeq   ','abco   '/)
!
!--one-fluid dust (small grains)
!
 real :: dustfrac(maxp_dustfrac)
 real :: dustevol(maxp_dustfrac)
 real :: deltav(3,maxp_dustfrac)
 character(len=*), parameter :: deltav_label(3) = &
   (/'deltavx','deltavy','deltavz'/)
!
!--sink particles
!
 integer, parameter :: nsinkproperties = 11
 integer, parameter :: ihacc  = 5 ! accretion radius
 integer, parameter :: ihsoft = 6 ! softening radius
 integer, parameter :: imacc  = 7 ! accreted mass
 integer, parameter :: ispinx = 8  ! spin angular momentum x
 integer, parameter :: ispiny = 9  ! spin angular momentum y
 integer, parameter :: ispinz = 10 ! spin angular momentum z
 integer, parameter :: i_tlast = 11 ! time of last injection
 real :: xyzmh_ptmass(nsinkproperties,maxptmass)
 real :: vxyz_ptmass(3,maxptmass)
 real :: fxyz_ptmass(4,maxptmass)
 integer :: nptmass = 0   ! zero by default
 real    :: epot_sinksink
 character(len=*), parameter :: xyzmh_ptmass_label(11) = &
  (/'x        ','y        ','z        ','m        ','h        ',&
    'hsoft    ','maccreted','spinx    ','spiny    ','spinz    ','tlast    '/)
 character(len=*), parameter :: vxyz_ptmass_label(3) = (/'vx','vy','vz'/)
!
!--self-gravity
!
 real(kind=4) :: poten(maxgrav)
!
!--Non-ideal MHD
!
 real         :: n_R(4,maxmhdni),n_electronT(maxmhdni)
 real(kind=4) :: ionfrac_eta(4,maxmhdni)
#ifdef NONIDEALMHD
 character(len=*), parameter :: ionfrac_eta_label(4) = (/'ne_by_n','eta_OR ','eta_HE ','eta_AD '/)
#endif
!
!--for analysis routines, do not allocate any more storage
!  than is strictly necessary
!
#ifdef ANALYSIS
 integer, parameter, private :: maxan = 0
 integer, parameter, private :: maxmhdan = 0
 integer, parameter, private :: maxdustan = 0
#else
 integer, parameter, private :: maxan = maxp
 integer, parameter, private :: maxmhdan = maxmhd
 integer, parameter, private :: maxdustan = maxp_dustfrac
#endif
!
!--lightcurves
!
 real(kind=4) :: luminosity(maxlum)
!
!--derivatives (only needed if derivs is called)
!
 real         :: fxyzu(maxvxyzu,maxan)
 real(kind=4) :: dBevol(maxBevol,maxmhdan)
 real(kind=4) :: divBsymm(maxmhdan)
 real         :: fext(3,maxan)
 real         :: ddustfrac(maxdustan)
!
!--storage associated with/dependent on timestepping
!
#ifdef IND_TIMESTEPS
 integer(kind=1)    :: ibin(maxan)
 integer(kind=1)    :: ibin_wake(maxan)
 integer(kind=1)    :: ibinold(maxan)
 integer(kind=1)    :: ibinsink(maxan)
 real(kind=4)       :: dt_in(maxan)
 real               :: twas(maxan)
#else
 integer(kind=1)    :: ibin_wake(1)
#endif
 integer, parameter :: maxphase = maxan
 integer, parameter :: maxgradh = maxan
 integer(kind=1)    :: iphase(maxphase)
 logical, public    :: all_active = .true.

 real(kind=4) :: gradh(ngradh,maxgradh)
!
!--storage associated with link list
!  (used for dead particle list also)
!
 integer :: ll(maxan)
!
!--size of the buffer required for transferring particle
!  information between MPI threads
!
 integer, parameter, private :: maxpd =  max(maxp,1) ! avoid divide by zero
 integer, parameter :: ipartbufsize = 4 &  ! xyzh
   +2*maxvxyzu                          &  ! vxyzu, fxyzu_prev
   +maxalpha/maxpd                      &  ! alphaind
   +ngradh*maxgradh/maxpd               &  ! gradh
   +maxphase/maxpd                      &  ! iphase
#ifdef IND_TIMESTEPS
   +2 +maxvxyzu                         &  ! ibin, divv, fxyzu
   +(maxmhd/maxpd)*maxBevol +3*(maxvecp/maxpd)  &  ! dB/dt, Bxyz
#endif
   +(maxmhd/maxpd)*                      &  ! (mhd quantities)
    (2*maxBevol                         &  ! Bevol, dBevol_prev
   +3*(maxvecp/maxpd))                     ! Bxyz

 real            :: hfact,Bextx,Bexty,Bextz
 integer         :: npart
 integer(kind=8) :: ntot
 integer         :: ideadhead = 0

 integer :: npartoftype(maxtypes)
 real    :: massoftype(maxtypes)
!
!--labels for each type
!  NOTE: If new particle is added, and it is allowed to be accreted onto
!        a sink particle, add it to the list in 'is_accretable' below
!  NOTE: set_boundaries_to_active = .true., but will be set to .false. at the
!        end of initial.  This will allow boundary particles to always be
!         initialised (even on restarts where not all arrays, e.g. gradh,
!         are not saved)
!
 integer, parameter :: igas  = 1
 integer, parameter :: idust = 2
 integer, parameter :: iboundary = 3
 integer, parameter :: istar = 4
 integer, parameter :: idarkmatter = 5
 integer, parameter :: ibulge = 6
 integer, parameter :: iunknown = 0
 logical            :: set_boundaries_to_active = .true.
 character(len=5), dimension(maxtypes), parameter :: &
    labeltype = (/'gas  ','dust ','bound','star ','darkm','bulge'/)
!
!--generic interfaces for routines
!
 interface hrho
    module procedure hrho4,hrho8,hrho4_pmass,hrho8_pmass,hrhomixed_pmass
 end interface hrho

 private :: hrho4,hrho8,hrho4_pmass,hrho8_pmass,hrhomixed_pmass

contains
!----------------------------------------------------------------
!+
!  this function determines the mass of the particle
!  where use_gas == .not. (maxphase==maxp)
!  currently used only in readwrite_dump
!+
!----------------------------------------------------------------
real function get_pmass(i,use_gas)
  integer, intent(in) :: i
  logical, intent(in) :: use_gas

  if (use_gas) then
    get_pmass = massoftype(igas)
  else
    if (iphase(i) /= 0) then
       get_pmass = massoftype(iamtype(iphase(i)))
    else
       get_pmass = massoftype(igas)
    endif
  endif

end function get_pmass
!
!----------------------------------------------------------------
!+
!  this function gives rho as a function of h
!  (used by output routines to get rho given that we only
!   store h)
!+
!----------------------------------------------------------------
pure real function rhoh(hi,pmassi)
 real, intent(in) :: hi,pmassi

 rhoh = pmassi*(hfact/abs(hi))**3

end function rhoh

!----------------------------------------------------------------
!+
!  this function gives dh/drho as a function of h
!+
!----------------------------------------------------------------
real function dhdrho(hi,pmassi)
 real, intent(in) :: hi,pmassi
 real :: rhoi

 rhoi   = rhoh(hi,pmassi)
 dhdrho = -hi/(3.*rhoi)

end function dhdrho
!----------------------------------------------------------------
!+
!  this subroutine does both of the above
!  (routine is an optimisation - also returns divisions)
!+
!----------------------------------------------------------------
subroutine rhoanddhdrho(hi,hi1,rhoi,rho1i,dhdrhoi,pmassi)
 real,         intent(in)  :: hi, pmassi
 real(kind=8), intent(out) :: hi1
 real,         intent(out) :: rhoi
 real,         intent(out) :: rho1i,dhdrhoi
 real, parameter :: third = 1./3.

 hi1 = 1./abs(hi)
 rhoi = pmassi*(hfact*hi1)**3
 rho1i = 1./rhoi
 dhdrhoi = -third*hi*rho1i

end subroutine rhoanddhdrho

!----------------------------------------------------------------
!+
!  this function gives h as a function of rho
!+
!----------------------------------------------------------------
real(kind=4) function hrho4(rhoi)
 real(kind=4), intent(in) :: rhoi

 hrho4 = real(hfact*(massoftype(1)/abs(rhoi))**(1./3.),kind=4)

end function hrho4

real(kind=8) function hrho8(rhoi)
 real(kind=8), intent(in) :: rhoi

 hrho8 = hfact*(massoftype(1)/abs(rhoi))**(1.d0/3.d0)

end function hrho8

real(kind=4) function hrho4_pmass(rhoi,pmassi)
 real(kind=4), intent(in) :: rhoi,pmassi

 hrho4_pmass = real(hfact*(pmassi/abs(rhoi))**(1./3.),kind=4)

end function hrho4_pmass

real(kind=8) function hrho8_pmass(rhoi,pmassi)
 real(kind=8), intent(in) :: rhoi,pmassi

 hrho8_pmass = hfact*(pmassi/abs(rhoi))**(1.d0/3.d0)

end function hrho8_pmass

real(kind=8) function hrhomixed_pmass(rhoi,pmassi)
 real(kind=8), intent(in) :: rhoi
 real(kind=4), intent(in) :: pmassi

 hrhomixed_pmass = hfact*(pmassi/abs(rhoi))**(1.d0/3.d0)

end function hrhomixed_pmass

!----------------------------------------------------------------
!+
!  query function returning whether or not a particle is dead
!  (currently indicated by having a negative smoothing length)
!+
!----------------------------------------------------------------
logical function isdead(i)
 integer, intent(in) :: i

 ! h = 0 indicates a dead particle
 if (abs(xyzh(4,i)) < tiny(xyzh)) then
    isdead = .true.
 else
    isdead = .false.
 endif

end function isdead

pure logical function isdeadh(hi)
 real, intent(in) :: hi

 ! h = 0 indicates a dead particle
 if (abs(hi) < tiny(hi)) then
    isdeadh = .true.
 else
    isdeadh = .false.
 endif

end function isdeadh

pure logical function isdead_or_accreted(hi)
 real, intent(in) :: hi

 ! h <= 0 indicates either dead or accreted
 if (hi < tiny(hi)) then
    isdead_or_accreted = .true.
 else
    isdead_or_accreted = .false.
 endif

end function isdead_or_accreted

!----------------------------------------------------------------
!+
! routine which kills a particle and adds it to the dead list
!+
!----------------------------------------------------------------
subroutine kill_particle(i)
 integer, intent(in) :: i

 xyzh(4,i) = 0.
!$omp critical
 ll(i) = ideadhead
 ideadhead = i
!$omp end critical

end subroutine kill_particle

!----------------------------------------------
!+
!  functions to deconstruct iphase
!  abs(iphase) is the particle type
!  sign(iphase) gives whether it is active/inactive
!+
!----------------------------------------------
pure integer(kind=1) function isetphase(itype,iactive)
 integer, intent(in) :: itype
 logical, intent(in) :: iactive

 if ((set_boundaries_to_active .and. itype==iboundary)   .or. &
     (iactive                  .and. itype/=iboundary) ) then
    isetphase = int(itype,kind=1)
 else
    isetphase = -int(abs(itype),kind=1)
 endif

end function isetphase

pure subroutine get_partinfo(iphasei,isactive,isdust,itype)
 integer(kind=1), intent(in)  :: iphasei
 logical,         intent(out) :: isactive,isdust
 integer,         intent(out) :: itype

! isactive = iactive(iphasei)
! itype = iamtype(iphasei)
! isdust = itype==idust

!--inline versions of above (for speed)
 if (iphasei > 0) then
    isactive = .true.
    itype    = iphasei
 else
    isactive = .false.
    itype    = -iphasei
 endif
#ifdef DUST
 isdust = itype==idust
#else
 isdust = .false.
#endif

 return
end subroutine get_partinfo

pure logical function iactive(iphasei)
 integer(kind=1), intent(in) :: iphasei

 iactive = (iphasei > 0)

end function iactive

pure elemental integer function iamtype(iphasei)
 integer(kind=1), intent(in) :: iphasei

 iamtype = abs(iphasei)

end function iamtype

pure elemental function iamtype_int1(iphasei)
 integer(kind=1), intent(in) :: iphasei
 integer(kind=1) :: iamtype_int1

 iamtype_int1 = abs(iphasei)

end function iamtype_int1

pure function iamtype_int11(iphasei)
 integer(kind=1), intent(in) :: iphasei
 integer(kind=1) :: iamtype_int11

 iamtype_int11 = abs(iphasei)

end function iamtype_int11

pure elemental logical function iamgas(iphasei)
 integer(kind=1), intent(in) :: iphasei
 integer :: itype

 itype = iamtype(iphasei)
 iamgas = int(itype)==igas

end function iamgas

pure elemental logical function iamdust(iphasei)
 integer(kind=1), intent(in) :: iphasei
 integer :: itype

 itype = iamtype(iphasei)
 iamdust = int(itype)==idust

end function iamdust

pure integer function get_ntypes(noftype)
 integer, intent(in) :: noftype(:)
 integer :: i

 get_ntypes = 0
 do i=1,size(noftype)
    if (noftype(i) > 0) get_ntypes = i
 enddo

end function get_ntypes

!-----------------------------------------------------------------------
!+
!  Determine if particle is of a type that is accretable
!    Modify the if-statement to include all the types of particles
!    that can be accreted onto the sink
!+
!-----------------------------------------------------------------------
pure logical function is_accretable(itype)
 integer, intent(in)  :: itype

 if (itype==igas .or. itype==idust) then
   is_accretable = .true.
 else
   is_accretable = .false.
 endif

end function is_accretable

!----------------------------------------------------------------
!+
!  utility function for setup routines to set initial value
!  of iphase (assumes particle is active)
!+
!----------------------------------------------------------------
subroutine set_particle_type(i,itype)
 use io, only:fatal
 integer, intent(in) :: i,itype

 if (maxphase==maxp) then
    iphase(i) = isetphase(itype,iactive=.true.)
 elseif (itype /= igas) then
    call fatal('set_particle_type','attempt to setup a particle of type > 1, but iphase not allocated')
 endif

end subroutine set_particle_type

!----------------------------------------------------------------
!+
! routine which copies a particle from one location to another
! (prior to a derivs evaluation - so no derivs required)
!+
!----------------------------------------------------------------
subroutine copy_particle(src, dst)
 integer, intent(in) :: src, dst

 xyzh(:,dst)  = xyzh(:,src)
 vxyzu(:,dst) = vxyzu(:,src)
 if (mhd) then
    Bevol(:,dst) = Bevol(:,src)
 endif
 if (ndivcurlv  > 0) divcurlv(:,dst)  = divcurlv(:,src)
 if (maxalpha ==maxp) alphaind(:,dst) = alphaind(:,src)
 if (maxgradh ==maxp) gradh(:,dst)    = gradh(:,src)
 if (maxphase ==maxp) iphase(dst)   = iphase(src)
 if (maxgrav  ==maxp) poten(dst) = poten(src)
#ifdef IND_TIMESTEPS
 ibin(dst) = ibin(src)
#endif

 return
end subroutine copy_particle

!----------------------------------------------------------------
!+
! routine which copies a particle from one location to another
! (copies everything which is stored on a particle)
!
! Note that link list information CANNOT be copied so link list
! must be rebuilt after a copy operation.
!+
!----------------------------------------------------------------
subroutine copy_particle_all(src,dst)
 integer, intent(in) :: src,dst

 xyzh(:,dst)  = xyzh(:,src)
 vxyzu(:,dst) = vxyzu(:,src)
 fxyzu(:,dst) = fxyzu(:,src)
 fext(:,dst) = fext(:,src)
 if (mhd) then
    Bevol(:,dst)  = Bevol(:,src)
    dBevol(:,dst) = dBevol(:,src)
    if (maxvecp==maxp) Bxyz(:,dst)   = Bxyz(:,src)
    divBsymm(dst) = divBsymm(src)
 endif
 if (ndivcurlv > 0) divcurlv(:,dst) = divcurlv(:,src)
 if (ndivcurlB > 0) divcurlB(:,dst) = divcurlB(:,src)
 if (maxalpha ==maxp) alphaind(:,dst) = alphaind(:,src)
 if (maxgradh ==maxp) gradh(:,dst) = gradh(:,src)
 if (maxphase ==maxp) iphase(dst) = iphase(src)
#ifdef IND_TIMESTEPS
 ibin(dst) = ibin(src)
#endif

 return
end subroutine copy_particle_all

!------------------------------------------------------------------
!+
! routine which reorders the particles according to an input list
! (prior to a derivs evaluation - so no derivs required)
! (allocates temporary arrays for each variable, so use with caution)
!+
!------------------------------------------------------------------
subroutine reorder_particles(iorder,np)
 integer, intent(in) :: iorder(:)
 integer, intent(in) :: np

 call copy_array(xyzh(:,1:np), iorder(1:np))
 call copy_array(vxyzu(:,1:np),iorder(1:np))
 call copy_array(fext(:,1:np),iorder(1:np))
 if (mhd) then
    call copy_arrayr4(Bevol(:,1:npart),iorder(1:np))
    !--also copy the Bfield here, as this routine is used in setup routines
    if (maxvecp        ==maxp) call copy_arrayr4(Bxyz(:,1:np),        iorder(1:np))
 endif
 if (ndivcurlv > 0)     call copy_arrayr4(divcurlv(:,1:np),iorder(1:np))
 if (maxalpha ==maxp) call copy_arrayr4(alphaind(:,1:np),  iorder(1:np))
 if (maxgradh ==maxp) call copy_arrayr4(gradh(:,1:np),  iorder(1:np))
 if (maxphase ==maxp) call copy_arrayint1(iphase(1:np),iorder(1:np))
 if (maxgrav  ==maxp) call copy_array1(poten(1:np),  iorder(1:np))
#ifdef IND_TIMESTEPS
 call copy_arrayint1(ibin(1:np),  iorder(1:np))
#endif

 return
end subroutine reorder_particles

!-----------------------------------------------------------------------
!+
!  routine to compactify the list of particles by removing dead
!  particles from the list
!  (could be openMP parallel if we sent in ndead)
!+
!-----------------------------------------------------------------------
subroutine shuffle_part(np)
 use io, only:fatal
 integer, intent(inout) :: np
 integer :: newpart

 do while (ideadhead /= 0)
    newpart = ideadhead
    if (newpart < np) then
       if (.not.isdead(np)) then
          !if (.not.isdead(newpart)) call fatal('shuffle','corrupted dead list',newpart)
          call copy_particle_all(np,newpart)
          ideadhead = ll(newpart)
       endif
       np = np - 1
    else
       ideadhead = ll(newpart)
    endif
    if (np <= 0) call fatal('shuffle','npart < 0')
 enddo

 return
end subroutine shuffle_part

!-----------------------------------------------------------------------
!+
!  routine to remove dead or accreted particles
!  uses the routines above for efficiency
!+
!-----------------------------------------------------------------------
subroutine delete_dead_or_accreted_particles(npart,npartoftype)
 integer, intent(inout) :: npart,npartoftype(:)
 integer :: i,itype

 do i=1,npart
    if (isdead_or_accreted(xyzh(4,i))) then
       ! get the type so we know how to decrement npartoftype
       if (maxphase==maxp) then
          itype = iamtype(iphase(i))
       else
          itype = igas
       endif
       npartoftype(itype) = npartoftype(itype) - 1
       call kill_particle(i)
    endif
 enddo
 call shuffle_part(npart)

 return
end subroutine delete_dead_or_accreted_particles

!----------------------------------------------------------------
!+
!   change the position and status of a dead particle
!
!+
!----------------------------------------------------------------

subroutine change_status_pos(npart,x,y,z,h,vx,vy,vz)

 integer, intent(in) :: npart
 real, intent (in) :: x,y,z,h
 real, intent (in) :: vx,vy,vz
 integer  :: i,ix

 ix=0

 do i=1,npart
    if (isdead_or_accreted(xyzh(4,i))) then
       ix=i
       exit
    endif
 enddo

 xyzh(1,ix)=x
 xyzh(2,ix)=y
 xyzh(3,ix)=z
 xyzh(4,ix)=h
 vxyzu(1,ix)=vx
 vxyzu(2,ix)=vy
 vxyzu(3,ix)=vz

 return

end subroutine change_status_pos

!----------------------------------------------------------------
!+
!  pack particle information into a contiguous buffer
!  to send to another processor
!+
!----------------------------------------------------------------
subroutine fill_sendbuf(i,xtemp)
 use io,       only:fatal
 use mpiutils, only:fill_buffer
 integer, intent(in)  :: i
 real,    intent(out) :: xtemp(ipartbufsize)
 integer :: nbuf
!
!--package particle information into one simple wrapper
!
 nbuf = 0
!--NB: could use MPI_PACK here...
 if (i > 0) then
    call fill_buffer(xtemp,xyzh(:,i),nbuf)
    call fill_buffer(xtemp,vxyzu(:,i),nbuf)
    !call fill_buffer(xtemp,fxyzu_prev(:,i),nbuf)
    if (maxalpha==maxp) then
       call fill_buffer(xtemp,alphaind(:,i),nbuf)
    endif
    if (maxgradh==maxp) then
       call fill_buffer(xtemp,gradh(:,i),nbuf)
    endif
    if (mhd) then
       call fill_buffer(xtemp,Bevol(:,i),nbuf)
       !call fill_buffer(xtemp,dBevol_prev(:,i),nbuf)
       if (maxvecp==maxp) then
          call fill_buffer(xtemp,Bxyz(:,i),nbuf)
       endif
    endif
    if (maxphase==maxp) then
       call fill_buffer(xtemp,iphase(i),nbuf)
    endif
#ifdef IND_TIMESTEPS
    call fill_buffer(xtemp,ibin(i),nbuf)
    !--inactive particles require derivs sent
    call fill_buffer(xtemp,fxyzu(:,i),nbuf)
    call fill_buffer(xtemp,fext(:,i),nbuf)
    if (ndivcurlv >= 1) then
       call fill_buffer(xtemp,divcurlv(1,i),nbuf)
    endif
    if (mhd) then
       call fill_buffer(xtemp,dBevol(:,i),nbuf)
       if (maxvecp==maxp) then
          call fill_buffer(xtemp,Bxyz(:,i),nbuf)
       endif
    endif
#endif
 endif
 if (nbuf /= ipartbufsize) call fatal('fill_sendbuf','error in send buffer size')

 return
end subroutine fill_sendbuf

!----------------------------------------------------------------
!+
!  unpack particle information from the send buffer
!  after receiving from another processor
!+
!----------------------------------------------------------------
subroutine unfill_buffer(ipart,xbuffer)
 integer, intent(in) :: ipart
 real,    intent(in) :: xbuffer(ipartbufsize)
 integer :: i1,i2

 i1 = 1
 i2 = 4
 xyzh(:,ipart) = xbuffer(i1:i2)
 i1 = i2+1
 i2 = i2+maxvxyzu
 vxyzu(:,ipart) = xbuffer(i1:i2)
 !i1 = i2 + 1
 !i2 = i2 + maxvxyzu
 !fxyzu_prev(:,ipart) = xbuffer(i1:i2)
 if (maxalpha==maxp) then
    i1 = i2 + 1
    i2 = i2 + nalpha
    alphaind(:,ipart) = real(xbuffer(i1:i2),kind=kind(alphaind))
 endif
 if (maxgradh==maxp) then
    i1 = i2 + 1
    i2 = i2 + ngradh
    gradh(:,ipart) = real(xbuffer(i1:i2),kind=kind(gradh))
 endif
 if (mhd) then
    i1 = i2 + 1
    i2 = i2 + maxBevol
    Bevol(:,ipart) = real(xbuffer(i1:i2),kind=kind(Bevol))
    !i1 = i2 + 1
    !i2 = i2 + maxBevol
    !dBevol_prev(:,ipart) = real(xbuffer(i1:i2),kind=kind(dBevol_prev))
    if (maxvecp==maxp) then
       i1 = i2 + 1
       i2 = i2 + 3
       Bxyz(:,ipart) = real(xbuffer(i1:i2),kind=kind(Bxyz))
    endif
 endif
 if (maxphase==maxp) then
    i1 = i2 + 1
    i2 = i2 + 1
    iphase(ipart) = nint(xbuffer(i1),kind=1)
 endif
#ifdef IND_TIMESTEPS
 i1 = i2 + 1
 i2 = i2 + 1
 ibin(ipart) = nint(xbuffer(i1),kind=1)

 !--receive derivs (strictly only necessary for inactive parts)
 i1 = i2 + 1
 i2 = i2 + maxvxyzu
 fxyzu(:,ipart) = xbuffer(i1:i2)
 i1 = i2 + 1
 i2 = i2 + 3
 fext(:,ipart) = xbuffer(i1:i2)
 if (ndivcurlv > 0) then
    i1 = i2 + 1
    i2 = i2 + 1
    divcurlv(1,ipart) = real(xbuffer(i1),kind=kind(divcurlv))
 endif
 if (mhd) then
    i1 = i2 + 1
    i2 = i2 + maxBevol
    dBevol(:,ipart) = real(xbuffer(i1:i2),kind=kind(dBevol))
    if (maxvecp==maxp) then
       i1 = i2 + 1
       i2 = i2 + 3
       Bxyz(:,ipart) = real(xbuffer(i1:i2),kind=kind(Bxyz))
    endif
 endif
#endif
!--just to be on the safe side, set other things to zero
 if (mhd) then
    divBsymm(ipart) = 0.
 endif

 return
end subroutine unfill_buffer

!----------------------------------------------------------------
!+
!  utility to reorder an array
!  (rank 2 arrays)
!+
!----------------------------------------------------------------

subroutine copy_array(array,ilist)
 real,    intent(inout) :: array(:,:)
 integer, intent(in)    :: ilist(:)
 real :: arraytemp(size(array(1,:)))
 integer :: i

 do i=1,size(array(:,1))
    arraytemp(:) = array(i,ilist(:))
    array(i,:) = arraytemp
 enddo

 return
end subroutine copy_array

!----------------------------------------------------------------
!+
!  utility to reorder an array
!  (real4, rank 2 arrays)
!+
!----------------------------------------------------------------

subroutine copy_arrayr4(array,ilist)
 real(kind=4), intent(inout) :: array(:,:)
 integer,      intent(in)    :: ilist(:)
 real(kind=4) :: arraytemp(size(array(1,:)))
 integer :: i

 do i=1,size(array(:,1))
    arraytemp(:) = array(i,ilist(:))
    array(i,:) = arraytemp
 enddo

 return
end subroutine copy_arrayr4

!----------------------------------------------------------------
!+
!  utility to reorder an array
!  (real4, rank 1 arrays)
!+
!----------------------------------------------------------------

subroutine copy_array1(array,ilist)
 real(kind=4), intent(inout) :: array(:)
 integer,      intent(in)    :: ilist(:)
 real(kind=4) :: arraytemp(size(array(:)))

 arraytemp(:) = array(ilist(:))
 array = arraytemp

 return
end subroutine copy_array1

!----------------------------------------------------------------
!+
!  utility to reorder an array
!  (real4, rank 1 arrays)
!+
!----------------------------------------------------------------

subroutine copy_arrayint1(iarray,ilist)
 integer(kind=1), intent(inout) :: iarray(:)
 integer,         intent(in)    :: ilist(:)
 integer(kind=1) :: iarraytemp(size(iarray(:)))

 iarraytemp(:) = iarray(ilist(:))
 iarray = iarraytemp

 return
end subroutine copy_arrayint1

!----------------------------------------------------------------
!+
!  Delete particles outside of a defined box
!+
!----------------------------------------------------------------
subroutine delete_particles_outside_box(xmin, xmax, ymin, ymax, zmin, zmax)
  real, intent(in) :: xmin, xmax, ymin, ymax, zmin, zmax

  integer :: i
  real :: x, y, z, h

  do i=1,npart
    x = xyzh(1,i)
    y = xyzh(2,i)
    z = xyzh(3,i)
    h = xyzh(4,i)
    if (x  <  xmin .or. x  >  xmax .or. y  <  ymin .or. y  >  ymax .or. z  <  zmin .or. z  >  zmax) then
      xyzh(4,i) = -abs(h)
    endif
  enddo
end subroutine

!----------------------------------------------------------------
!+
!  Delete particles outside of a defined sphere
!+
!----------------------------------------------------------------
subroutine delete_particles_outside_sphere(center, radius)
  real, intent(in) :: center(3), radius

  integer :: i
  real :: r(3), radius_squared

  radius_squared = radius**2
  do i=1,npart
    r = xyzh(1:3,i) - center
    if (dot_product(r,r)  >  radius_squared) then
     xyzh(4,i) = -abs(xyzh(4,i))
    endif
  enddo
end subroutine

!----------------------------------------------------------------
!+
!  Delete particles outside of a defined cylinder
!+
!----------------------------------------------------------------
subroutine delete_particles_outside_cylinder(center, radius, zmax)
  real, intent(in) :: center(3), radius, zmax

  integer :: i
  real :: x, y, z, rcil

  do i=1,npart
    x = xyzh(1,i)
    y = xyzh(2,i)
    z = xyzh(3,i)
    rcil=sqrt((x-center(1))**2+(y-center(2))**2)

    if (rcil>radius .or. abs(z)>zmax) then
       call kill_particle(i)
    endif
  enddo
end subroutine

end module part
