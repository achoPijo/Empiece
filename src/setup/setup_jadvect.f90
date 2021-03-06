!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: setup
!
!  DESCRIPTION:
!   Setup for the MHD current loop advection problem
!
!  REFERENCES:
!    Gardiner and Stone (2005)
!    Rosswog and Price (2007)
!    Stone et al. (2008)
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: boundary, dim, io, mpiutils, part, physcon, prompting,
!    setup_params, unifdis
!+
!--------------------------------------------------------------------------
module setup
 implicit none
 public :: setpart

 private

contains

!----------------------------------------------------------------
!+
!  setup for MHD current loop advection problem
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use dim,          only:maxp,maxvxyzu
 use setup_params, only:rhozero,ihavesetupB
 use unifdis,      only:set_unifdis
 use boundary,     only:set_boundary,xmin,ymin,zmin,xmax,ymax,zmax,dxbound,dybound,dzbound
 use part,         only:Bevol,mhd
 use io,           only:master,real4
 use prompting,    only:prompt
 use mpiutils,     only:bcast_mpi
 use physcon,      only:pi
 integer,           intent(in)    :: id
 integer,           intent(out)   :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma
 real,              intent(in)    :: hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real :: deltax,totmass,dz,rcyl
 integer :: i,nx
 real :: vzero,przero,uuzero
 real :: Azero,rloop,gam1,costheta,sintheta
!
!--general parameters
!
 time = 0.
 gamma = 5./3.
!
!--setup parameters
!
 vzero = sqrt(5.)
 przero = 1.0
 rhozero = 1.0
 Azero = 1.e-3
 rloop = 0.3

 gam1 = gamma - 1.
 uuzero = przero/(gam1*rhozero)

 print "(/,a)",' Three dimensional current loop advection problem '
 print 10, Azero,rloop,rhozero,przero
10 format(/,' Azero   = ',f6.3,', radius of loop = ',f6.3,/,&
            ' density = ',f6.3,',       pressure = ',f6.3,/)

 if (maxvxyzu < 4) then
    polyk = przero/rhozero**gamma
 else
    polyk = 0.
 endif

 nx = 128
 if (id==master) call prompt('Enter resolution (number of particles in x)',nx,8)
 call bcast_mpi(nx)
!
!--boundaries
!
 dz = 4.*sqrt(6.)/nx
 call set_boundary(-1.,1.,-0.5,0.5,-dz,dz)
 deltax = dxbound/nx

 costheta = dxbound/sqrt(dxbound**2 + dybound**2)
 sintheta = dybound/sqrt(dxbound**2 + dybound**2)

 call set_unifdis('closepacked',id,master,xmin,xmax,ymin,ymax,zmin,zmax,deltax,hfact,npart,xyzh)
 npartoftype(:) = 0
 npartoftype(1) = npart

 totmass = rhozero*dxbound*dybound*dzbound
 massoftype = totmass/npart
 print*,'npart = ',npart,' particle mass = ',massoftype(1)

 do i=1,npart
    vxyzu(1,i) = vzero*costheta
    vxyzu(2,i) = vzero*sintheta
    vxyzu(3,i) = 0.1*vzero
    if (maxvxyzu >= 4) vxyzu(4,i) = uuzero
    if (mhd) then
       Bevol(:,i) = 0.
       rcyl = sqrt(xyzh(1,i)**2 + xyzh(2,i)**2)
       if (rcyl < rloop) then
          Bevol(1,i) = real4(Azero*(-xyzh(2,i)/rcyl))
          Bevol(2,i) = real4(Azero*(xyzh(1,i)/rcyl))
       endif
    endif
 enddo

 if (mhd) ihavesetupB = .true.

end subroutine setpart

end module setup

