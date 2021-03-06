!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  PROGRAM: phantom
!
!  DESCRIPTION: The Phantom SPH code, by Daniel Price.
!
!  This code is designed to be an ultra-sleek, ultra-low-memory,
!  code for high resolution SPH simulations
!
!  The requirements mean we need to store as few quantities as possible
!  (aim is to be able to run 10^7 particles in under 1Gb)
!  and to use the fastest possible implementation
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  USAGE: phantom infilename
!
!  DEPENDENCIES: dim, evolve, initial, io, mpiutils, test
!+
!--------------------------------------------------------------------------
program phantom
 use dim,             only:tagline
 use mpiutils,        only:init_mpi, finalise_mpi
 use initial,         only:initialise,startrun,endrun
 use io,              only:id,master,nprocs,set_io_unit_numbers,die
 use evolve,          only:evol
 use test,            only:testsuite
 implicit none
 integer :: nargs,i
 character(len=120) :: infile,logfile,evfile,dumpfile

 id = 0

 call init_mpi(id,nprocs)

 call set_io_unit_numbers
!
!--get name of run from the command line
!
 nargs = command_argument_count()
 if (nargs < 1) then
    if (id==master) then
       print "(a,/)",trim(tagline)
       print "(a)",' Usage: phantom infilename '
    endif
    call die
 endif
 call get_command_argument(1,infile)

 if (trim(infile)=='test') then
    call initialise()
    if (nargs >= 2) then
       do i=2,nargs
          call get_command_argument(i,infile)
          call testsuite(trim(infile),(i==2),(i==nargs))
       enddo
    else
       call testsuite('all',.true.,.true.)
    endif
 else
    if (index(infile,'.in')==0) then
       infile = trim(infile)//'.in'
    endif
    call startrun(infile,logfile,evfile,dumpfile)
    call evol(infile,logfile,evfile,dumpfile)
    if (id==master) call endrun()
 endif

 call finalise_mpi()

end program phantom
