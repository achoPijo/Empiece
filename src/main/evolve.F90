!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: evolve
!
!  DESCRIPTION:
!   Evolves the simulation through all timesteps
!   This subroutine contains the main timestepping loop and calls
!   the output routines at the appropriate times
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: analysis, centreofmass, dim, energies, evwrite,
!    externalforces, fileutils, forcing, inject, io, io_summary, mf_write,
!    mpiutils, options, part, ptmass, quitdump, readwrite_dumps,
!    readwrite_infile, sort_particles, step_lf_global, supertimestep,
!    timestep, timestep_ind, timestep_sts, timing
!+
!--------------------------------------------------------------------------
module evolve
 implicit none
 public :: evol

 private

contains

subroutine evol(infile,logfile,evfile,dumpfile)
 use io,               only:iprint,iwritein,id,master,iverbose,flush_warnings,nprocs,fatal
 use timestep,         only:time,tmax,dt,dtmax,nmax,nout,nsteps,dtextforce
 use evwrite,          only:write_evfile,write_evlog
 use energies,         only:etot,totmom,angtot,mdust,get_erot_com
 use dim,              only:calc_erot,maxvxyzu,mhd,use_dustfrac
 use fileutils,        only:getnextfilename
 use options,          only:nfulldump,twallmax,dtwallmax,nmaxdumps,iexternalforce,&
                            icooling,ieos,ipdv_heating,ishock_heating,iresistive_heating
 use readwrite_infile, only:write_infile
 use readwrite_dumps,  only:write_smalldump,write_fulldump
 use step_lf_global,   only:step
 use timing,           only:get_timings,print_time,timer,reset_timer,increment_timer,&
                            timer_dens,timer_force,timer_link
 use mpiutils,         only:reduce_mpi,reduceall_mpi
#ifdef SORT
 use sort_particles,   only:sort_part
#endif
#ifdef IND_TIMESTEPS
 use dim,              only:maxp
 use part,             only:maxphase,ibin,iphase
 use timestep_ind,     only:istepfrac,nbinmax,set_active_particles,update_time_per_bin,&
                            write_binsummary,change_nbinmax,nactive,nactivetot,maxbins,&
                            decrease_dtmax
 use timestep,         only:restartonshortest,dtdiff
 use timestep_sts,     only:sts_get_dtau_next,sts_init_step
 use io,               only:fatal,warning
 use step_lf_global,   only:init_step
#else
 use timestep,         only:dtforce,dtcourant,dterr
#endif
 use timestep_sts,     only: use_sts
 use supertimestep,    only: step_sts
#ifdef DRIVING
 use forcing,          only:write_forcingdump
#endif
#ifdef CORRECT_BULK_MOTION
 use centreofmass,     only:correct_bulk_motion
#endif
#ifdef MPI
 use part,             only:ideadhead,shuffle_part
#endif
#ifdef INJECT_PARTICLES
 use inject,           only:inject_particles
 use part,             only:npartoftype
#ifdef IND_TIMESTEPS
 use part,             only:twas
 use timestep_ind,     only:get_dt
#endif
#endif
#ifdef LIVE_ANALYSIS
 use analysis,         only:do_analysis
 use part,             only:xyzh,vxyzu,massoftype,igas
 use fileutils,        only:numfromfile
 use io,               only:ianalysis
#endif
 use part,             only:npart,nptmass,xyzh,vxyzu,fxyzu,fext,divcurlv,massoftype, &
                            xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,gravity,iboundary,npartoftype
 use quitdump,         only:quit
 use ptmass,           only:icreate_sinks,ptmass_create,ipart_rhomax,pt_write_sinkev
 use io_summary,       only:iosum_nreal,summary_counter,summary_printout,summary_printnow
 use externalforces,   only:iext_spiral
#ifdef MFLOW
 use mf_write,         only:mflow_write
#endif
#ifdef VMFLOW
use mf_write,          only:vmflow_write
#endif
#ifdef BINPOS
use mf_write,          only:binpos_write
#endif

 character(len=*), intent(in)    :: infile
 character(len=*), intent(inout) :: logfile,evfile,dumpfile
 integer         :: noutput,noutput_dtmax,nsteplast,ncount_fulldumps
 real            :: dtnew,dtlast,timecheck,etot_in,angtot_in,totmom_in,mdust_in
 real            :: tprint,tzero,dtmaxold
 real(kind=4)    :: t1,t2,tcpu1,tcpu2,tstart,tcpustart
 real(kind=4)    :: twalllast,tcpulast,twallperdump,twallused
#ifdef IND_TIMESTEPS
 integer         :: i,nalive,inbin,iamtypei
 integer(kind=1) :: nbinmaxprev
 integer(kind=8) :: nmovedtot,ntot
 real            :: tlast,fracactive,speedup,tcheck,dtau
 real(kind=4)    :: tall
 real(kind=4)    :: timeperbin(0:maxbins)
 logical         :: dt_changed
#ifdef INJECT_PARTICLES
 integer         :: iloop,npart_old
#endif
#else
 real            :: dtprint
 integer         :: nactive
 logical, parameter :: dt_changed = .false.
#endif
 logical         :: fulldump,abortrun,at_dump_time,update_tzero
 logical         :: should_conserve_energy,should_conserve_momentum,should_conserve_angmom
 logical         :: should_conserve_dustmass
 integer         :: nskip,nskipped,nevwrite_threshold,nskipped_sink,nsinkwrite_threshold
 type(timer)     :: timer_fromstart,timer_lastdump,timer_step,timer_ev,timer_io

 tprint    = 0.
 nsteps    = 0
 nsteplast = 0
 tzero     = time
 dtlast    = 0.
 etot_in   = etot
 angtot_in = angtot
 totmom_in = totmom
 mdust_in  = mdust
 should_conserve_energy = (maxvxyzu==4 .and. ieos==2 .and. icooling==0 .and. &
                           ipdv_heating==1 .and. ishock_heating==1 &
                           .and. (.not.mhd .or. iresistive_heating==1))
 should_conserve_momentum = (npartoftype(iboundary)==0)
 should_conserve_angmom   = (npartoftype(iboundary)==0)
 should_conserve_dustmass = use_dustfrac

 noutput          = 1
 noutput_dtmax    = 1
 ncount_fulldumps = 0
 tprint           = tzero + dtmax
 update_tzero     = .false.
#ifdef IND_TIMESTEPS
 istepfrac = 0
 tlast = tzero
 dt = dtmax/2**nbinmax
 nmovedtot = 0
 tall = 0.
 tcheck = time
 timeperbin(:) = 0.
 dt_changed = .false.
!
! first time through, move all particles on shortest timestep
! then allow them to gradually adjust levels.
! Keep boundary particles on level 0 since forces are never calculated
! and to prevent boundaries from limiting the timestep
!
 if (time < tiny(time) .or. restartonshortest) then
    !$omp parallel do schedule(static) private(i,iamtypei)
    do i=1,npart
       ibin(i) = nbinmax
       if (maxphase==maxp) then
          if (abs(iphase(i))==iboundary) ibin(i) = 0
       endif
    enddo
 endif
 if (use_sts) then
    call sts_get_dtau_next(dtau,dt,dtmax,dtdiff,nbinmax)
    call sts_init_step(npart,time,dtmax,dtau)
 else
    call init_step(npart,time,dtmax)
 endif

#else
 nskip = npart
 nactive = npart
 if (dt >= (tprint-time)) dt = tprint-time   ! reach tprint exactly
#endif
!
! threshold for writing to .ev file, to avoid repeatedly computing energies
! for all the particles which would add significantly to the cpu time
!
 nskipped = 0
 if (iexternalforce==iext_spiral) then
    nevwrite_threshold = int(4.99*npart) ! every 5 full steps
 else
    nevwrite_threshold = int(1.99*npart) ! every 2 full steps
 endif
 nskipped_sink = 0
 nsinkwrite_threshold  = int(0.99*npart)
!
! timing between dumps
!
 call get_timings(twalllast,tcpulast)
 tstart    = twalllast
 tcpustart = tcpulast
 call reset_timer(timer_fromstart,'all')
 call reset_timer(timer_lastdump,'last')
 call reset_timer(timer_step,'step')
 call reset_timer(timer_io,'write_dump')
 call reset_timer(timer_ev,'write_ev')
 call reset_timer(timer_dens,'density')
 call reset_timer(timer_force,'force')
 call reset_timer(timer_link,'link')

 call flush(iprint)
#ifdef LIVE_ANALYSIS
 if (id==master) then
    call do_analysis(dumpfile,numfromfile(dumpfile),xyzh,vxyzu, &
                     massoftype(igas),npart,time,ianalysis)
 endif
#endif
!
! --------------------- main loop ----------------------------------------
!
 timestepping: do while ((time < tmax).and.((nsteps < nmax) .or.  (nmax < 0)))

#ifdef INJECT_PARTICLES
    !
    ! injection of new particles into simulation
    !
#ifdef IND_TIMESTEPS
    npart_old=npart
#endif
    call inject_particles(time,dtlast,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,npart,npartoftype)
#ifdef IND_TIMESTEPS
    do iloop=npart_old+1,npart
       ibin(iloop) = nbinmax
       twas(iloop) = time + 0.5*get_dt(dtmax,ibin(iloop))
    enddo
#endif
#endif

    dtmaxold    = dtmax
#ifdef IND_TIMESTEPS
    istepfrac   = istepfrac + 1
    nbinmaxprev = nbinmax
    !--determine if dt needs to be decreased; if so, then this will be done
    !  in step the next time it is called;
    !  for global timestepping, this is called in the block where at_dump_time==.true.
    if (istepfrac==2**nbinmax) then
       twallperdump = timer_lastdump%wall
       call check_dtmax_for_decrease(iprint,dtmax,twallperdump,update_tzero)
    endif

    !--sanity check on istepfrac...
    if (istepfrac > 2**nbinmax) then
       write(iprint,*) 'ERROR: istepfrac = ',istepfrac,' / ',2**nbinmax
       call fatal('evolve','error in individual timesteps')
    endif

    !--flag particles as active or not for this timestep
    call set_active_particles(npart,nactive,nalive,iphase,ibin,xyzh)
    nactivetot = nactive
    ntot = nalive
    nskip = int(nactivetot)

    !--print summary of timestep bins
    if (iverbose >= 2) call write_binsummary(npart,nbinmax,dtmax,timeperbin,iphase,ibin,xyzh)
#endif

    if (gravity .and. icreate_sinks > 0 .and. ipart_rhomax > 0) then
       !
       ! creation of new sink particles
       !
       call ptmass_create(nptmass,npart,ipart_rhomax,xyzh,vxyzu,fxyzu,fext,divcurlv,&
                          massoftype,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,time)
    endif

    nsteps = nsteps + 1
!
!--evolve data for one timestep
!  for individual timesteps this is the shortest timestep
!
    call get_timings(t1,tcpu1)
!
    if ( use_sts ) then
       call step_sts(npart,nactive,time,dt,dtextforce,dtnew,iprint)
    else
       call step(npart,nactive,time,dt,dtextforce,dtnew)
    endif
    dtlast = dt

    !--timings for step call
    call get_timings(t2,tcpu2)
    call increment_timer(timer_step,t2-t1,tcpu2-tcpu1)
    call summary_counter(iosum_nreal,t2-t1)

#ifdef IND_TIMESTEPS
    tcheck = tcheck + dt

    !--update time in way that is free of round-off errors
    time = tlast + istepfrac/real(2**nbinmaxprev)*dtmaxold

    !--print efficiency of partial timestep
    if (id==master .and. iverbose >= 0 .and. ntot > 0) then
       if (nactivetot==ntot) then
          tall = t2-t1
       elseif (tall > 0.) then
          fracactive = nactivetot/real(ntot)
          speedup = (t2-t1)/tall
          if (iverbose >= 2) &
             write(iprint,"(1x,'(',3(a,f6.2,'%'),')')") &
                  'moved ',100.*fracactive,' of particles in ',100.*speedup, &
                  ' of time, efficiency = ',100.*fracactive/speedup
       endif
    endif
    call update_time_per_bin(tcpu2-tcpu1,istepfrac,nbinmaxprev,timeperbin,inbin)
    nmovedtot = nmovedtot + nactivetot

    !--check that time is as it should be, may indicate error in individual timestep routines
    if (abs(tcheck-time) > 1.e-4) call warning('evolve','time out of sync',var='error',val=abs(tcheck-time))

    if (id==master .and. (iverbose >= 1 .or. inbin <= 3)) &
       call print_dtlog_ind(iprint,istepfrac,2**nbinmaxprev,time,dt,nactivetot,tcpu2-tcpu1)

    !--if total number of bins has changed, adjust istepfrac and dt accordingly
    !  (ie., decrease or increase the timestep)
    if (nbinmax /= nbinmaxprev) then
       call change_nbinmax(nbinmax,nbinmaxprev,istepfrac,dtmax,dt)
       dt_changed = .true.
    endif

#else
    time = time + dt
!
!--set new timestep from Courant/forces condition
!
    ! constraint from time to next printout, must reach this exactly
    ! Following redefinitions are to avoid crashing if dtprint = 0 & to reach next output while avoiding round-off errors
    dtprint = min(tprint,tmax) - time + epsilon(dtmax)
    if (dtprint <= epsilon(dtmax) .or. dtprint >= (1.0-1e-8)*dtmax ) dtprint = dtmax + epsilon(dtmax)
    dt = min(dtforce,dtcourant,dterr,dtmax+epsilon(dtmax),dtprint)
!
!--write log every step (NB: must print after dt has been set in order to identify timestep constraint)
!
    if (id==master) call print_dtlog(iprint,time,dt,dtforce,dtcourant,dterr,dtmax,dtprint)
#endif

!    if (abs(dt) < 1e-8*dtmax) then
!       write(iprint,*) 'main loop: timestep too small, dt = ',dt
!       call quit   ! also writes dump file, safe here because called by all threads
!    endif


!   check that MPI threads are synchronised in time
    timecheck = reduceall_mpi('+',time)
    if (abs(timecheck/nprocs - time) > 1.e-13) then
       call fatal('evolve','time differs between MPI threads',var='time',val=timecheck/nprocs)
    endif
!
!--Update timer from last dump to see if dtmax needs to be reduced
!
    call get_timings(t2,tcpu2)
    call increment_timer(timer_lastdump,t2-t1,tcpu2-tcpu1)

!
!--Determine if this is the correct time to write to the data file
!
    at_dump_time = (time >= tmax).or.((mod(nsteps,nout)==0).and.(nout > 0)) &
                   .or.((nsteps >= nmax).and.(nmax >= 0))
#ifdef IND_TIMESTEPS
    if (istepfrac==2**nbinmax) at_dump_time = .true.
#else
    if (time >= tprint) at_dump_time = .true.
#endif
!
!--Calculate total energy etc and write to ev file
!  For individual timesteps, we do not want to do this every step, but we want
!  to do this as often as possible without a performance hit. The criteria
!  here is that it is done once >10% of particles (cumulatively) have been evolved.
!  That is, either >10% are being stepped, or e.g. 1% have moved 10 steps.
!  Perform this prior to writing the dump files so that diagnostic values calculated
!  in energies can be correctly included in the dumpfiles
!
    nskipped = nskipped + nskip
    if (nskipped >= nevwrite_threshold .or. at_dump_time .or. dt_changed) then
       ! update centre of mass for calculation of rotational energy (infrequent for optimisation)(
       if (at_dump_time .and. calc_erot) then
          call get_erot_com(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)
       endif
       nskipped = 0
       call get_timings(t1,tcpu1)
       call write_evfile(time,dt)
       if (should_conserve_momentum) call check_conservation_error(totmom,totmom_in,1.e-1,'linear momentum')
       if (should_conserve_angmom)   call check_conservation_error(angtot,angtot_in,1.e-1,'angular momentum')
       if (should_conserve_energy)   call check_conservation_error(etot,etot_in,1.e-1,'energy')
       if (should_conserve_dustmass) call check_conservation_error(mdust,mdust_in,1.e-1,'dust mass')
       !--write with the same ev file frequency also mass flux and binary position
#ifdef MFLOW
       call mflow_write(time,dt)
#endif
#ifdef VMFLOW
       call vmflow_write(time,dt)
#endif
#ifdef BINPOS
       call binpos_write(time,dt)
#endif
       !--timings for step call
       call get_timings(t2,tcpu2)
       call increment_timer(timer_ev,t2-t1,tcpu2-tcpu1)
    endif
!-- Print out the sink particle properties & reset dt_changed
    nskipped_sink = nskipped_sink + nskip
    if (nskipped_sink >= nsinkwrite_threshold .or. at_dump_time .or. dt_changed) then
       nskipped_sink = 0
       call pt_write_sinkev(nptmass,time,xyzmh_ptmass,vxyz_ptmass)
#ifdef IND_TIMESTEPS
       dt_changed = .false.
#endif
    endif
!
!--write to data file if time is right
!
    if (at_dump_time) then
       !--modify evfile and logfile names with new number
       if (noutput==1) then
          evfile  = getnextfilename(evfile)
          logfile = getnextfilename(logfile)
       endif
       dumpfile = getnextfilename(dumpfile)

#ifdef MPI
       !--do not dump dead particles into dump files
       if (ideadhead > 0) call shuffle_part(npart)
#endif
!
!--get timings since last dump and overall code scaling
!  (get these before writing the dump so we can check whether or not we
!   need to write a full dump based on the wall time;
!   move timer_lastdump outside at_dump_time block so that dtmax can
!   be reduced it too long between dumps)
!
       call increment_timer(timer_fromstart,t2-tstart,tcpu2-tcpustart)
       !call increment_timer(timer_lastdump,t2-twalllast,tcpu2-tcpulast)
       timer_fromstart%cpu  = reduce_mpi('+',timer_fromstart%cpu)
       timer_lastdump%cpu   = reduce_mpi('+',timer_lastdump%cpu)
       timer_step%cpu = reduce_mpi('+',timer_step%cpu)
       timer_ev%cpu   = reduce_mpi('+',timer_ev%cpu)

       fulldump = (mod(noutput,nfulldump)==0)
!
!--if max wall time is set (> 1 sec) stop the run at the last full dump
!  that will fit into the walltime constraint, based on the wall time between
!  the last two dumps added to the current total walltime used.  The factor of three for
!  chaning to full dumps is to account for the possibility that the next step will take longer.
!  If we are about to write a small dump but it looks like we won't make the next dump,
!  dump a full dump instead and stop the run
!
       abortrun = .false.
       if (twallmax > 1.) then
          twallused    = timer_fromstart%wall
          twallperdump = timer_lastdump%wall
          if (fulldump) then
             if ((twallused + abs(nfulldump)*twallperdump) > twallmax) then
                abortrun = .true.
             endif
          else
             if ((twallused + 3.0*twallperdump) > twallmax) then
                fulldump = .true.
                write(iprint,"(1x,a)") '>> PROMOTING DUMP TO FULL DUMP BASED ON WALL TIME CONSTRAINTS... '
                nfulldump = 1  !  also set all future dumps to be full dumps (otherwise gets confusing)
                if ((twallused + twallperdump) > twallmax) abortrun = .true.
             endif
          endif
       endif
!
!--if max wall time between dumps is set (> 1 sec), switch to writing full dumps for all subsequent dumps
!
       if (dtwallmax > 1.0 .and. .not.fulldump) then
          twallperdump = timer_lastdump%wall
          if (twallperdump > dtwallmax) then
            fulldump = .true.
            write(iprint,"(1x,a)") '>> PROMOTING DUMP TO FULL DUMP BASED ON DT WALL TIME CONSTRAINTS... '
            nfulldump = 1  !  also set all future dumps to be full dumps (otherwise gets confusing)
          endif
       endif
!
!--Promote to full dump if this is the final dump
!
       if ( (time >= tmax) .or. ( (nmax > 0) .and. (nsteps >= nmax) ) ) fulldump = .true.

#ifndef IND_TIMESTEPS
!
!--Decrease dtmax if requested (this is done in step for individual timesteps)
       twallperdump = timer_lastdump%wall
       call check_dtmax_for_decrease(iprint,dtmax,twallperdump,update_tzero,.true.)
#endif
!
!--flush any buffered warnings to the log file
!
       if (id==master) call flush_warnings()
!
!--write dump file
!
       call get_timings(t1,tcpu1)
       if (fulldump) then
          call write_fulldump(time,dumpfile)
          if (id==master) then
             call write_infile(infile,logfile,evfile,dumpfile,iwritein,iprint)
#ifdef DRIVING
             call write_forcingdump(time,dumpfile)
#endif
          endif
          ncount_fulldumps = ncount_fulldumps + 1

#ifndef IND_TIMESTEPS
#ifdef SORT
          if (time < tmax) call sort_part()
#endif
#endif
       else
          call write_smalldump(time,dumpfile)
       endif
       call get_timings(t2,tcpu2)
       call increment_timer(timer_io,t2-t1,tcpu2-tcpu1)
       timer_io%cpu   = reduce_mpi('+',timer_io%cpu)

#ifdef LIVE_ANALYSIS
       if (id==master) then
          call do_analysis(dumpfile,numfromfile(dumpfile),xyzh,vxyzu, &
                           massoftype(igas),npart,time,ianalysis)
       endif
#endif

       if (id==master) then
          call print_timinginfo(iprint,nsteps,nsteplast,timer_fromstart,timer_lastdump,timer_step,timer_ev,timer_io,&
                                             timer_dens,timer_force,timer_link)
          !--Write out summary to log file
          call summary_printout(iprint,nptmass)
       endif
#ifdef IND_TIMESTEPS
       !--print summary of timestep bins
       if (iverbose >= 0 .and. id==master .and. abs(tall) > tiny(tall) .and. ntot > 0) then
          fracactive = nmovedtot/real(ntot)
          speedup = timer_lastdump%wall/(tall + tiny(tall))
          write(iprint,"(/,a,f6.2,'%')") ' IND TIMESTEPS efficiency = ',100.*fracactive/speedup
          if (iverbose >= 1) then
             write(iprint,"(a,1pe14.2,'s')") '  wall time per particle (last full step) : ',tall/real(ntot)
             write(iprint,"(a,1pe14.2,'s')") '  wall time per particle (ave. all steps) : ',timer_lastdump%wall/real(nmovedtot)
          endif
       endif
       if (iverbose >= 0) then
          call write_binsummary(npart,nbinmax,dtmax,timeperbin,iphase,ibin,xyzh)
          timeperbin(:) = 0.
       endif
       tlast = tprint
       istepfrac = 0
       nmovedtot = 0
#endif
       !
       !--if twallmax > 1s stop the run at the last full dump that will fit into the walltime constraint,
       !  based on the wall time between the last two dumps added to the current total walltime used.
       !
       if (abortrun) then
          call print_time(t2-tstart,'>> WALL TIME = ',iprint)
          call print_time(twallmax,'>> NEXT DUMP WILL TRIP OVER MAX WALL TIME: ',iprint)
          write(iprint,"(1x,a)") '>> ABORTING... '
          return
       endif

       if (nmaxdumps > 0 .and. ncount_fulldumps >= nmaxdumps) then
          write(iprint,"(a)") '>> reached maximum number of full dumps as specified in input file, stopping...'
          return
       endif

       twalllast = t2
       tcpulast = tcpu2
       call reset_timer(timer_fromstart)
       call reset_timer(timer_lastdump)
       call reset_timer(timer_step)
       call reset_timer(timer_ev)
       call reset_timer(timer_io)
       call reset_timer(timer_dens)
       call reset_timer(timer_force)
       call reset_timer(timer_link)

       noutput_dtmax = noutput_dtmax + 1
       noutput       = noutput + 1
       tprint        = tzero + noutput_dtmax*dtmaxold
       nsteplast     = nsteps
       if (update_tzero) then
          tzero         = tprint - dtmaxold
          tprint        = tzero  + dtmax
          noutput_dtmax = 1
       endif

       !  print summary of energies and other useful values to the log file
       if (id==master) call write_evlog(iprint)

    endif

#ifdef CORRECT_BULK_MOTION
    call correct_bulk_motion()
#endif

#ifndef IND_TIMESTEPS
#ifdef INJECT_PARTICLES
    !call dt_spheres_sync(time, dt)
#endif
#endif

    if (iverbose >= 1 .and. id==master) write(iprint,*)
    call flush(iprint)
    !--Write out log file prematurely (if requested based upon nstep, walltime)
    if ( summary_printnow() ) call summary_printout(iprint,nptmass)


 enddo timestepping

!------------------------------------------------------------------------

 return
end subroutine evol

!----------------------------------------------------------------
!+
!  This will determine if dtmax needs to be decreased
!  This subroutine is called at different times depending on
!  whether individual or global timesteps are being used.
!+
!----------------------------------------------------------------
subroutine check_dtmax_for_decrease(iprint,dtmax,twallperdump,update_tzero,update_dtmax)
 use timestep, only: dtmax_rat,dtmax_rat0,dtwall_dtthresh, &
                     mod_dtmax,mod_dtmax_now,mod_dtmax_in_step
 integer,      intent(in)    :: iprint
 real,         intent(inout) :: dtmax
 real(kind=4), intent(in)    :: twallperdump
 logical,      intent(out)   :: update_tzero
 logical,      intent(in), optional :: update_dtmax

 ! initialise variables
 dtmax_rat    = 0
 update_tzero = .false.

 ! modify dtmax based upon wall time constraint, if requested
 if ( dtwall_dtthresh > 0.0 ) then
    if (twallperdump > dtwall_dtthresh) then
       dtmax_rat = int(2**(int(log(real(twallperdump/dtwall_dtthresh))/log(2.0))+1))
       write(iprint,'(1x,a,2(es10.3,a))') &
          "modifying dtmax: ",dtmax," --> ",dtmax/dtmax_rat," due to wall time constraint"
       dtmax_rat0 = max(0,dtmax_rat0-dtmax_rat)
    endif
 endif

 ! modify dtmax based upon density, if requested, and print info to the logfile
 if (mod_dtmax_now) then
    dtmax_rat         = dtmax_rat0 + dtmax_rat
    mod_dtmax_now     = .false. ! reset variable
    mod_dtmax         = .false. ! prevent any additional modifications of dtmax
    write(iprint,'(1x,a,2(es10.3,a))') &
       "modifying dtmax: ",dtmax," --> ",dtmax/dtmax_rat," due to user requested decrease"
 endif

 ! dtmax will be modified; update all required logicals such that this will happen
 if (dtmax_rat > 1) then
    update_tzero      = .true.  ! to reset tzero and tlast if individual dt
    mod_dtmax_in_step = .true.  ! to tell step.f to decrease dtmax on its next call for individual dt
    if (present(update_dtmax)) then
       if (update_dtmax) dtmax = dtmax/dtmax_rat
    endif
 endif

end subroutine check_dtmax_for_decrease

!----------------------------------------------------------------
!+
!  routine to check conservation errors during the calculation
!  and stop if it is too large
!+
!----------------------------------------------------------------
subroutine check_conservation_error(val,ref,tol,label)
 use io, only:error,fatal,iverbose
 real, intent(in) :: val,ref,tol
 character(len=*), intent(in) :: label
 real :: err
 character(len=20) :: string

 if (abs(ref) > 1.e-6) then
    err = abs(val - ref)/abs(ref)
 else
    err = abs(val - ref)
 endif
 if (err > tol) then
    call error('evolve','Large error in '//trim(label)//' conservation ',var='err',val=err)
    call get_environment_variable('I_WILL_NOT_PUBLISH_CRAP',string)
    if (.not. (trim(string)=='yes')) then
       print "(2(/,a))",' You can ignore this error and continue by setting the ',&
                        ' environment variable I_WILL_NOT_PUBLISH_CRAP=yes to continue'
       call fatal('evolve',' Conservation errors too large to continue simulation')
    endif
 else
    if (iverbose >= 2) print "(a,es10.3)",trim(label)//' error is ',err
 endif

end subroutine check_conservation_error

!----------------------------------------------------------------
!+
!  routine to print out the timing information at each full dump
!+
!----------------------------------------------------------------
subroutine print_timinginfo(iprint,nsteps,nsteplast,&
           timer_fromstart,timer_lastdump,timer_step,timer_ev,timer_io,&
           timer_dens,timer_force,timer_link)
 use io,     only:formatreal
 use timing, only:timer
 integer,      intent(in) :: iprint,nsteps,nsteplast
 type(timer),  intent(in) :: timer_fromstart,timer_lastdump,timer_step,timer_ev,timer_io,&
                             timer_dens,timer_force,timer_link
 real                     :: dfrac,fracinstep
 real(kind=4)             :: time_fullstep
 character(len=20)        :: string,string1,string2,string3

 !write(iprint,*) 'Timings:'
 write(string,"(i12)") nsteps
 call formatreal(real(timer_fromstart%wall),string1)
 call formatreal(real(timer_fromstart%cpu),string2)
 call formatreal(real(timer_fromstart%cpu/(timer_fromstart%wall+epsilon(0._4))),string3)
 write(iprint,"(1x,'Since code start: ',a,' timesteps, wall: ',a,'s cpu: ',a,'s cpu/wall: ',a)") &
       trim(adjustl(string)),trim(string1),trim(string2),trim(string3)

 write(string,"(i12)") nsteps-nsteplast
 call formatreal(real(timer_lastdump%wall),string1)
 call formatreal(real(timer_lastdump%cpu),string2)
 call formatreal(real(timer_lastdump%cpu/(timer_lastdump%wall+epsilon(0._4))),string3)
 write(iprint,"(1x,'Since last dump : ',a,' timesteps, wall: ',a,'s cpu: ',a,'s cpu/wall: ',a)") &
       trim(adjustl(string)),trim(string1),trim(string2),trim(string3)

 time_fullstep = timer_lastdump%wall + timer_io%wall
 write(iprint,"(/,16x,a)") ' wall        cpu    cpu/wall   frac'
 call print_timer(iprint,timer_step%label,timer_step, time_fullstep)
 call print_timer(iprint,"step (force)",  timer_force,time_fullstep)
 call print_timer(iprint,"step (dens) ",  timer_dens, time_fullstep)
 call print_timer(iprint,"step (link) ",  timer_link, time_fullstep)
 call print_timer(iprint,timer_ev%label,  timer_ev,   time_fullstep)
 call print_timer(iprint,timer_io%label,  timer_io,   time_fullstep)

 dfrac = 1./(timer_lastdump%wall + epsilon(0._4))
 fracinstep = timer_step%wall*dfrac
 if (fracinstep < 0.99) then
    write(iprint,"(1x,a,f6.2,a)") 'WARNING: spending ',100.*(1.-fracinstep),'% of time in routines EXTERNAL to step'
 endif

end subroutine print_timinginfo

!-----------------------------------------------
!+
!  Pretty-print timing entries in a nice table
!+
!-----------------------------------------------
subroutine print_timer(lu,label,my_timer,time_between_dumps)
 use timing, only:timer
 integer,          intent(in) :: lu
 real(kind=4),     intent(in) :: time_between_dumps
 character(len=*), intent(in) :: label
 type(timer),      intent(in) :: my_timer

 if (my_timer%wall > epsilon(0._4)) then
    if (time_between_dumps > 7200.0) then
       write(lu,"(1x,a12,':',f7.2,'h   ',f7.2,'h   ',f6.2,'   ',f6.2,'%')")  &
            label,my_timer%wall/3600.,my_timer%cpu/3600.,my_timer%cpu/my_timer%wall, &
            my_timer%wall/time_between_dumps*100.
    else if (time_between_dumps > 120.0) then
       write(lu,"(1x,a12,':',f7.2,'min ',f7.2,'min ',f6.2,'   ',f6.2,'%')")  &
            label,my_timer%wall/60.,my_timer%cpu/60.,my_timer%cpu/my_timer%wall, &
            my_timer%wall/time_between_dumps*100.
    else
       write(lu,"(1x,a12,':',f7.2,'s   ',f7.2,'s   ',f6.2,'   ',f6.2,'%')")  &
            label,my_timer%wall,my_timer%cpu,my_timer%cpu/my_timer%wall, &
            my_timer%wall/time_between_dumps*100.
    endif
 endif

end subroutine print_timer

#ifdef IND_TIMESTEPS
!-----------------------------------------------------------------
!+
!  routine to print out the timestep information to the log file
!  this version handles individual timesteps
!+
!-----------------------------------------------------------------
subroutine print_dtlog_ind(iprint,ifrac,nfrac,time,dt,nactive,tcpu)
 use io, only:formatreal,formatint
 integer,         intent(in) :: iprint,ifrac,nfrac
 real,            intent(in) :: time,dt
 integer(kind=8), intent(in) :: nactive
 real(kind=4),    intent(in) :: tcpu
 character(len=120) :: string
 character(len=14) :: tmp

 call formatint(ifrac,tmp)
 string = '> step '//tmp
 call formatint(nfrac,tmp)
 string = trim(string)//' / '//tmp
 write(tmp,"(g14.7)") time
 string = trim(string)//' t = '//trim(adjustl(tmp))
 call formatreal(dt,tmp)
 string = trim(string)//' dt = '//tmp
 call formatint(nactive,tmp)
 string = trim(string)//' moved '//tmp
 call formatreal(real(tcpu),tmp)
 string = trim(string)//' in '//trim(tmp)//' cpu-s <'

 write(iprint,"(a)") trim(string)
! write(iprint,5) ifrac,2**nbinmaxprev,time,dt,nactivetot,tcpu2-tcpu1
!5   format('> step ',i6,' /',i6,2x,'t = ',es14.7,1x,'dt = ',es10.3,' moved ',i10,' in ',f8.2,' cpu-s <')

end subroutine print_dtlog_ind
#else
!-----------------------------------------------------------------
!+
!  routine to print out the timestep information to the log file
!+
!-----------------------------------------------------------------
subroutine print_dtlog(iprint,time,dt,dtforce,dtcourant,dterr,dtmax,dtprint)
 integer, intent(in) :: iprint
 real,    intent(in) :: time,dt,dtforce,dtcourant,dterr,dtmax
 real,    intent(in), optional :: dtprint

 if (abs(dt-dtforce) < tiny(dt)) then
    write(iprint,10) time,dt,'(force)'
 elseif (abs(dt-dtcourant) < tiny(dt)) then
    write(iprint,10) time,dt,'(courant)'
 elseif (abs(dt-dterr) < tiny(dt)) then
    write(iprint,10) time,dt,'(tolv)'
 elseif (abs(dt-dtmax) <= epsilon(dt)) then
    write(iprint,10) time,dt,'(dtmax)'
 elseif (present(dtprint) .and. abs(dt-dtprint) < tiny(dt)) then
    write(iprint,10) time,dt,'(dtprint)'
 else
    !print*,dt,dtforce,dtcourant,dterr,dtmax
    write(iprint,10) time,dt,'(unknown)'
 endif
10 format(' t = ',g12.5,' dt = ',es10.3,1x,a)

end subroutine print_dtlog
#endif

end module evolve
