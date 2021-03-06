!----------------------------------------------------------------------!
!                               N I C I L                              !                                  
!           Non-Ideal mhd Coefficient and Ionisation Library           !
!         Example programme to test various parameters of NICIL        !
!                                                                      !
!                 Copyright (c) 2015-2016 James Wurster                !
!        See LICENCE file for usage and distribution conditions        !
!----------------------------------------------------------------------!
!+
! This is a test programme to calculate grain charge, electron number 
! densities and the non-ideal MHD coefficients for a given range of 
! densities.  Results from this programme can be directly compared to 
! the author's output found in the data folder.
!+
!----------------------------------------------------------------------!
program nicil_ex_eta
 use nicil,  only:nicil_initialise,nicil_get_ion_n,nicil_get_eta,nicil_translate_error
 use nicil,  only:nelements_max,nelements,nlevels
 use etasup, only:Bconst,which_Bfield,get_Bfield_code,get_temperature
 use etasup, only:write_data_header,write_data_to_file,fatal
 use etasup, only:iprint,iprintrho,iprintbaro,iprinttemp,iprintwarn
 implicit none
 !--Input parameters for fixed temperature (varying density)
 integer, parameter  :: nlogd           = 10000           ! Number of density points to test
 real,    parameter  :: logrho_min      = -22.0           ! min(log(rho)) to test [g/cm^3]
 real,    parameter  :: logrho_max      =   0.5           ! max(log(rho)) to test [g/cm^3]
 real,    parameter  :: temperature     =  30.0           ! Temperture of gas [K]
 !--Input parameters for fixed density (varying temperature)
 integer, parameter  :: nlogT           = 10000           ! Number of temperature points to test
 real,    parameter  :: temp_min        = 1.0d1           ! minimum temperature to test [K]
 real,    parameter  :: temp_max        = 2.0d5           ! maximum temperature to test [K]
 real,    parameter  :: rho_in          = 1.0d-13         ! density to test [g/cm^3]
 !--Input parameters for using a fixed magnetic field
 real,    parameter  :: B_input         = 1.583d-4        ! user input magnetic field [G]
 logical             :: use_input_B     = .false.         ! use user's B (true) or use pre-defined function for B (false)
 !--Input parameters (other)
 real,    parameter  :: mu0             =  2.38095236     ! Mean molecular mass; calculated with default values in NICIL
 !--Physical Constants
 real,    parameter  :: fourpi          =  12.566370614d0 ! 4pi
 real,    parameter  :: kboltz          = 1.38066d-16     ! Boltzmann constant  [erg/K]
 real,    parameter  :: mass_proton_cgs = 1.67262158d-24  ! Proton mass [g]
 real,    parameter  :: cgsmu0          = fourpi          ! Vacuum permeability [cm/s]
 !--Local variables
 integer             :: i,j,ierr
 real                :: utime,udist,umass
 real                :: unit_density,unit_ndensity,unit_charge,unit_Bfield,unit_eta,mump1
 real                :: T,rho,Bfield,temp
 real                :: dlogrho,dlogT
 real                :: eta_ohm,eta_hall,eta_ambi
 real                :: n_R(4,nlogd),n_electronT(nlogd)
 real                :: data_out(17+nelements_max*nlevels-3)
 character(len=1)    :: Btype
 !
 !
 !--Initialise parameters
 !  Code Units
 udist         = 1.000d16                                 ! = 1 code length unit 
 umass         = 1.989d33                                 ! = 1 M_sun = 1 code mass unit
 utime         = 8.681d10                                 ! = 1 code time unit (Chosen such that G=1)
 unit_density  = umass/udist**3                           ! = 1 code density unit
 unit_ndensity = 1.0/udist**3                             ! = 1 code number density unit
 unit_charge   = sqrt(umass*udist/cgsmu0)                 ! = 1 code charge unit
 unit_Bfield   = umass/(utime*unit_charge)                ! = 1 code magentic field unit *may be defined differently in user's code*
 unit_eta      = udist**2/utime                           ! = 1 code eta unit
 mump1         = 1.0/(mu0*mass_proton_cgs)                ! inverse mean mass [CGS]
 !  Zero Arrays
 n_R           = 0.0
 n_electronT   = 0.0
 !  Set pseudo-grid spacing
 dlogrho       = (logrho_max  - logrho_min )/nlogd        ! spacing of density points
 dlogT         = (log10(temp_max) - log10(temp_min))/nlogT! spacing of temperature points
 !  Rename the input magnetic field
 Bconst        = B_input
 !
 !--Call the initialisation routine.  
 !  Required only once to calculate and save all the values required for this library
 call nicil_initialise(utime,udist,umass,unit_Bfield,ierr,iprint,iprintwarn)
 if (ierr/=0) call fatal(ierr)                            ! Abort programme if there are errors in the setup
 !
 !--Open output files & write header
 open(unit=iprintrho, file="data/eta_density.dat")
 open(unit=iprintbaro,file="data/eta_barotropic.dat")
 open(unit=iprinttemp,file="data/eta_temperature.dat")
 open(unit=iprintwarn,file="data/eta_warning.log")
 do i = iprintrho,iprinttemp
   call write_data_header(i)
 end do
 write(iprintwarn,'(a)') "NICIL: ETA TEST: WARNINGS LOG"
 write(iprint,'(a)')     "NICIL: ETA TEST"
 !
 !--This will calculate the coefficients, grain charge and electron number density
 !  for a range of densities at a fixed temperatures;
 !  magnetic field and sound speed are related to the density by a given prescription
 write(iprintwarn,'(a)') "NICIL: ETA TEST: properties vs number density (with constant T)"
 Btype = which_Bfield(use_input_B,"P")
 do i = 1,nlogd
   rho    = 10**(logrho_min +(i-1)*dlogrho )             ! density to test [CGS]
   Bfield = get_Bfield_code(rho*mump1,unit_Bfield,Btype) ! Magnetic field [code units]
   rho    = rho/unit_density                             ! density to test [code units]
   !
   ! Call NICIL to get the eta coefficients.
   ! data_out is an optional out variable that we will pass through to track the values
   call nicil_get_ion_n(rho,temperature,n_R(1:4,i),n_electronT(i),ierr)
   if (ierr/=0) then
     call nicil_translate_error(ierr)
     call fatal(ierr,rho*unit_density,temperature)
   end if
   call nicil_get_eta(eta_ohm,eta_hall,eta_ambi,Bfield,rho,temperature,n_R(1:4,i),n_electronT(i),ierr,data_out)
   if (ierr/=0) then
     call nicil_translate_error(ierr)
     call fatal(ierr,rho*unit_density,temperature)
   end if
   !
   ! Write values to file for testing purposes
   call write_data_to_file(iprintrho,rho*unit_density,temperature,Bfield,eta_ohm,eta_hall,eta_ambi &
                          ,data_out,unit_eta,unit_Bfield,unit_density,unit_ndensity)
   !
 end do
 write(iprint,'(a)') "NICIL: ETA TEST: properties vs number density (with constant T) written to data/eta_density.dat"
 !
 !--This will calculate the coefficients, grain charge and electron number density
 !  assuming a barotropic equation of state
 write(iprintwarn,'(a)') "NICIL: ETA TEST: properties vs {number density, temperature} using barotropic EOS"
 n_R         = 0.0                                     ! re-initialise the array since this is a new simulation
 n_electronT = 0.0                                     ! re-initialise the array since this is a new simulation
 Btype       = which_Bfield(use_input_B,"U")
 do i = 1,nlogd
   rho    = 10**(logrho_min +(i-1)*dlogrho )             ! density to test [CGS]
   Bfield = get_Bfield_code(rho*mump1,unit_Bfield,Btype) ! Magnetic field [code units]
   temp   = get_temperature(rho*mump1)                   ! temperature
   rho    = rho/unit_density                             ! density to test [code units]
   !
   ! Call NICIL to get the eta coefficients.
   ! data_out is an optional out variable that we will pass through to track the values
   call nicil_get_ion_n(rho,temp,n_R(1:4,i),n_electronT(i),ierr)
   if (ierr/=0) then
     call nicil_translate_error(ierr)
     call fatal(ierr,rho*unit_density,temperature)
   end if
   call nicil_get_eta(eta_ohm,eta_hall,eta_ambi,Bfield,rho,temp,n_R(1:4,i),n_electronT(i),ierr,data_out)
   if (ierr/=0) then
     call nicil_translate_error(ierr)
     call fatal(ierr,rho*unit_density,temperature)
   end if
   !
   ! Write values to file for testing purposes
   call write_data_to_file(iprintbaro,rho*unit_density,temp,Bfield,eta_ohm,eta_hall,eta_ambi &
                          ,data_out,unit_eta,unit_Bfield,unit_density,unit_ndensity)
   !
 end do
 write(iprint,'(a)') &
   "NICIL: ETA TEST: properties vs {number density, temperature} using barotropic EOS written to data/eta_barotropic.dat"
 !
 !--This will calculate the coefficients, grain charge and electron number density
 !  for a range of temperatures at a fixed density and magnetic field strength
 write(iprintwarn,'(a)') "NICIL: ETA TEST: properties vs temperature (with constant density)"
 Btype       = which_Bfield(use_input_B,"P")
 Bfield      = get_Bfield_code(rho_in*mump1,unit_Bfield,Btype) ! Magnetic field [code units]
 rho         = rho_in/unit_density                             ! density to test [code units]
 n_R         = 0.0                                             ! re-initialise the array since this is a new simulation
 n_electronT = 0.0                                             ! re-initialise the array since this is a new simulation
 do i = 1,nlogT
   temp = 10**(log10(temp_min) +(i-1)*dlogT )                  ! temperature to test [CGS]
   !
   ! Call NICIL to get the eta coefficients.
   ! data_out is an optional out variable that we will pass through to track the values
   call nicil_get_ion_n(rho,temp,n_R(1:4,i),n_electronT(i),ierr)
   if (ierr/=0) then
     call nicil_translate_error(ierr)
     call fatal(ierr,rho*unit_density,temp)
   end if
   call nicil_get_eta(eta_ohm,eta_hall,eta_ambi,Bfield,rho,temp,n_R(1:4,i),n_electronT(i),ierr,data_out)
   if (ierr/=0) then
     call nicil_translate_error(ierr)
     call fatal(ierr,rho*unit_density,temp)
   end if
   !
   ! Write values to file for testing purposes;
   call write_data_to_file(iprinttemp,rho_in,temp,Bfield,eta_ohm,eta_hall,eta_ambi &
                          ,data_out,unit_eta,unit_Bfield,unit_density,unit_ndensity)
   !
 end do
 write(iprint,'(a)') "NICIL: ETA TEST: properties vs temperature (with constant density) written to data/eta_temperature.dat"
 !
 do i = iprintrho,iprintwarn
   close(i)
 end do
 !
!----------------------------------------------------------------------!
end program nicil_ex_eta
