      SUBROUTINE degenerate2(p)
!========================================================================
!     THIS   SUBROUTINE   ALLOWS   TO  SELECT   BETWEEN  THE DIFFERENT  
!     EQUATIONS  OF  STATE  OF  THE  CODE. IT'S IMPORTANT  TO  REALIZE  
!     THAT  EVERY  EOS  SUBROUTINE MUST PROVIDE (AT LEAST):
!
!     * TEMPERATURE, PRESSURE, Cv AND SPEED OF SOUND
!
!     Last revision: 30/April/2016
!=========================================================================
!
!--Load modules
!
      USE mod_essentials
      USE mod_parameters, ONLY : uden, unm, uen, up, uv, unl, tmin, tmax
      USE mod_commons,    ONLY : rho, xss, aion, zion, vxyzut, press,   &
                                 css, cvs, dPdT, cps, error, partype
!
!--Force to declare EVERYTHING
!
      IMPLICIT NONE
!
!--Helmholtz EOS definitions
!
      INCLUDE 'vector_eos.dek'
!
!--I/O variables
!
      INTEGER, INTENT(IN) :: p
!
!--Local variables
!
      REAL(8), DIMENSION(1) :: ewant_row
      REAL(8), PARAMETER :: eos_tol=1.0d-8
      REAL(8)  :: abar, zbar, tnew, densp, temp, ewantp, errorp, etotp, &
                  detp, cvp, rel
      INTEGER  :: k, newton, eosflag
      INTEGER, PARAMETER :: max_newton=50
      LOGICAL  :: eosmal
!
!--Set initial data for the EOS
!
      IF (vxyzut(5,p) > tmax) vxyzut(5,p) = tmax
      IF (vxyzut(5,p) < tmin) vxyzut(5,p) = tmin
!
      ewant_row(1) = vxyzut(4,p)*(uen/unm) ! EOS works in cgs units!!!
      temp_row(1)  = vxyzut(5,p)
      den_row(1)   = rho(p)*uden    ! EOS works in cgs units!!!
!
      abar = 0.0
      zbar = 0.0
      DO k=1,nel-1
         abar = abar + xss(k,p)/aion(k)
         zbar = zbar + xss(k,p)*zion(k)/aion(k)
      ENDDO
      abar_row(1) = 1.0/abar
      zbar_row(1) = zbar/abar
!
!--Set the number of elements sent into the EOS
!
      jlo_eos = 1
      jhi_eos = 1
!
!--Call EOS to obtain a first estimate for energy, pressure and derivatives
!
      CALL helmeos
      IF ((rank == MASTER).AND.(eosfail.EQV..true.)) PRINT*,'EOScorr fail'
!
!--Finish EOS in case of a relaxation
!
      IF ((RELFLAG.EQV..true.).OR.(vxyzut(5,p) == tmin).OR.             &
          (vxyzut(5,p) == tmax)) THEN
         vxyzut(4,p) = etot_row(1)*(unm/uen)
         press(p)    = ptot_row(1)/up
         css(p)      = cs_row(1)/uv
         cvs(p)      = cv_row(1)*(unm/uen)
         dPdT(p)     = dpt_row(1)*(unl**3/uen)
         cps(p)      = cp_row(1)/(uen/unm)
         RETURN
      ENDIF
!
!--Create the initial condition. Calculate new temperature
!
      tnew = temp_row(1) - (etot_row(1) - ewant_row(1))/det_row(1)
!
!--Calculate the error
!
      errorp = ABS(tnew-temp_row(1))/temp_row(1)
!
!--Store new temperature
!
      temp_row(1) = tnew
!
!--If freezing, keep temperature inside EOS tables
!
      IF (temp_row(1) < tmin) THEN
          temp_row(1) = tmin
          errorp = 0.1d0*eos_tol
      ENDIF
!
!--Loop over particles doing the Newton-Raphson iteration
!
      newton = 0
      DO WHILE ((errorp > eos_tol).AND.(newton < max_newton))
         newton = newton + 1
!
!--Call EOS
!
         jlo_eos = 1
         jhi_eos = 1
         CALL helmeos
         IF ((rank == MASTER).AND.(eosmal.EQV..true.)) PRINT*,'EOS in   &
             particle: ',p,' failed'
!
!--New temperature
!
         tnew = temp_row(1) - (etot_row(1)-ewant_row(1))/det_row(1)
!
!--Do no allow temperature to change more than an order of magnitude per
!  iteration
!
         IF (tnew > 10.d0*temp_row(1)) tnew = 10.0d0*temp_row(1)
         IF (tnew < 0.1d0*temp_row(1)) tnew = 0.1d0*temp_row(1)
!
!--Calculate the error
!
         errorp = ABS(tnew-temp_row(1))/temp_row(1)
! 
!--Store new temperature
!
         temp_row(1) = tnew
!
!--If freezing, keep temperature inside EOS tables
!
         IF (temp_row(1) < tmin) THEN
             temp_row(1) = tmin
             errorp      = 0.1d0*eos_tol
         ENDIF
!
!--If too hot, keep temperature inside EOS tables
!
         IF (temp_row(1) > tmax) THEN
             temp_row(1) = tmax
             errorp      = 0.1d0*eos_tol
         ENDIF
30    ENDDO    ! End Newton-Raphson loop
!
!--If the Newton-Rapshon fails to find a valid temperature, keep it 
!  constant. Check also if temperature and the temperature predicted
!  from the internal energy differe more than a 5%
      rel = DABS(temp_row(1)-vxyzut(5,p))/vxyzut(5,p)
      IF ((newton >= max_newton).OR.(rel > 0.05)) THEN
          eosflag = 2
      ELSE
          eosflag = 1
      END IF
      !IF (partype(p).EQ.1) eosflag=1
!
!--If eosflag=1 we store temp_row, whereas for eosflag=2 we store 
!  etot_row
!
      IF (eosflag == 2) THEN  ! Temperature as input
          vxyzut(4,p) = etot_row(1)*(unm/uen)
          error(p)    = 0.0d0
      ELSE                    ! Internal energy as input
          vxyzut(5,p) = temp_row(1)
          error(p)    = errorp
      END IF
      press(p) = ptot_row(1)/up
      css(p)   = cs_row(1)/uv
      cvs(p)   = cv_row(1)*(unm/uen)
      dPdT(p)  = dpt_row(1)*(unl**3/uen)
      cps(p)   = cp_row(1)/(uen/unm)
!
      END SUBROUTINE degenerate2
