! Helper function:
! fill vec with N elements from parms, starting at position ii
!==========================================================================

      SUBROUTINE mbd_fill1d (vec, DIMP, parms, II)
      IMPLICIT NONE
      INTEGER DIMP, II, I
      DOUBLE PRECISION vec(DIMP), parms(*)
      II = II
        DO I = 1, DIMP
          II = II + 1
          vec(I) = parms(II)
        ENDDO

      END SUBROUTINE mbd_fill1d

!==========================================================================
! module with declarations
!==========================================================================

      MODULE mbd_dimmod

      ! length of the vector -  decided in R-code
      INTEGER  :: N

      ! 1 parameter vectors with unknown length
      DOUBLE PRECISION, ALLOCATABLE  :: P(:)

      ! Boolean: will become TRUE if the parameters have a value
      LOGICAL :: initialised = .FALSE.

      END MODULE mbd_dimmod

!==========================================================================
!==========================================================================
! Initialisation: name of this function as passed by "initfunc" argument
! Sets the fixed parameter vector, and allocates memory
!==========================================================================
!==========================================================================

      SUBROUTINE mbd_initmod (steadyparms)
      USE mbd_dimmod

      IMPLICIT NONE
      EXTERNAL steadyparms

      INTEGER, PARAMETER :: nparsmall = 1  ! constant-length parameters

      DOUBLE PRECISION parms(nparsmall)
      COMMON /XCBPar/parms                 ! common block

! Set the fixed parameters obtained from R
      CALL steadyparms(nparsmall, parms)

! first parameter has the length of the vector
      N = INT(parms(1) + 1e-6)

! Allocate variable size arrays (state variables, derivatives and parameters)

      IF (ALLOCATED(P)) DEALLOCATE(P)
      ALLOCATE(P(N ** 2))

      initialised = .FALSE.

      END SUBROUTINE mbd_initmod

!==========================================================================
!==========================================================================
! Dynamic routine: name of this function as passed by "func" argument
! variable parameter values are passed via yout
!==========================================================================
!==========================================================================

      SUBROUTINE mbd_runmod (neq, t, Conc, dConc, yout, ip)
      USE mbd_dimmod
      IMPLICIT NONE

!......................... declaration section.............................
      INTEGER           :: neq, ip(*), i, ii
      DOUBLE PRECISION  :: t, Conc(N), yout(*)
      REAL(16)          :: dConc(N)

! parameters - named here
      DOUBLE PRECISION rn(1)
      COMMON /XCBPar/rn

! local variables
      CHARACTER(len=100) msg

!............................ statements ..................................

      IF (.NOT. Initialised) THEN
        ! check memory allocated to output variables
        IF (ip(1) < 1) CALL rexit("nout not large enough")

        ! save parameter values in yout
        ii = ip(1)   ! Start of parameter values
        CALL mbd_fill1d(P, N ** 2, yout, ii)   ! ii is updated in fill1d
        Initialised = .TRUE.          ! to prevent from initialising more than once
      ENDIF

! dynamics

!mbd_loglik_rhs <- function(t, x, params) {
!  list(params %*% x)
!}

      DO I = 1, N
        dConc(I) = 0
        DO II = 1, N
          dConc(I) = dConc(I) + P((I - 1) * N + II) * Conc(II)
        ENDDO
      ENDDO

      END SUBROUTINE mbd_runmod

