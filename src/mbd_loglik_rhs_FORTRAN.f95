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

      ! length of the vector - decided in R-code
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
! Initialisation: name of this function as passed by "initfunc" argument
! Sets the fixed parameter vector, and allocates memory
!==========================================================================
!==========================================================================

      SUBROUTINE mbd_initmodpc (steadyparms)
      USE mbd_dimmod

      IMPLICIT NONE
      EXTERNAL steadyparms

      INTEGER, PARAMETER :: nparsmall = 1  ! constant-length parameters

      DOUBLE PRECISION parms(nparsmall)
      COMMON /XCBPar/parms                 ! common block

! Set the fixed parameters obtained from R
      CALL steadyparms(nparsmall, parms)

! first parameter has the length of the vector
      N = INT(sqrt(parms(1)) + 1e-1)

! Allocate variable size arrays (state variables, derivatives and parameters)

      IF (ALLOCATED(P)) DEALLOCATE(P)
      ALLOCATE(P(3 * N ** 2))

      initialised = .FALSE.

      END SUBROUTINE mbd_initmodpc



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
      DOUBLE PRECISION  :: t, Conc(N), dConc(N), yout(*)
      !REAL(16)          :: V(N)

! parameters - named here
      DOUBLE PRECISION rn
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

      !DO I = 1, N
      !  V(I) = 0
      !  DO II = 1, N
      !    V(I) = V(I) + P((II - 1) * N + I) * Conc(II)
      !  ENDDO
      !  dConc(I) = V(I)
      !ENDDO

      dConc = MATMUL(RESHAPE(P,(/N,N/), order = (/1,2/)),Conc)

      END SUBROUTINE mbd_runmod


!==========================================================================
!==========================================================================
! Dynamic routine: name of this function as passed by "func" argument
! variable parameter values are passed via yout
!==========================================================================
!==========================================================================

      SUBROUTINE mbd_runmodpcp (neq, t, Conc, dConc, yout, ip)
      USE mbd_dimmod
      IMPLICIT NONE

!......................... declaration section.............................
      INTEGER           :: neq, ip(*), i, ii, lx
      DOUBLE PRECISION  :: t, Conc(N ** 2), dConc(N ** 2), yout(*)
      DOUBLE PRECISION  :: vec(N)
      DOUBLE PRECISION  :: dp(N,N), V(N, N), V2(N + 2, N + 2)
      DOUBLE PRECISION  :: nu_q_mat(N,N), m1(N, N), m2(N,N)

! parameters - named here
      DOUBLE PRECISION rn
      COMMON /XCBPar/rn

! local variables
      CHARACTER(len=100) msg

!............................ statements ..................................

      IF (.NOT. Initialised) THEN
        ! check memory allocated to output variables
        IF (ip(1) < 1) CALL rexit("nout not large enough")

        ! save parameter values in yout
        ii = ip(1)   ! Start of parameter values
        CALL mbd_fill1d(P, 3 + N ** 2, yout, ii)   ! ii is updated in fill1d
        Initialised = .TRUE.          ! to prevent from initialising more than once
      ENDIF

! dynamics

   V = RESHAPE(Conc,(/N,N/), order = (/1,2/))
   V2 = 0
   V2(2:(N+1),2:(N+1)) = V
   nu_q_mat = RESHAPE(P((3 + 1):(3 + N ** 2)),(/N,N/), order = (/1,2/))
   !m1 = RESHAPE(P((3 + N ** 2 + 1):(3 + 2 * N ** 2)),(/N,N/), order = (/1,2/))
   vec = (/(I, I = 0, N - 1, 1)/)
   DO I = 1, N
     m1(I,:) = vec
   ENDDO
   m2 = TRANSPOSE(m1)

   dp=P(1)*((m1-1)*V2(2:(N+1),1:N)+(m2-1)*V2(1:N,2:(N+1))-(m1+m2)*V)
   dp=dp+P(2)*((m1+1)*V2(2:(N+1),3:(N+2))+(m2+1)*V2(3:(N+2),2:(N+1))-(m1+m2)*V)
   dp=dp+P(3)*(MATMUL(MATMUL(nu_q_mat,V),TRANSPOSE(nu_q_mat)) - V)

   dConc = RESHAPE(dp,(/N ** 2/))

   END SUBROUTINE mbd_runmodpcp


!==========================================================================
!==========================================================================
! Dynamic routine: name of this function as passed by "func" argument
! variable parameter values are passed via yout
!==========================================================================
!==========================================================================

      SUBROUTINE mbd_runmodpcq (neq, t, Conc, dConc, yout, ip)
      USE mbd_dimmod
      IMPLICIT NONE

!......................... declaration section.............................
      INTEGER           :: neq, ip(*), i, ii, lx, I1, J1, n1
      DOUBLE PRECISION  :: t, Conc(N ** 2), dConc(N ** 2), yout(*)
      DOUBLE PRECISION  :: vec(N)
      DOUBLE PRECISION  :: dq(N,N), V(N, N), V2(N + 2, N + 2)
      DOUBLE PRECISION  :: nu_q_mat(N,N), m1(N, N), m2(N,N)

! parameters - named here
      DOUBLE PRECISION rn
      COMMON /XCBPar/rn

! local variables
      CHARACTER(len=100) msg

!............................ statements ..................................

      IF (.NOT. Initialised) THEN
        ! check memory allocated to output variables
        IF (ip(1) < 1) CALL rexit("nout not large enough")

        ! save parameter values in yout
        ii = ip(1)   ! Start of parameter values
        CALL mbd_fill1d(P, 3 + N ** 2, yout, ii)   ! ii is updated in fill1d
        Initialised = .TRUE.          ! to prevent from initialising more than once
      ENDIF

! dynamics

   V = RESHAPE(Conc,(/N,N/), order = (/1,2/))
   V2 = 0
   V2(2:(N+1),2:(N+1)) = V
   nu_q_mat = RESHAPE(P((3 + 1):(3 + N ** 2)),(/N,N/), order = (/1,2/))
   !m1 = RESHAPE(P((3 + N ** 2 + 1):(3 + 2 * N ** 2)),(/N,N/), order = (/1,2/))
   vec = (/(I, I = 0, N - 1, 1)/)
   DO I = 1, N
     m1(I,:) = vec
   ENDDO
   m2 = TRANSPOSE(m1)

   dq=P(1)*((m1+1)*V2(2:(N+1),1:N)+(m2+1)*V2(1:N,2:(N+1))-(m1+m2+2)*V)
   dq=dq+P(2)*((m1+1)*V2(2:(N+1),3:(N+2))+(m2+1)*V2(3:(N+2),2:(N+1))-(m1+m2+2)*V)
   dq=dq+P(3)*(MATMUL(MATMUL(nu_q_mat,V),TRANSPOSE(nu_q_mat)) - V)

   dConc = RESHAPE(dq,(/N ** 2/))


!    DO I = 0, N - 1
!      DO II = 0, N - 1
!        I1 = I + 1
!        J1 = II + 1
!        dp1=J1*Conc2(I1+1,J1)+I1*Conc2(I1,J1+1)-(I1+J1)*Conc2(I1+1,J1+1)
!        dp2=J1*Conc2(I1+1,J1+2)+I1*Conc2(I1+2,J1+1)-(I1+J1)*Conc2(I1+1,J1+1)
!        dp3 = -Conc((I1 - 1) * N + J1)
!        DO n1 = 1, N
!          dp3 = dp3 + aux1(I1,n1) * P(3 + (n1 - 1) * N + J1)
!        ENDDO
!        dConc((I1 - 1)*N + J1) = P(1)*dp1 + P(2)*dp2 + P(3)*dp3
!      ENDDO
!    ENDDO

   END SUBROUTINE mbd_runmodpcq
