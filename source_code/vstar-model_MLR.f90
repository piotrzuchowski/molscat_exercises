! VINIT/VSTAR/VSTAR1/VSTAR2 for the 8-parameter model potential of
! PotentialGenerator/MODEL_POTENTIAL_MLR_MORSE.md (atom + linear rotor, ITYPE = 1):
!
!   V(R,theta) = V0(R) + V1(R) P1(cos theta) + V2(R) P2(cos theta)
!
!   V0      : MLR (N = 0, r_ref = r_e, p = 6) from De, C6, rho, a
!   V_lam   : s_lam A0 exp(-b R) - q_lam De exp(-(b/2)(R - r_e)),  lam = 1, 2
!             (A0, b: Born-Mayer fit of V0's wall through V0 = 100 De and V0 = De)
!   R < R_lo: V0 continues; V_lam = V0 * [V_lam/V0](R_lo)
!
! Parameters are read from the namelist file model_MLR.inp in the working directory:
!
!   &MODEL_MLR  ID = 'p000001', DE = 30.0, C6 = 12.0, RHO = 1.5, A = 6.5,
!               S1 = 0.2, Q1 = -0.15, S2 = -0.25, Q2 = 0.2 /
!
!   DE [cm-1], C6 [Eh bohr^6], RHO, A, S1, Q1, S2, Q2 dimensionless.
!   S1, Q1, S2, Q2 default to 0 (e.g. S1 = Q1 = 0 for a homonuclear partner).
!
! VINIT sets RM = 1 bohr (in Angstrom) and EPSIL = 1 cm-1, so VSTAR receives R in bohr and
! returns V in cm-1. The order of each term is taken from LAMBDA(I) (module potential), so
! &POTL may list any subset of LAMBDA = 0, 1, 2 with NTERM = -1 for each term.
! VSTAR1 and VSTAR2 return analytic first and second radial derivatives.

      module model_mlr_data
      implicit none
      double precision, parameter :: HARTREE = 219474.63d0          ! cm-1 per Eh (as in the .md)
      double precision, parameter :: BOHR = 0.529177210903d0        ! Angstrom
      character(len=40) :: id = ' '
      double precision  :: de = -1d0, c6 = -1d0, rho = -1d0, a = -1d0
      double precision  :: s(2) = 0d0, q(2) = 0d0
      ! derived quantities (bohr, cm-1, bohr-1)
      double precision  :: re, binf, b0, dlt, a0, b, rlo, rhi, sigma, clo(2)
      logical           :: loaded = .false.
      contains

      subroutine v0all(rin, v, v1, v2)
!     V0 and its first and second derivatives (cm-1, cm-1/bohr, cm-1/bohr^2) at rin (bohr)
      double precision, intent(in)  :: rin
      double precision, intent(out) :: v, v1, v2
      double precision :: r, r6, re6, d, y, y1, y2, g, g1, g2, u, l1, l2, u1, u2
      r = max(rin, 1d-6)
      r6 = r**6
      re6 = re**6
      d = r6 + re6
      y = (r6 - re6)/d
      y1 = 12d0*re6*r**5/d**2
      y2 = 12d0*re6*r**4*(5d0*d - 12d0*r6)/d**3
!     beta(y)*y = b0*y + (binf - b0)*y**2
      g = b0*y + dlt*y*y
      g1 = (b0 + 2d0*dlt*y)*y1
      g2 = 2d0*dlt*y1*y1 + (b0 + 2d0*dlt*y)*y2
!     u = (re/R)**6 * exp(-beta*y);  V0 = De*(1 - u)**2 - De
      u = (re/r)**6*exp(-g)
      l1 = -6d0/r - g1
      l2 = 6d0/(r*r) - g2
      u1 = u*l1
      u2 = u*(l1*l1 + l2)
      v = de*(1d0 - u)**2 - de
      v1 = -2d0*de*(1d0 - u)*u1
      v2 = 2d0*de*(u1*u1 - (1d0 - u)*u2)
      end subroutine v0all

      subroutine vlamall(k, r, v, v1, v2)
!     V_k (k = 1, 2) and its derivatives at r (bohr)
      integer, intent(in) :: k
      double precision, intent(in)  :: r
      double precision, intent(out) :: v, v1, v2
      double precision :: e1, e2, w, w1, w2
      if (r .ge. rlo) then
        e1 = s(k)*a0*exp(-b*r)
        e2 = q(k)*de*exp(-0.5d0*b*(r - re))
        v = e1 - e2
        v1 = -b*e1 + 0.5d0*b*e2
        v2 = b*b*e1 - 0.25d0*b*b*e2
      else
        call v0all(r, w, w1, w2)
        v = clo(k)*w
        v1 = clo(k)*w1
        v2 = clo(k)*w2
      endif
      end subroutine vlamall

      double precision function v0root(target, rl, rr)
!     R in [rl, rr] with V0(R) = target, by bisection (V0 decreasing on the inner wall)
      double precision, intent(in) :: target, rl, rr
      double precision :: x1, x2, xm, v, v1, v2
      integer :: it
      x1 = rl
      x2 = rr
      call v0all(x1, v, v1, v2)
      if (v .le. target) then
        write(6,'(/2A,1PG14.6,A,0PF8.4,A)') '  *** MODEL_MLR: V0 at the lower bracket is below ', &
       'the target', target, ' cm-1 at R =', x1, ' bohr'
        stop
      endif
      do it = 1, 200
        xm = 0.5d0*(x1 + x2)
        call v0all(xm, v, v1, v2)
        if (v .gt. target) then
          x1 = xm
        else
          x2 = xm
        endif
        if (x2 - x1 .le. 1d-15*rr) exit
      enddo
      v0root = 0.5d0*(x1 + x2)
      end function v0root

      subroutine load
      double precision :: s1, q1, s2, q2, v, v1, v2, w
      integer :: ios, k
      namelist /model_mlr/ id, de, c6, rho, a, s1, q1, s2, q2
      s1 = 0d0
      q1 = 0d0
      s2 = 0d0
      q2 = 0d0
      open(801, file='model_MLR.inp', status='old', action='read', iostat=ios)
      if (ios .ne. 0) then
        write(6,'(/A)') '  *** MODEL_MLR: cannot open model_MLR.inp in the working directory'
        stop
      endif
      read(801, nml=model_mlr, iostat=ios)
      close(801)
      if (ios .ne. 0) then
        write(6,'(/A)') '  *** MODEL_MLR: error reading namelist &MODEL_MLR from model_MLR.inp'
        stop
      endif
      if (de .le. 0d0 .or. c6 .le. 0d0 .or. rho .le. 0d0 .or. a .le. 0d0) then
        write(6,'(/A)') '  *** MODEL_MLR: DE, C6, RHO and A must all be given and positive'
        stop
      endif
      s = (/ s1, s2 /)
      q = (/ q1, q2 /)

      re = (c6*HARTREE/(rho*de))**(1d0/6d0)
      binf = log(2d0/rho)
      b0 = 2d0*(a - 6d0)/6d0
      dlt = binf - b0
      rlo = v0root(100d0*de, 0.2d0*re, re)
      rhi = v0root(de, rlo, re)
      sigma = v0root(0d0, rhi, re)
      b = log(100d0)/(rhi - rlo)
      a0 = de*exp(b*rhi)
      call v0all(rlo, w, v1, v2)
      do k = 1, 2
        call vlamall(k, rlo, v, v1, v2)
        clo(k) = v/w
      enddo
      loaded = .true.

      write(6,'(/A/A,A)') '  MODEL_MLR potential (PotentialGenerator/MODEL_POTENTIAL_MLR_MORSE.md)', &
       '    ID   = ', trim(id)
      write(6,'(A,1PG16.8,A)') '    De   =', de,   ' cm-1', &
       '    C6   =', c6,   ' Eh bohr^6', &
       '    rho  =', rho,  ' ', &
       '    a    =', a,    ' '
      write(6,'(A,4F10.5)') '    s1, q1, s2, q2 =', s(1), q(1), s(2), q(2)
      write(6,'(A)') '  derived quantities:'
      write(6,'(A,1PG16.8,A)') '    r_e     =', re,    ' bohr', &
       '    sigma   =', sigma, ' bohr  (V0 = 0)', &
       '    R_lo    =', rlo,   ' bohr  (V0 = 100 De)', &
       '    R_hi    =', rhi,   ' bohr  (V0 = De)', &
       '    b       =', b,     ' bohr-1', &
       '    A0      =', a0,    ' cm-1', &
       '    k       =', 2d0*de*a*a/(re*re), ' cm-1 bohr-2  (V0 curvature at r_e)', &
       '    beta_inf=', binf,  ' ', &
       '    beta_0  =', b0,    ' '
      write(6,'(A)') '  for R < R_lo the ratios V1/V0 and V2/V0 are frozen at their R_lo values'
      end subroutine load
      end module model_mlr_data


      SUBROUTINE VINIT(I, RM, EPSIL)
      USE model_mlr_data, ONLY: loaded, load, BOHR
      USE potential, ONLY: LAMBDA, RMNAME, EPNAME
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: I
      DOUBLE PRECISION, INTENT(INOUT) :: RM, EPSIL
      IF (.NOT. loaded) CALL load
      IF (LAMBDA(I).LT.0 .OR. LAMBDA(I).GT.2) THEN
        WRITE(6,'(/A,I4,A,I4)') '  *** MODEL_MLR: only LAMBDA = 0, 1, 2 are defined; term', &
       I, ' has LAMBDA =', LAMBDA(I)
        STOP
      ENDIF
      RM = BOHR
      EPSIL = 1.D0
      RMNAME = 'BOHR'
      EPNAME = 'CM-1'
      END SUBROUTINE VINIT

      SUBROUTINE VSTAR(I, R, V)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: I
      DOUBLE PRECISION, INTENT(IN) :: R
      DOUBLE PRECISION, INTENT(OUT) :: V
      CALL MLRDER(I, R, V, 0)
      END SUBROUTINE VSTAR

      SUBROUTINE VSTAR1(I, R, V)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: I
      DOUBLE PRECISION, INTENT(IN) :: R
      DOUBLE PRECISION, INTENT(OUT) :: V
      CALL MLRDER(I, R, V, 1)
      END SUBROUTINE VSTAR1

      SUBROUTINE VSTAR2(I, R, V)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: I
      DOUBLE PRECISION, INTENT(IN) :: R
      DOUBLE PRECISION, INTENT(OUT) :: V
      CALL MLRDER(I, R, V, 2)
      END SUBROUTINE VSTAR2

      SUBROUTINE MLRDER(I, R, V, IDER)
!     IDER-th radial derivative of the term with LAMBDA(I) at R (bohr), in cm-1/bohr**IDER
      USE model_mlr_data, ONLY: v0all, vlamall
      USE potential, ONLY: LAMBDA
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: I, IDER
      DOUBLE PRECISION, INTENT(IN) :: R
      DOUBLE PRECISION, INTENT(OUT) :: V
      DOUBLE PRECISION :: D(0:2)
      IF (LAMBDA(I).EQ.0) THEN
        CALL v0all(R, D(0), D(1), D(2))
      ELSE
        CALL vlamall(LAMBDA(I), R, D(0), D(1), D(2))
      ENDIF
      V = D(IDER)
      END SUBROUTINE MLRDER
