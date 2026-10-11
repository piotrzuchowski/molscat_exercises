! VINIT/VSTAR/VSTAR1/VSTAR2 for the ionic model potential of
! PotentialGenerator/MODEL_POTENTIAL_MLR_MORSE_ION.md (molecular ion + He or para-H2, ITYPE = 1):
!
!   V(R,theta) = V0(R) + V1(R) P1(cos theta) + V2(R) P2(cos theta)
!
!   V0      : MLR (N = 1, p = 6, q = 3, r_ref = r_e) with u(R) = C4/R^4 + C6/R^6;
!             beta_1 is set by the wall position sigma = sigma4p (1 + delta)
!   V1      : s1 A0 exp(-b R) - q1 De exp(-(b/2)(R - r_e)) - f5(b R) C5/R^5
!   V2      : s2 A0 exp(-b R) - q2 De exp(-(b/2)(R - r_e))
!             (A0, b: Born-Mayer fit of V0's wall through V0 = 100 De and V0 = De;
!              f5: Tang-Toennies damping; C5 = g1 De r_e^5)
!   R < R_lo: V0 continues; V_lam = V0 * [V_lam/V0](R_lo)
!
! Parameters are read from the namelist file model_MLR_ion.inp in the working directory:
!
!   &MODEL_MLR_ION  ID = 'i000001', DE = 300.0, C4 = 0.69, K6 = 0.3, RHO = 1.8, A = 4.2,
!                   DELTA = 0.01, S1 = 0.1, Q1 = -0.1, G1 = 0.7, S2 = -0.15, Q2 = 0.12 /
!
!   DE [cm-1], C4 [Eh bohr^4]; K6, RHO, A, DELTA, S1, Q1, G1, S2, Q2 dimensionless.
!   DELTA, S1, Q1, G1, S2, Q2 default to 0 (e.g. S1 = Q1 = G1 = 0 for a homonuclear ion).
!
! VINIT sets RM = 1 bohr (in Angstrom) and EPSIL = 1 cm-1, so VSTAR receives R in bohr and
! returns V in cm-1. The order of each term is taken from LAMBDA(I) (module potential), so
! &POTL may list any subset of LAMBDA = 0, 1, 2 with NTERM = -1 for each term.
! VSTAR1 and VSTAR2 return analytic first and second radial derivatives.

      module model_mlr_ion_data
      implicit none
      double precision, parameter :: HARTREE = 219474.63d0          ! cm-1 per Eh (as in the .md)
      double precision, parameter :: BOHR = 0.529177210903d0        ! Angstrom
      integer, parameter :: NP = 6, NQ = 3
      character(len=40) :: id = ' '
      double precision  :: de = -1d0, c4 = -1d0, k6 = -1d0, rho = -1d0, a = -1d0, delta = 0d0
      double precision  :: s(2) = 0d0, q(2) = 0d0, g1 = 0d0
      ! derived quantities (bohr, cm-1, bohr-1, au)
      double precision  :: re, c6, ue, meff, binf, b0, beta1, sig4p, sigma
      double precision  :: a0, b, rlo, rhi, c5, clo(2)
      logical           :: loaded = .false.
      contains

      subroutine yfun(n, r, y, y1, y2)
!     y_n(R) = (R^n - re^n)/(R^n + re^n) and its first two derivatives
      integer, intent(in) :: n
      double precision, intent(in)  :: r
      double precision, intent(out) :: y, y1, y2
      double precision :: rn, ren, d
      rn = r**n
      ren = re**n
      d = rn + ren
      y = (rn - ren)/d
      y1 = 2d0*n*ren*r**(n-1)/d**2
      y2 = 2d0*n*ren*r**(n-2)*((n-1)*d - 2d0*n*rn)/d**3
      end subroutine yfun

      subroutine v0all(rin, bet1, v, v1, v2)
!     V0 and its first and second derivatives (cm-1, cm-1/bohr, cm-1/bohr^2) at rin (bohr)
      double precision, intent(in)  :: rin, bet1
      double precision, intent(out) :: v, v1, v2
      double precision :: r, yp, yp1, yp2, yq, yq1, yq2, bb, bb1, bb2, bt, bt1, bt2
      double precision :: g, g1d, g2d, u, u1, u2, e, w, w1, w2
      r = max(rin, 1d-6)
      call yfun(NP, r, yp, yp1, yp2)
      call yfun(NQ, r, yq, yq1, yq2)
!     beta(R) = binf*yp + (1 - yp)*B,  B = b0 + beta1*yq;  exponent G = beta*yp
      bb = b0 + bet1*yq
      bb1 = bet1*yq1
      bb2 = bet1*yq2
      bt = binf*yp + (1d0 - yp)*bb
      bt1 = binf*yp1 - yp1*bb + (1d0 - yp)*bb1
      bt2 = binf*yp2 - yp2*bb - 2d0*yp1*bb1 + (1d0 - yp)*bb2
      g = bt*yp
      g1d = bt1*yp + bt*yp1
      g2d = bt2*yp + 2d0*bt1*yp1 + bt*yp2
!     u(R) = H (C4/R^4 + C6/R^6);  w = u/u(re) * exp(-G);  V0 = De (1 - w)^2 - De
      u = HARTREE*(c4/r**4 + c6/r**6)
      u1 = HARTREE*(-4d0*c4/r**5 - 6d0*c6/r**7)
      u2 = HARTREE*(20d0*c4/r**6 + 42d0*c6/r**8)
      e = exp(-g)/ue
      w = e*u
      w1 = e*(u1 - u*g1d)
      w2 = e*(u2 - 2d0*u1*g1d + u*g1d*g1d - u*g2d)
      v = de*(1d0 - w)**2 - de
      v1 = -2d0*de*(1d0 - w)*w1
      v2 = 2d0*de*(w1*w1 - (1d0 - w)*w2)
      end subroutine v0all

      subroutine c5all(r, v, v1, v2)
!     -f5(bR) H C5/R^5 and its derivatives; f5(x) = 1 - exp(-x) sum_{k=0..5} x^k/k!,
!     f5'(x) = exp(-x) x^5/5!,  f5''(x) = exp(-x) (5 x^4 - x^5)/5!
      double precision, intent(in)  :: r
      double precision, intent(out) :: v, v1, v2
      double precision :: x, ex, sm, term, f, f1, f2, c
      integer :: k
      x = b*r
      ex = exp(-x)
      sm = 1d0
      term = 1d0
      do k = 1, 5
        term = term*x/k
        sm = sm + term
      enddo
      f = 1d0 - ex*sm
      f1 = ex*x**5/120d0
      f2 = ex*(5d0*x**4 - x**5)/120d0
      c = HARTREE*c5
      v = -c*f/r**5
      v1 = -c*(b*f1/r**5 - 5d0*f/r**6)
      v2 = -c*(b*b*f2/r**5 - 10d0*b*f1/r**6 + 30d0*f/r**7)
      end subroutine c5all

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
        if (k .eq. 1 .and. c5 .ne. 0d0) then
          call c5all(r, w, w1, w2)
          v = v + w
          v1 = v1 + w1
          v2 = v2 + w2
        endif
      else
        call v0all(r, beta1, w, w1, w2)
        v = clo(k)*w
        v1 = clo(k)*w1
        v2 = clo(k)*w2
      endif
      end subroutine vlamall

      double precision function bisect(f, x1in, x2in, what)
!     root of f in [x1in, x2in] by bisection; f must change sign on the interval
      interface
        double precision function f(x)
        double precision, intent(in) :: x
        end function f
      end interface
      double precision, intent(in) :: x1in, x2in
      character(len=*), intent(in) :: what
      double precision :: x1, x2, xm, f1, fm
      integer :: it
      x1 = x1in
      x2 = x2in
      f1 = f(x1)
      if (f1*f(x2) .gt. 0d0) then
        write(6,'(/3A,2(1PG14.6))') '  *** MODEL_MLR_ION: no sign change for ', what, &
          ' on the interval', x1, x2
        stop
      endif
      do it = 1, 300
        xm = 0.5d0*(x1 + x2)
        fm = f(xm)
        if (sign(1d0, fm) .eq. sign(1d0, f1)) then
          x1 = xm
          f1 = fm
        else
          x2 = xm
        endif
        if (abs(x2 - x1) .le. 1d-15*max(abs(x1), abs(x2), 1d0)) exit
      enddo
      bisect = 0.5d0*(x1 + x2)
      end function bisect

      subroutine load
      double precision :: s1, q1, s2, q2, v, v1, v2, w
      double precision :: target, sigtry
      integer :: ios, k
      namelist /model_mlr_ion/ id, de, c4, k6, rho, a, delta, s1, q1, g1, s2, q2
      s1 = 0d0
      q1 = 0d0
      s2 = 0d0
      q2 = 0d0
      open(801, file='model_MLR_ion.inp', status='old', action='read', iostat=ios)
      if (ios .ne. 0) then
        write(6,'(/A)') '  *** MODEL_MLR_ION: cannot open model_MLR_ion.inp in the working directory'
        stop
      endif
      read(801, nml=model_mlr_ion, iostat=ios)
      close(801)
      if (ios .ne. 0) then
        write(6,'(/A)') '  *** MODEL_MLR_ION: error reading namelist &MODEL_MLR_ION from model_MLR_ion.inp'
        stop
      endif
      if (de .le. 0d0 .or. c4 .le. 0d0 .or. k6 .lt. 0d0 .or. rho .le. 0d0 .or. a .le. 0d0) then
        write(6,'(/A)') '  *** MODEL_MLR_ION: DE, C4, RHO, A (> 0) and K6 (>= 0) must all be given'
        stop
      endif
      s = (/ s1, s2 /)
      q = (/ q1, q2 /)

      re = (c4*(1d0 + k6)*HARTREE/(rho*de))**0.25d0
      c6 = k6*c4*re**2
      ue = HARTREE*(c4/re**4 + c6/re**6)
      meff = (4d0 + 6d0*k6)/(1d0 + k6)
      binf = log(2d0/rho)
      b0 = 2d0*(a - meff)/NP
!     wall position: zero of V0 with beta1 = 0, scaled by (1 + delta), then beta1 to match it
      sig4p = bisect(v0zero_b0, 0.3d0*re, re, 'sigma4p (V0 = 0, beta1 = 0)')
      sigma = sig4p*(1d0 + delta)
      sigtry = sigma
      beta1 = bisect(v0_at_sigma, -30d0, 30d0, 'beta1 (V0(sigma) = 0)')
!     Born-Mayer fit of the wall
      target = 100d0*de
      rlo = bisect(v0_minus_target, 0.05d0*re, re, 'R_lo (V0 = 100 De)')
      target = de
      rhi = bisect(v0_minus_target, rlo, re, 'R_hi (V0 = De)')
      b = log(100d0)/(rhi - rlo)
      a0 = de*exp(b*rhi)
      c5 = g1*de*re**5/HARTREE
      call v0all(rlo, beta1, w, v1, v2)
      do k = 1, 2
        call vlamall(k, rlo, v, v1, v2)
        clo(k) = v/w
      enddo
      loaded = .true.

      write(6,'(/A/A,A)') '  MODEL_MLR_ION potential (PotentialGenerator/MODEL_POTENTIAL_MLR_MORSE_ION.md)', &
        '    ID   = ', trim(id)
      write(6,'(A,1PG16.8,A)') '    De    =', de,    ' cm-1', &
        '    C4    =', c4,    ' Eh bohr^4', &
        '    k6    =', k6,    ' ', &
        '    rho   =', rho,   ' ', &
        '    a     =', a,     ' ', &
        '    delta =', delta, ' '
      write(6,'(A,5F10.5)') '    s1, q1, g1, s2, q2 =', s(1), q(1), g1, s(2), q(2)
      write(6,'(A)') '  derived quantities:'
      write(6,'(A,1PG16.8,A)') '    r_e     =', re,    ' bohr', &
        '    C6      =', c6,    ' Eh bohr^6', &
        '    C5      =', c5,    ' Eh bohr^5', &
        '    sigma   =', sigma, ' bohr  (V0 = 0)', &
        '    sigma4p =', sig4p, ' bohr  (V0 = 0 with beta1 = 0)', &
        '    beta_1  =', beta1, ' ', &
        '    R_lo    =', rlo,   ' bohr  (V0 = 100 De)', &
        '    R_hi    =', rhi,   ' bohr  (V0 = De)', &
        '    b       =', b,     ' bohr-1', &
        '    A0      =', a0,    ' cm-1', &
        '    k       =', 2d0*de*a*a/(re*re), ' cm-1 bohr-2  (V0 curvature at r_e)', &
        '    m_eff   =', meff,  ' ', &
        '    beta_inf=', binf,  ' ', &
        '    beta_0  =', b0,    ' '
      write(6,'(A)') '  for R < R_lo the ratios V1/V0 and V2/V0 are frozen at their R_lo values'

      contains
        double precision function v0zero_b0(x)
        double precision, intent(in) :: x
        double precision :: d1, d2
        call v0all(x, 0d0, v0zero_b0, d1, d2)
        end function v0zero_b0

        double precision function v0_at_sigma(z)
        double precision, intent(in) :: z
        double precision :: d1, d2
        call v0all(sigtry, z, v0_at_sigma, d1, d2)
        end function v0_at_sigma

        double precision function v0_minus_target(x)
        double precision, intent(in) :: x
        double precision :: d1, d2
        call v0all(x, beta1, v0_minus_target, d1, d2)
        v0_minus_target = v0_minus_target - target
        end function v0_minus_target
      end subroutine load
      end module model_mlr_ion_data


      SUBROUTINE VINIT(I, RM, EPSIL)
      USE model_mlr_ion_data, ONLY: loaded, load, BOHR
      USE potential, ONLY: LAMBDA, RMNAME, EPNAME
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: I
      DOUBLE PRECISION, INTENT(INOUT) :: RM, EPSIL
      IF (.NOT. loaded) CALL load
      IF (LAMBDA(I).LT.0 .OR. LAMBDA(I).GT.2) THEN
        WRITE(6,'(/A,I4,A,I4)') '  *** MODEL_MLR_ION: only LAMBDA = 0, 1, 2 are defined; term', &
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
      USE model_mlr_ion_data, ONLY: v0all, vlamall, beta1
      USE potential, ONLY: LAMBDA
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: I, IDER
      DOUBLE PRECISION, INTENT(IN) :: R
      DOUBLE PRECISION, INTENT(OUT) :: V
      DOUBLE PRECISION :: D(0:2)
      IF (LAMBDA(I).EQ.0) THEN
        CALL v0all(R, beta1, D(0), D(1), D(2))
      ELSE
        CALL vlamall(LAMBDA(I), R, D(0), D(1), D(2))
      ENDIF
      V = D(IDER)
      END SUBROUTINE MLRDER
