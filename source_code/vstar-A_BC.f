      SUBROUTINE  VINIT(II,RM,RESULT)
C -------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NANGMAX=50,NRADMAX=250)
      SAVE
      DOUBLE PRECISION, ALLOCATABLE:: ANGLS(:),ALPHAS(:,:),BETA(:)
      DOUBLE PRECISION, ALLOCATABLE:: POL(:,:), POLN(:,:),WT(:)
      DOUBLE PRECISION, ALLOCATABLE:: PTST(:)
      DOUBLE PRECISION, ALLOCATABLE:: RPT(:,:), VPT(:,:), RR(:),P(:)
      INTEGER, ALLOCATABLE:: NRA(:)
      INTEGER LAM(2)
      DATA LAM /0,2/
      COMMON /PES_DATA/ LMAX

C
C  MGNH - POTENTIAL OBTAINED IN AVQZ (AVTZ ON CR) BASIS
C  SET WITH  MIDBOND

      write(*,*) ' ii ',II  
      IF (II.NE.1) GOTO 15

C  READ THE FILE POT.DATA
C  THE STRUCTURE OF FILE IS FOLLOWING:
C  NANGLE  ICASEM: NUMBER OF ANGULAR CUTS,
C          ICASE = IF THE MORPHING IS SWITCHED ON
C  ANGLE1
C      XX1  EINT1 EINT2  EEX  EIND  EDISP
C      ...
C      XXNR1   EINT1 EINT2  EEX  EIND  EDISP
C  ANGLE2
C  ETC...
      OPEN(9,FILE='pot-A+BC.data',
     O     FORM='FORMATTED',STATUS='OLD',ERR=2010)
      READ(9,*) NANGLES, AA,BB,  C60,C80,C62,C82
      write(*,*) ' nangles  ',nangles  
      ALLOCATE(ANGLS(NANGLES))
      ALLOCATE(NRA(NANGLES))
      C60=C60*219474.63
      C62=C62*219474.63
      C80=C80*219474.63
      C82=C82*219474.63
C  DETERMINE MAXIMUM NRA(I) TO ALLOCATE PROPERLY ALPHA
C  AND VPT ARRAYS
      NRAMAX = 0
      DO I=1,NANGLES
        READ(9,*) ANG_I,NRA_I
        ANGLS(I) = ANG_I
        NRA(I) = NRA_I
        DO J=1,NRA(I)
          READ(9,*)
          NRAMAX = MAX(NRAMAX,J)
        ENDDO
      ENDDO
      REWIND(9)
      READ(9,*)
      ALLOCATE(VPT(NRAMAX,NANGLES))
      ALLOCATE(ALPHAS(NRAMAX,NANGLES))
      ALLOCATE(RPT(NRAMAX,NANGLES))
      ALLOCATE(P(NANGLES))
      ALLOCATE(RR(NANGLES))

      DO I=1,NANGLES
      READ(9,*) ANGLS(I), NRA(I)
        DO J=1,NRA(I)
           READ(9,*)  RPT0, EINT
           RPT(J,I) = RPT0
           VPT(J,I) = EINT
        ENDDO
      ENDDO
      CLOSE(9)

      WRITE(6,*) ' NUMBER OF ANGLES: ',NANGLES
c      WRITE(6,'(20F8.2)') (ANGLS(I),I=1,NANGLES)
      WRITE(6,*) ' NUMBER OF CALCULATED POINTS FOR EACH ANGLE: '
c      WRITE(6,'(20I8)') (NRA(I),I=1,NANGLES)
      LMAX = NANGLES-1

C  INITIALIZE ARRAY WITH PN FOR EACH ANGLE
C  CALCULATE ASSOCIATED LEGENDRE POLYNOMIALS
C  IT WILL BE NEEDED LATER
C
      ALLOCATE(POL(NANGLES,NANGLES))
      ALLOCATE(POLN(NANGLES,NANGLES)) ! NOT USED HERE
      PI = ACOS(-1.0D0)
      LAMMAX = LMAX
      DO IA=1,NANGLES
        COSTHETA = COS(ANGLS(IA)*PI/180.0D0)
        CALL CALCPOL(NANGLES, COSTHETA,POL(1,IA),POLN(1,IA))
      ENDDO


      NP = 3
      MP = 5
      LP = 1
      NFIX = 0
      ASYM = 0.0D0
      BETA1= 3.D0/56.D0
      BETA2=-4.D0*BETA1/3.D0
      BETA3= 7.D0*BETA1/15.D0
      ALLOCATE(BETA(NP))
      BETA(1) =  BETA1
      BETA(2) =  BETA2
      BETA(3) =  BETA3


C  CALCULATE RKHS EXPANSION COEFFICIENTS 'ALPHA'
C  NOTE THAT FOR THIS SURFACE WE USE N=3, M=5
C  THE ASYMPTOTICS IS NOT FORCED
C  CALCULATE RKHS EXPANSION COEFFICIENTS 'ALPHA'
C  NOTE THAT FOR THIS SURFACE WE USE N=3, M=5

      DO I=1,NANGLES
       WRITE(*,*) (RPT(JJ,I),JJ=1,NRA(I))
        CALL RK_INIT(NRA(I),RPT(1,I),VPT(1,I),ALPHAS(1,I),BETA,
     O               NP, MP, LP, NFIX, ASYM )
      ENDDO

       IP=0
       MXLAM=LMAX
       ALLOCATE(WT(NANGLES))
       DO 70 IA=1,NANGLES
          WT(IA)=2.D0/(DBLE(NANGLES*(NANGLES-1))
     1           *POL(NANGLES,IA)*POL(NANGLES,IA))
   70  CONTINUE


  15  RETURN

C -------------------------------------------------------
      ENTRY VSTAR(II,X,TOTAL)
C -------------------------------------------------------
      write(*,*) '...... ',II
c     IF (LAM(II).GT.LMAX) GOTO 17
      
      TOTALS = 0.0D0
      TOTALL = 0.0D0
      TOTAL  = 0.0D0
C----------------------------------
      if (X.LE.BB) THEN 
      CALL GET_RADIAL(X, RR, VPT, NANGLES, NRA, RPT, WT,
     O                ALPHAS, BETA, NP, MP, LP, NRAMAX )
C
C  RR - OUTPUT - RR FOR EACH ANGLE AT GIVEN X
C

C
C  INTEGRATION - LOOP OVER QUADRATURES WT TIMES VALUE
Ctest        Write(*,*) 'LAM ',LAM(II),(POL(LAM(II)+1,IA),IA=1,NANGLES) 
        DO IA=1, NANGLES
          TOTALS = TOTALS + POL(LAM(II)+1,IA)*RR(IA)*WT(IA)
        ENDDO
        TOTALS = TOTALS * (2.0D0*(LAM(II)) + 1.0D0 )/2.0D0
      ELSE
        TOTALS=0.0d0
      ENDIF
c      write(95+LAM(II),'(I5,F7.2,F21.12)')
c     o LAM(II),X,TOTALS
      if (X.GT.AA) THEN 
         if (LAM(II).EQ.0) THEN
            TOTALL = - C60*X**(-6) -  C80*X**(-8)  
         ELSEIF (LAM(II).EQ.2) THEN
            TOTALL = - C62*X**(-6) -  C82*X**(-8)  
!            TOTALL = TOTALL*5.0/2.0 
         ELSE 
           TOTALL = 0d0
         ENDIF 
      ELSE
        TOTALL = 0.0d0
      ENDIF
      
C----------------------------------
      IF (X.LE.AA) W=0.0d0      
      IF (X.GT.BB) W=1.0d0      
      if (X.GT.AA.AND.X.LE.BB) THEN
       w=fswitch(X,AA,BB)
      ENDIF
      TOTAL = TOTALS*(1.0-w) + TOTALL*w   
      write(89+LAM(II),'(I5,F7.2,4F21.12)')
     o LAM(II),X,TOTALS,TOTALL,TOTAL,w

   16 RETURN
   17 STOP 'WRONG LAMBDA IN VSTAR ROUTINE '
      RETURN
C  --------------------------------------------------------
C  DUMMY ROUTINES
C  --------------------------------------------------------
      ENTRY VSTAR1(II,X,TOTAL)
C --------------------------------------------------------
C --------------------------------------------------------
      ENTRY VSTAR2(II,X,TOTAL)
C --------------------------------------------------------
      WRITE(6,*) 'VSTAR: DERIVATIVES NOT IMPLEMENTED'
      RETURN
 2010 WRITE(6,*) ' ERROR OPENING THE FILE POT.DATA !!! '
      END
C ----------------------------------------------------
      SUBROUTINE RK_INIT(NR,RPT,VPT,ALPHA,BETA,N,M,L,NFIX,ASYM)
C ----------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION IPIV(:),WORK(:),Q(:,:),RPT(NR),VPT(NR)
      DIMENSION BETA(1),ALPHA(NR)
      DIMENSION WORK2(:)
      ALLOCATABLE IPIV,WORK,Q,WORK2


      ALLOCATE(IPIV(NR))
      ALLOCATE(WORK(NR))
      ALLOCATE(WORK2(NR))
      ALLOCATE(Q(NR,NR))

      DO I=1,NR
        DO J=1,NR
          RS = MIN( RPT(I),RPT(J))
          RL = MAX( RPT(I),RPT(J))

          Q(I,J)=0.0D0
          DO K=0,N-1
            Q(I,J)=Q(I,J)+BETA(K+1)*(RS/RL)**K
          ENDDO
          Q(I,J)=Q(I,J)*(RL)**(-(M+1))
        ENDDO
      ENDDO
      CALL DCOPY (NR,VPT,1,WORK2,1)
      CALL DGESV (NR, 1, Q, NR, IPIV, VPT, NR, INFO)
      CALL DCOPY (NR,VPT,1,ALPHA,1)
      CALL DCOPY (NR,WORK2,1,VPT,1)
      DEALLOCATE(IPIV,WORK,WORK2,Q)
      END
C ----------------------------------------------------
      DOUBLE PRECISION FUNCTION V_RKHS(R,N,M,L,NR,BETA,RPT,ALPHA)
C ----------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION BETA(1),RPT(NR),ALPHA(NR)

      V = 0.0D0
      DO I=1,NR
        RS = MIN(R,RPT(I))
        RL = MAX(R,RPT(I))
        Q  = 0.0D0
        DO K=1,N
          Q = Q+ BETA(K)*(RS/RL)**(K-1)
        ENDDO
        Q = Q*(1.0D0/RL)**(L*(M+1))
        V = V + Q*ALPHA(I)
      ENDDO
      V_RKHS = V
      END
C ----------------------------------------------------
      SUBROUTINE CALCPOL(NANGLE,THETA,P,PN)
C ----------------------------------------------------
      IMPLICIT NONE
      INTEGER I, NANGLE, MU
      DOUBLE PRECISION THETA,PI
      DOUBLE PRECISION P(NANGLE)
      DOUBLE PRECISION PN(NANGLE)

      P(1)=1.D0
      P(2)=THETA
      DO I=1,NANGLE-2
        P(I+2)=(P(I+1)*P(2)*DBLE(2*I+1)-P(I)*DBLE(I))/DBLE(I+1)
      ENDDO

      DO I=0,NANGLE-1
        PN(I+1)=P(I+1)*(DBLE(2*I+1)/(4.D0*ACOS(-1.D0)))**0.5D0
      ENDDO
      END
C ----------------------------------------------------
       SUBROUTINE GET_RADIAL(R, RR, V, NANG, NRA, RPT, WT,
     O                       ALPHA, BETA, N, M, L, NRMAX)
C ----------------------------------------------------
       IMPLICIT DOUBLE PRECISION(A-H,O-Z)
       DIMENSION V(NRMAX, NANG), NRA(NANG), ALPHA(NRMAX, NANG)
       DIMENSION RPT(NRMAX,NANG), BETA(N), RR(NANG), WT(NANG),P(NANG)
       LOGICAL FLAT

       DO I=1,NANG
         RR(I) = V_RKHS(R,N,M,L,NRA(I),BETA,RPT(1,I),ALPHA(1,I))
       ENDDO
       END
      double precision function fswitch(r,a,b)
      implicit real*8(a-h,o-z)
      parameter (pi2=1.570796326794897)
      if(r.le.a) then 
         fswitch=0d0
      elseif(r.ge.b) then
         fswitch=1d0
      else
         x=((r-b)-(a-r))/(b-a)
         fswitch=0.5d0+0.25d0*dsin(pi2*x)*(3d0-(dsin(pi2*x))**2)
      endif
      end  

