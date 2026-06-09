      SUBROUTINE  VINIT(II,RM,RESULT)
C  Copyright (C) 2019 J. M. Hutson & C. R. Le Sueur
C  Distributed under the GNU General Public License, version 3
C -------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      DOUBLE PRECISION, ALLOCATABLE:: ALPHA0(:),ALPHA2(:),BETA(:)
      DOUBLE PRECISION, ALLOCATABLE:: V0(:), V2(:), RR(:)
      INTEGER, ALLOCATABLE:: NRA(:)
      LOGICAL LFIRST
      CHARACTER(80) FILNAM
      DIMENSION LAM(2)
      DATA LAM /0,2/

      DATA LFIRST /.TRUE./
      DATA FILNAM /'data/pot-Sr_O2.data'/


      IF (II.NE.1) GOTO 15

C  READ THE FILE data/pot-Mg_NH.data
C  THE STRUCTURE OF FILE IS FOLLOWING:
C  NANGLE  ICASEM: NUMBER OF ANGULAR CUTS,
C          ICASE = IF THE MORPHING IS SWITCHED ON
C  ANGLE1
C      XX1  EINT1 EINT2  EEX  EIND  EDISP
C      ...
C      XXNR1   EINT1 EINT2  EEX  EIND  EDISP
C  ANGLE2
C  ETC...

      IF (.NOT.LFIRST) DEALLOCATE(V0,V2,ALPHA0,ALPHA2,RR,
     1                            BETA)
      LFIRST=.FALSE.
      OPEN(9,FILE=FILNAM,FORM='FORMATTED',STATUS='OLD',ERR=2010)
      READ(9,*) NR   
C  DETERMINE MAXIMUM NRA(I) TO ALLOCATE PROPERLY ALPHA
C  AND VPT ARRAYS
      NANGS=2
      ALLOCATE(V0(NR))
      ALLOCATE(V2(NR))
      ALLOCATE(ALPHA0(NR))
      ALLOCATE(ALPHA2(NR))
      ALLOCATE(RR(NR))

        DO J=1,NR 
           READ(9,*)  RR(J), EL, ET
           V0(J) = (EL + 2.0d0*ET)/3.0d0     
           V2(J) = 2.0d0*(EL - ET)/3.0d0    
        ENDDO
      CLOSE(9)

      LMAX = 2       

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
      CALL RK_INIT(NR, RR, V0, ALPHA0, BETA, NP, MP, LP, NFIX, ASYM)
      CALL RK_INIT(NR, RR, V2, ALPHA2, BETA, NP, MP, LP, NFIX, ASYM)

  15  RETURN

C -------------------------------------------------------
      ENTRY VSTAR(II,X,TOTAL)
C -------------------------------------------------------

      TOTAL = 0.0D0
      IF (LAM(II).EQ.0) THEN
       TOTAL = VRKHS(X, Np, Mp, Lp, NR, BETA, RR, ALPHA0)  
      ENDIF
 
      IF (LAM(II).EQ.2) THEN
       TOTAL = VRKHS(X, Np, Mp, Lp, NR, BETA, RR, ALPHA2)  
      ENDIF

      IF (II.GT.LMAX) THEN
       WRITE (*,*) ' lambda = ',II,' not good ... ' 
       GOTO 17
      ENDIF
   16 RETURN

   17 STOP 'WRONG LAMBDA IN VSTAR ROUTINE '
      RETURN
C --------------------------------------------------------
C  DUMMY ROUTINES
C --------------------------------------------------------
      ENTRY VSTAR1(II,X,TOTAL)
C --------------------------------------------------------
C --------------------------------------------------------
      ENTRY VSTAR2(II,X,TOTAL)
C --------------------------------------------------------
      WRITE(6,*) ' VSTAR: DERIVATIVES NOT IMPLEMENTED'
      RETURN
 2010 WRITE(6,*) ' *** ERROR: '//TRIM(FILNAM)//' NOT FOUND'
      STOP
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
      DOUBLE PRECISION FUNCTION VRKHS(R,N,M,L,NR,BETA,RPT,ALPHA)
C ----------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
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
      VRKHS = V

      END
