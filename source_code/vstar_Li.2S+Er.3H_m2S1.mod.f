      SUBROUTINE VINIT(I,RM,EPSIL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C-----------------------------------------------------------------------
C  DATE (LAST UPDATE): 22/08/14       STATUS: FINISHED                 |
C  AUTHOR: MAYKEL LEONARDO GONZALEZ MARTINEZ                           |
C-----------------------------------------------------------------------
C      THIS SUBROUTINE INCLUDES TO BOUND & MOLSCAT, VIA THE VSTAR      |
C      MECHANISM, THE V^S_0--10 TERMS FOR Li+Er                        |
C       (AT INPUT, X IN ANGSTROMS.  AT OUTPUT, ENERGIES IN CM-1)       |
C-----------------------------------------------------------------------
C USES 
C PARAMETERS
      PARAMETER (MXL=100,IXMX=200,IEXMX=200,MXDIM=3)
      PARAMETER (A02AA=0.52917726D0,EV2EH=3.674932379D-2)
      PARAMETER (BETA=2.D0*A02AA,X0=12.D0*A02AA)
C 'EXTRA' VARIABLES:
      DIMENSION LAMBDA(MXL),NPTS(MXDIM),NTERM(MXL),NPOWER(IXMX)
      DIMENSION A(IXMX),E(IEXMX)
      LOGICAL LVRTP
C NAMELISTS
C . POTL
      NAMELIST /POTL/ RM,EPSIL,MXLAM,MXSYM,LAMBDA,NPOTL,NTERM,
     &   A,NPOWER,E,CFLAG,LVRTP,NPT,NPS,
     &   IHOMO,ICNSYM,LMAX,L1MAX,L2MAX,MMAX,IVMIN,IVMAX,
     &   NPTS,IHOMO2,ICNSY2
C
      IF (I.NE.1) RETURN
      REWIND(5)
      READ(5,POTL)
C INITIALISES THE V^{S = 1/2}_{K = 0, 2,..., 8} 
      CALL VL0_INIT
      CALL VL2_INIT
      CALL VL4_INIT
      CALL VL6_INIT
      CALL VL8_INIT
      CALL VL10_INIT
C INITIALISES THE V^{S = 3/2}_{K = 0, 2,..., 8} 
      CALL VH0_INIT
      CALL VH2_INIT
      CALL VH4_INIT
      CALL VH6_INIT
      CALL VH8_INIT
      CALL VH10_INIT
C
      RETURN
C-----------------------------------------------------------------------
      ENTRY VSTAR(I,X,SUM)
C-----------------------------------------------------------------------
      XA0=X/A02AA ! CONVERTS R/X FROM ANGSTROMS TO BOHRS
      LMBD=LAMBDA(I)
      IF (LMBD.EQ.-1) THEN
C CP1 ELECTRON SPIN-CP2 ELECTRON SPIN R-DEPENDENCE
         SUM=XA0**(-3) ! X HAS BEEN CONVERTED FROM ANGSTROMS TO BOHRS
         RETURN
      ELSE IF (LMBD.GE.0) THEN
! NOTE THE CHANGE IN THE INTERPRETATION OF LAMBDA!: MLGM 07/07/2014
C FOR SEVERAL SPIN MULTIPLICITIES, LAMBDA IS MAPPED AS:
C  LAMBDA = (I2S1-1)*100 + K, WITH
C  (1) I2S1 => ~ (I2S1+1)TH SPIN MULTIPLICITY (MAX. 10 DIFFERENT I2S1)
C  (2) K => 0, 2, ..., 2*L2 TERM IN THE POTENTIAL EXPANSION, GIVEN I2S1
         I2S1=INT(LMBD/100)+1
         K=LMBD-100*(I2S1-1)
C USES VK FROM Li+Er---EXCEPT FOR V0, FROM Li+Yb
         IF (I2S1.EQ.1) THEN
            SELECT CASE (K)
             CASE (0)
               CALL VL0_EVAL(X,VL0)
               SUM=VL0
             CASE (2)
               CALL VL2_EVAL(X,VL2)
               SUM=VL2
             CASE (4)
               CALL VL4_EVAL(X,VL4)
               SUM=VL4*SWITCH(BETA,X0,X)
             CASE (6)
               CALL VL6_EVAL(X,VL6)
               SUM=VL6*SWITCH(BETA,X0,X)
             CASE (8)
               CALL VL8_EVAL(X,VL8)
               SUM=VL8*SWITCH(BETA,X0,X)
             CASE (10)
               CALL VL10_EVAL(X,VL10)
               SUM=VL10*SWITCH(BETA,X0,X)
            END SELECT
         ELSE IF (I2S1.EQ.2) THEN
            SELECT CASE (K)
             CASE (0)
               CALL VH0_EVAL(X,VH0)
               SUM=VH0
             CASE (2)
               CALL VH2_EVAL(X,VH2)
               SUM=VH2
             CASE (4)
               CALL VH4_EVAL(X,VH4)
               SUM=VH4*SWITCH(BETA,X0,X)
             CASE (6)
               CALL VH6_EVAL(X,VH6)
               SUM=VH6*SWITCH(BETA,X0,X)
             CASE (8)
               CALL VH8_EVAL(X,VH8)
               SUM=VH8*SWITCH(BETA,X0,X)
             CASE (10)
               CALL VH10_EVAL(X,VH10)
               SUM=VH10*SWITCH(BETA,X0,X)
            END SELECT
         ELSE
            PRINT*, '* * * ERROR * * * SUBROUTINE VSTAR'
            STOP    ' ONLY 2S + 1 = 2, 4 INCLUDED IN THIS VERSION'
         END IF
C
         RETURN
      END IF
C DUMMY ROUTINES
      ENTRY VSTAR1(I,X,SUM)
      ENTRY VSTAR2(I,X,SUM)
      PRINT*, '* * * ERROR * * * SUBROUTINE VSTAR'
      STOP    ' VSTAR: DERIVATIVES NOT IMPLEMENTED'
C
      RETURN
      END SUBROUTINE VINIT
C******************************* VL0_INIT ******************************
      SUBROUTINE VL0_INIT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C-----------------------------------------------------------------------
C  DATE (LAST UPDATE): 07/07/14       STATUS: FINISHED, UNDER TEST     |
C  AUTHOR: MAYKEL LEONARDO GONZALEZ MARTINEZ                           |
C-----------------------------------------------------------------------
C  THIS SUBROUTINE EVALUATES THE |
C-----------------------------------------------------------------------
C USES VMORSE, RKHS_INIT, V_RKHS
C PARAMETERS
      PARAMETER (A02AA=0.52917726D0)
C 'EXTRA' VARIABLES:
      INTEGER S
      CHARACTER FILE_IN*27,LINE*100
      LOGICAL ASYMP,CSN_PRINT
      ALLOCATABLE RR(:),VV(:),CSN(:),ANMK(:),ALPHA(:)
C
C OPENS, READS & CLOSES INPUT/PARAMETERS FILE
      FILE_IN='V0.LiEr.MCSCF.S_0.5+1.5.dat'
      OPEN(9,FILE=FILE_IN,STATUS='OLD',ERR=2010)
      WRITE(6,'(/" V^S=1/2_k=0 FROM FILE: ",A30)') FILE_IN
C . GLOBAL PARAMETERS FOR THE PEC
      READ(9,'(A100)') LINE
      WRITE(6,'(" VL0_INIT: ",A100)') LINE
      READ(9,*) SCL,NPT,NPTMX,ESHIFT
      WRITE(6,'(13X,F8.5,2(2X,I4),2X,F10.6)') SCL,NPT,NPTMX,ESHIFT
C . INPUT PARAMETERS FOR RKHS INTER-/EXTRAPOLATION 
      READ(9,'(A100)') LINE
      WRITE(6,'(" VL0_INIT: ",A100)') LINE
      READ(9,*) N,M,S,ASYMP,NCSN,CSN_PRINT
      WRITE(6,'(13X,3(I2,2X),L2,2X,I3,2X,L2)')N,M,S,ASYMP,NCSN,CSN_PRINT
      IF (ASYMP) THEN
         ALLOCATE (CSN(NCSN))
         ALLOCATE (RR(NPT+NCSN))
         READ(9,'(A100)') LINE
         WRITE(6,'(" VL0_INIT: ",A100)') LINE
         READ(9,*) RA,CSN,RR(NPT+1:NPT+NCSN)
         WRITE(6,'(13X,F5.2,2(2X,F14.2),2(2X,F6.2))') RA,CSN,
     &                                                RR(NPT+1:NPT+NCSN)
      ELSE
         read(9,*)
         read(9,*)
         ALLOCATE (RR(NPTMX))
         NCSN=0
         RA=0.D0
      END IF
      ALLOCATE (VV(NPTMX))
      read(9,*)
      DO I=1,NPTMX
         READ(9,*) RR(I),VV(I)
         RR(I)=RR(I)*A02AA ! CONVERTS FROM a0 TO AA
         VV(I)=VV(I)*SCL
         write(20,*) rr(i),vv(i)
      END DO
      CLOSE(9)
C COMPUTES RKHS COEFFICIENTS ANMK AND ALPHA
      NN=NPT+NCSN
      ALLOCATE (ANMK(N))
      ALLOCATE (ALPHA(NN))
      CALL RKHS_INIT(RR,VV,NPT,N,M,S,ASYMP,CSN,NCSN,RA,ANMK,ALPHA,
     & CSN_PRINT)
C
      IF (ASYMP) DEALLOCATE (CSN)
C
      RETURN
C------------------------------- VL0_EVAL ------------------------------
      ENTRY VL0_EVAL(X,VL0)
C-----------------------------------------------------------------------
      VL0=V_RKHS(X,RR,NN,N,M,S,ANMK,ALPHA)
C
      RETURN
 2010 PRINT*, " ERROR IN VL0_INIT: COULD NOT OPEN FILE ", FILE_IN
      STOP
      END SUBROUTINE VL0_INIT
C******************************* VH0_INIT ******************************
      SUBROUTINE VH0_INIT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C-----------------------------------------------------------------------
C  DATE (LAST UPDATE): 07/07/14       STATUS: FINISHED, UNDER TEST     |
C  AUTHOR: MAYKEL LEONARDO GONZALEZ MARTINEZ                           |
C-----------------------------------------------------------------------
C  THIS SUBROUTINE EVALUATES THE |
C-----------------------------------------------------------------------
C USES VMORSE, RKHS_INIT, V_RKHS
C PARAMETERS
      PARAMETER (A02AA=0.52917726D0)
C 'EXTRA' VARIABLES:
      INTEGER S
      CHARACTER FILE_IN*27,LINE*100
      LOGICAL ASYMP,CSN_PRINT
      ALLOCATABLE RR(:),VV(:),CSN(:),ANMK(:),ALPHA(:)
C
C OPENS, READS & CLOSES INPUT/PARAMETERS FILE
      FILE_IN='V0.LiEr.MCSCF.S_0.5+1.5.dat'
      OPEN(9,FILE=FILE_IN,STATUS='OLD',ERR=2011)
      WRITE(6,'(/" V^S=3/2_k=0 FROM FILE: ",A30)') FILE_IN
C . GLOBAL PARAMETERS FOR THE PEC
      READ(9,'(A100)') LINE
      WRITE(6,'(" VH0_INIT: ",A100)') LINE
      READ(9,*) SCL,NPT,NPTMX,ESHIFT
      WRITE(6,'(13X,F8.5,2(2X,I4),2X,F10.6)') SCL,NPT,NPTMX,ESHIFT
C . INPUT PARAMETERS FOR RKHS INTER-/EXTRAPOLATION 
      READ(9,'(A100)') LINE
      WRITE(6,'(" VH0_INIT: ",A100)') LINE
      READ(9,*) N,M,S,ASYMP,NCSN,CSN_PRINT
      WRITE(6,'(13X,3(I2,2X),L2,2X,I3,2X,L2)')N,M,S,ASYMP,NCSN,CSN_PRINT
      IF (ASYMP) THEN
         ALLOCATE (CSN(NCSN))
         ALLOCATE (RR(NPT+NCSN))
         READ(9,'(A100)') LINE
         WRITE(6,'(" VH0_INIT: ",A100)') LINE
         READ(9,*) RA,CSN,RR(NPT+1:NPT+NCSN)
         WRITE(6,'(13X,F5.2,2(2X,F14.2),2(2X,F6.2))') RA,CSN,
     &                                                RR(NPT+1:NPT+NCSN)
      ELSE
         read(9,*)
         read(9,*)
         ALLOCATE (RR(NPTMX))
         NCSN=0
         RA=0.D0
      END IF
      ALLOCATE (VV(NPTMX))
      read(9,*)
      DO I=1,NPTMX
         READ(9,*) RR(I),DUMMY,VV(I)
         VV(I)=VV(I)*SCL
         RR(I)=RR(I)*A02AA ! CONVERTS FROM a0 TO AA
         write(21,*) rr(i),vv(i)
      END DO
      CLOSE(9)
C COMPUTES RKHS COEFFICIENTS ANMK AND ALPHA
      NN=NPT+NCSN
      ALLOCATE (ANMK(N))
      ALLOCATE (ALPHA(NN))
      CALL RKHS_INIT(RR,VV,NPT,N,M,S,ASYMP,CSN,NCSN,RA,ANMK,ALPHA,
     & CSN_PRINT)
C
      IF (ASYMP) DEALLOCATE (CSN)
C
      RETURN
C------------------------------- VH0_EVAL ------------------------------
      ENTRY VH0_EVAL(X,VH0)
C-----------------------------------------------------------------------
      VH0=V_RKHS(X,RR,NN,N,M,S,ANMK,ALPHA)
C
      RETURN
 2011 PRINT*, " ERROR IN VH0_INIT: COULD NOT OPEN FILE ", FILE_IN
      STOP
      END SUBROUTINE VH0_INIT
C******************************* VL2_INIT ******************************
      SUBROUTINE VL2_INIT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C-----------------------------------------------------------------------
C  DATE (LAST UPDATE): 07/07/14       STATUS: FINISHED, UNDER TEST     |
C  AUTHOR: MAYKEL LEONARDO GONZALEZ MARTINEZ                           |
C-----------------------------------------------------------------------
C  THIS SUBROUTINE EVALUATES THE |
C-----------------------------------------------------------------------
C USES VMORSE, RKHS_INIT, V_RKHS
C PARAMETERS
      PARAMETER (A02AA=0.52917726D0)
C 'EXTRA' VARIABLES:
      INTEGER S
      CHARACTER FILE_IN*27,LINE*100
      LOGICAL ASYMP,CSN_PRINT
      ALLOCATABLE RR(:),VV(:),CSN(:),ANMK(:),ALPHA(:)
C
C OPENS, READS & CLOSES INPUT/PARAMETERS FILE
      FILE_IN='V2.LiEr.MCSCF.S_0.5+1.5.dat'
      OPEN(9,FILE=FILE_IN,STATUS='OLD',ERR=2012)
      WRITE(6,'(/" V^S=1/2_k=2 FROM FILE: ",A30)') FILE_IN
C . GLOBAL PARAMETERS FOR THE PEC
      READ(9,'(A100)') LINE
      WRITE(6,'(" VL2_INIT: ",A100)') LINE
      READ(9,*) SCL,NPT,NPTMX,ESHIFT
      WRITE(6,'(13X,F8.5,2(2X,I4),2X,F10.6)') SCL,NPT,NPTMX,ESHIFT
C . INPUT PARAMETERS FOR RKHS INTER-/EXTRAPOLATION 
      READ(9,'(A100)') LINE
      WRITE(6,'(" VL2_INIT: ",A100)') LINE
      READ(9,*) N,M,S,ASYMP,NCSN,CSN_PRINT
      WRITE(6,'(13X,3(I2,2X),L2,2X,I3,2X,L2)')N,M,S,ASYMP,NCSN,CSN_PRINT
      IF (ASYMP) THEN
         ALLOCATE (CSN(NCSN))
         ALLOCATE (RR(NPT+NCSN))
         READ(9,'(A100)') LINE
         WRITE(6,'(" VL2_INIT: ",A100)') LINE
         READ(9,*) RA,CSN,RR(NPT+1:NPT+NCSN)
         WRITE(6,'(13X,F5.2,2(2X,F14.2),2(2X,F6.2))') RA,CSN,
     &                                                RR(NPT+1:NPT+NCSN)
      ELSE
         read(9,*)
         read(9,*)
         ALLOCATE (RR(NPTMX))
         NCSN=0
         RA=0.D0
      END IF
      ALLOCATE (VV(NPTMX))
      read(9,*)
      DO I=1,NPTMX
         READ(9,*) RR(I),VV(I)
         RR(I)=RR(I)*A02AA ! CONVERTS FROM a0 TO AA
!         write(22,*) rr(i),vv(i)
      END DO
      CLOSE(9)
C COMPUTES RKHS COEFFICIENTS ANMK AND ALPHA
      NN=NPT+NCSN
      ALLOCATE (ANMK(N))
      ALLOCATE (ALPHA(NN))
      CALL RKHS_INIT(RR,VV,NPT,N,M,S,ASYMP,CSN,NCSN,RA,ANMK,ALPHA,
     & CSN_PRINT)
C
      IF (ASYMP) DEALLOCATE (CSN)
C
      RETURN
C------------------------------- VL2_EVAL ------------------------------
      ENTRY VL2_EVAL(X,VL2)
C-----------------------------------------------------------------------
      VL2=V_RKHS(X,RR,NN,N,M,S,ANMK,ALPHA)
C
      RETURN
 2012 PRINT*, " ERROR IN VL2_INIT: COULD NOT OPEN FILE ", FILE_IN
      STOP
      END SUBROUTINE VL2_INIT
C******************************* VH2_INIT ******************************
      SUBROUTINE VH2_INIT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C-----------------------------------------------------------------------
C  DATE (LAST UPDATE): 07/07/14       STATUS: FINISHED, UNDER TEST     |
C  AUTHOR: MAYKEL LEONARDO GONZALEZ MARTINEZ                           |
C-----------------------------------------------------------------------
C  THIS SUBROUTINE EVALUATES THE |
C-----------------------------------------------------------------------
C USES VMORSE, RKHS_INIT, V_RKHS
C PARAMETERS
      PARAMETER (A02AA=0.52917726D0)
C 'EXTRA' VARIABLES:
      INTEGER S
      CHARACTER FILE_IN*27,LINE*100
      LOGICAL ASYMP,CSN_PRINT
      ALLOCATABLE RR(:),VV(:),CSN(:),ANMK(:),ALPHA(:)
C
C OPENS, READS & CLOSES INPUT/PARAMETERS FILE
      FILE_IN='V2.LiEr.MCSCF.S_0.5+1.5.dat'
      OPEN(9,FILE=FILE_IN,STATUS='OLD',ERR=2013)
      WRITE(6,'(/" V^S=3/2_k=2 FROM FILE: ",A30)') FILE_IN
C . GLOBAL PARAMETERS FOR THE PEC
      READ(9,'(A100)') LINE
      WRITE(6,'(" VH2_INIT: ",A100)') LINE
      READ(9,*) SCL,NPT,NPTMX,ESHIFT
      WRITE(6,'(13X,F8.5,2(2X,I4),2X,F10.6)') SCL,NPT,NPTMX,ESHIFT
C . INPUT PARAMETERS FOR RKHS INTER-/EXTRAPOLATION 
      READ(9,'(A100)') LINE
      WRITE(6,'(" VH2_INIT: ",A100)') LINE
      READ(9,*) N,M,S,ASYMP,NCSN,CSN_PRINT
      WRITE(6,'(13X,3(I2,2X),L2,2X,I3,2X,L2)')N,M,S,ASYMP,NCSN,CSN_PRINT
      IF (ASYMP) THEN
         ALLOCATE (CSN(NCSN))
         ALLOCATE (RR(NPT+NCSN))
         READ(9,'(A100)') LINE
         WRITE(6,'(" VH2_INIT: ",A100)') LINE
         READ(9,*) RA,CSN,RR(NPT+1:NPT+NCSN)
         WRITE(6,'(13X,F5.2,2(2X,F14.2),2(2X,F6.2))') RA,CSN,
     &                                                RR(NPT+1:NPT+NCSN)
      ELSE
         read(9,*)
         read(9,*)
         ALLOCATE (RR(NPTMX))
         NCSN=0
         RA=0.D0
      END IF
      ALLOCATE (VV(NPTMX))
      read(9,*)
      DO I=1,NPTMX
         READ(9,*) RR(I),DUMMY,VV(I)
         RR(I)=RR(I)*A02AA ! CONVERTS FROM a0 TO AA
!         write(23,*) rr(i),vv(i)
      END DO
      CLOSE(9)
C COMPUTES RKHS COEFFICIENTS ANMK AND ALPHA
      NN=NPT+NCSN
      ALLOCATE (ANMK(N))
      ALLOCATE (ALPHA(NN))
      CALL RKHS_INIT(RR,VV,NPT,N,M,S,ASYMP,CSN,NCSN,RA,ANMK,ALPHA,
     & CSN_PRINT)
C
      IF (ASYMP) DEALLOCATE (CSN)
C
      RETURN
C------------------------------- VH2_EVAL ------------------------------
      ENTRY VH2_EVAL(X,VH2)
C-----------------------------------------------------------------------
      VH2=V_RKHS(X,RR,NN,N,M,S,ANMK,ALPHA)
C
      RETURN
 2013 PRINT*, " ERROR IN VH2_INIT: COULD NOT OPEN FILE ", FILE_IN
      STOP
      END SUBROUTINE VH2_INIT
C******************************* VL4_INIT ******************************
      SUBROUTINE VL4_INIT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C-----------------------------------------------------------------------
C  DATE (LAST UPDATE): 07/07/14       STATUS: FINISHED, UNDER TEST     |
C  AUTHOR: MAYKEL LEONARDO GONZALEZ MARTINEZ                           |
C-----------------------------------------------------------------------
C  THIS SUBROUTINE EVALUATES THE |
C-----------------------------------------------------------------------
C USES VMORSE, RKHS_INIT, V_RKHS
C PARAMETERS
      PARAMETER (A02AA=0.52917726D0)
C 'EXTRA' VARIABLES:
      INTEGER S
      CHARACTER FILE_IN*26,LINE*100
      LOGICAL ASYMP,CSN_PRINT
      ALLOCATABLE RR(:),VV(:),CSN(:),ANMK(:),ALPHA(:)
C
C OPENS, READS & CLOSES INPUT/PARAMETERS FILE
      FILE_IN='V4-10.LiEr.MCSCF.S_0.5.dat'
      OPEN(9,FILE=FILE_IN,STATUS='OLD',ERR=2014)
      WRITE(6,'(/" V^S=1/2_k=4 FROM FILE: ",A30)') FILE_IN
C . GLOBAL PARAMETERS FOR THE PEC
      READ(9,'(A100)') LINE
      WRITE(6,'(" VL4_INIT: ",A100)') LINE
      READ(9,*) SCL,NPT,NPTMX,ESHIFT
      WRITE(6,'(13X,F8.5,2(2X,I4),2X,F10.6)') SCL,NPT,NPTMX,ESHIFT
C . INPUT PARAMETERS FOR RKHS INTER-/EXTRAPOLATION 
      READ(9,'(A100)') LINE
      WRITE(6,'(" VL4_INIT: ",A100)') LINE
      READ(9,*) N,M,S,ASYMP,NCSN,CSN_PRINT
      WRITE(6,'(13X,3(I2,2X),L2,2X,I3,2X,L2)')N,M,S,ASYMP,NCSN,CSN_PRINT
      IF (ASYMP) THEN
         ALLOCATE (CSN(NCSN))
         ALLOCATE (RR(NPT+NCSN))
         READ(9,'(A100)') LINE
         WRITE(6,'(" VL4_INIT: ",A100)') LINE
         READ(9,*) RA,CSN,RR(NPT+1:NPT+NCSN)
         WRITE(6,'(13X,F5.2,2(2X,F14.2),2(2X,F6.2))') RA,CSN,
     &                                                RR(NPT+1:NPT+NCSN)
      ELSE
         read(9,*)
         read(9,*)
         ALLOCATE (RR(NPTMX))
         NCSN=0
         RA=0.D0
      END IF
      ALLOCATE (VV(NPTMX))
      read(9,*)
      DO I=1,NPTMX
         READ(9,*) RR(I),VV(I)
         RR(I)=RR(I)*A02AA ! CONVERTS FROM a0 TO AA
!         write(24,*) rr(i),vv(i)
      END DO
      CLOSE(9)
C COMPUTES RKHS COEFFICIENTS ANMK AND ALPHA
      NN=NPT+NCSN
      ALLOCATE (ANMK(N))
      ALLOCATE (ALPHA(NN))
      CALL RKHS_INIT(RR,VV,NPT,N,M,S,ASYMP,CSN,NCSN,RA,ANMK,ALPHA,
     & CSN_PRINT)
C
      IF (ASYMP) DEALLOCATE (CSN)
C
      RETURN
C------------------------------- VL4_EVAL ------------------------------
      ENTRY VL4_EVAL(X,VL4)
C-----------------------------------------------------------------------
      VL4=V_RKHS(X,RR,NN,N,M,S,ANMK,ALPHA)
C
      RETURN
 2014 PRINT*, " ERROR IN VL4_INIT: COULD NOT OPEN FILE ", FILE_IN
      STOP
      END SUBROUTINE VL4_INIT
C******************************* VH4_INIT ******************************
      SUBROUTINE VH4_INIT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C-----------------------------------------------------------------------
C  DATE (LAST UPDATE): 07/07/14       STATUS: FINISHED, UNDER TEST     |
C  AUTHOR: MAYKEL LEONARDO GONZALEZ MARTINEZ                           |
C-----------------------------------------------------------------------
C  THIS SUBROUTINE EVALUATES THE |
C-----------------------------------------------------------------------
C USES VMORSE, RKHS_INIT, V_RKHS
C PARAMETERS
      PARAMETER (A02AA=0.52917726D0)
C 'EXTRA' VARIABLES:
      INTEGER S
      CHARACTER FILE_IN*26,LINE*100
      LOGICAL ASYMP,CSN_PRINT
      ALLOCATABLE RR(:),VV(:),CSN(:),ANMK(:),ALPHA(:)
C
C OPENS, READS & CLOSES INPUT/PARAMETERS FILE
      FILE_IN='V4-10.LiEr.MCSCF.S_1.5.dat'
      OPEN(9,FILE=FILE_IN,STATUS='OLD',ERR=2015)
      WRITE(6,'(/" V^S=3/2_k=4 FROM FILE: ",A30)') FILE_IN
C . GLOBAL PARAMETERS FOR THE PEC
      READ(9,'(A100)') LINE
      WRITE(6,'(" VH4_INIT: ",A100)') LINE
      READ(9,*) SCL,NPT,NPTMX,ESHIFT
      WRITE(6,'(13X,F8.5,2(2X,I4),2X,F10.6)') SCL,NPT,NPTMX,ESHIFT
C . INPUT PARAMETERS FOR RKHS INTER-/EXTRAPOLATION 
      READ(9,'(A100)') LINE
      WRITE(6,'(" VH4_INIT: ",A100)') LINE
      READ(9,*) N,M,S,ASYMP,NCSN,CSN_PRINT
      WRITE(6,'(13X,3(I2,2X),L2,2X,I3,2X,L2)')N,M,S,ASYMP,NCSN,CSN_PRINT
      IF (ASYMP) THEN
         ALLOCATE (CSN(NCSN))
         ALLOCATE (RR(NPT+NCSN))
         READ(9,'(A100)') LINE
         WRITE(6,'(" VH4_INIT: ",A100)') LINE
         READ(9,*) RA,CSN,RR(NPT+1:NPT+NCSN)
         WRITE(6,'(13X,F5.2,2(2X,F14.2),2(2X,F6.2))') RA,CSN,
     &                                                RR(NPT+1:NPT+NCSN)
      ELSE
         read(9,*)
         read(9,*)
         ALLOCATE (RR(NPTMX))
         NCSN=0
         RA=0.D0
      END IF
      ALLOCATE (VV(NPTMX))
      read(9,*)
      DO I=1,NPTMX
         READ(9,*) RR(I),VV(I)
         RR(I)=RR(I)*A02AA ! CONVERTS FROM a0 TO AA
!         write(25,*) rr(i),vv(i)
      END DO
      CLOSE(9)
C COMPUTES RKHS COEFFICIENTS ANMK AND ALPHA
      NN=NPT+NCSN
      ALLOCATE (ANMK(N))
      ALLOCATE (ALPHA(NN))
      CALL RKHS_INIT(RR,VV,NPT,N,M,S,ASYMP,CSN,NCSN,RA,ANMK,ALPHA,
     & CSN_PRINT)
C
      IF (ASYMP) DEALLOCATE (CSN)
C
      RETURN
C------------------------------- VH4_EVAL ------------------------------
      ENTRY VH4_EVAL(X,VH4)
C-----------------------------------------------------------------------
      VH4=V_RKHS(X,RR,NN,N,M,S,ANMK,ALPHA)
C
      RETURN
 2015 PRINT*, " ERROR IN VH4_INIT: COULD NOT OPEN FILE ", FILE_IN
      STOP
      END SUBROUTINE VH4_INIT
C******************************* VL6_INIT ******************************
      SUBROUTINE VL6_INIT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C-----------------------------------------------------------------------
C  DATE (LAST UPDATE): 07/07/14       STATUS: FINISHED, UNDER TEST     |
C  AUTHOR: MAYKEL LEONARDO GONZALEZ MARTINEZ                           |
C-----------------------------------------------------------------------
C  THIS SUBROUTINE EVALUATES THE |
C-----------------------------------------------------------------------
C USES VMORSE, RKHS_INIT, V_RKHS
C PARAMETERS
      PARAMETER (A02AA=0.52917726D0)
C 'EXTRA' VARIABLES:
      INTEGER S
      CHARACTER FILE_IN*26,LINE*100
      LOGICAL ASYMP,CSN_PRINT
      ALLOCATABLE RR(:),VV(:),CSN(:),ANMK(:),ALPHA(:)
C
C OPENS, READS & CLOSES INPUT/PARAMETERS FILE
      FILE_IN='V4-10.LiEr.MCSCF.S_0.5.dat'
      OPEN(9,FILE=FILE_IN,STATUS='OLD',ERR=2016)
      WRITE(6,'(/" V^S=1/2_k=6 FROM FILE: ",A30)') FILE_IN
C . GLOBAL PARAMETERS FOR THE PEC
      READ(9,'(A100)') LINE
      WRITE(6,'(" VL6_INIT: ",A100)') LINE
      READ(9,*) SCL,NPT,NPTMX,ESHIFT
      WRITE(6,'(13X,F8.5,2(2X,I4),2X,F10.6)') SCL,NPT,NPTMX,ESHIFT
C . INPUT PARAMETERS FOR RKHS INTER-/EXTRAPOLATION 
      READ(9,'(A100)') LINE
      WRITE(6,'(" VL6_INIT: ",A100)') LINE
      READ(9,*) N,M,S,ASYMP,NCSN,CSN_PRINT
      WRITE(6,'(13X,3(I2,2X),L2,2X,I3,2X,L2)')N,M,S,ASYMP,NCSN,CSN_PRINT
      IF (ASYMP) THEN
         ALLOCATE (CSN(NCSN))
         ALLOCATE (RR(NPT+NCSN))
         READ(9,'(A100)') LINE
         WRITE(6,'(" VL6_INIT: ",A100)') LINE
         READ(9,*) RA,CSN,RR(NPT+1:NPT+NCSN)
         WRITE(6,'(13X,F5.2,2(2X,F14.2),2(2X,F6.2))') RA,CSN,
     &                                                RR(NPT+1:NPT+NCSN)
      ELSE
         read(9,*)
         read(9,*)
         ALLOCATE (RR(NPTMX))
         NCSN=0
         RA=0.D0
      END IF
      ALLOCATE (VV(NPTMX))
      read(9,*)
      DO I=1,NPTMX
         READ(9,*) RR(I),DUMMY,VV(I)
         RR(I)=RR(I)*A02AA ! CONVERTS FROM a0 TO AA
!         write(26,*) rr(i),vv(i)
      END DO
      CLOSE(9)
C COMPUTES RKHS COEFFICIENTS ANMK AND ALPHA
      NN=NPT+NCSN
      ALLOCATE (ANMK(N))
      ALLOCATE (ALPHA(NN))
      CALL RKHS_INIT(RR,VV,NPT,N,M,S,ASYMP,CSN,NCSN,RA,ANMK,ALPHA,
     & CSN_PRINT)
C
      IF (ASYMP) DEALLOCATE (CSN)
C
      RETURN
C------------------------------- VL6_EVAL ------------------------------
      ENTRY VL6_EVAL(X,VL6)
C-----------------------------------------------------------------------
      VL6=V_RKHS(X,RR,NN,N,M,S,ANMK,ALPHA)
C
      RETURN
 2016 PRINT*, " ERROR IN VL6_INIT: COULD NOT OPEN FILE ", FILE_IN
      STOP
      END SUBROUTINE VL6_INIT
C******************************* VH6_INIT ******************************
      SUBROUTINE VH6_INIT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C-----------------------------------------------------------------------
C  DATE (LAST UPDATE): 07/07/14       STATUS: FINISHED, UNDER TEST     |
C  AUTHOR: MAYKEL LEONARDO GONZALEZ MARTINEZ                           |
C-----------------------------------------------------------------------
C  THIS SUBROUTINE EVALUATES THE |
C-----------------------------------------------------------------------
C USES VMORSE, RKHS_INIT, V_RKHS
C PARAMETERS
      PARAMETER (A02AA=0.52917726D0)
C 'EXTRA' VARIABLES:
      INTEGER S
      CHARACTER FILE_IN*26,LINE*100
      LOGICAL ASYMP,CSN_PRINT
      ALLOCATABLE RR(:),VV(:),CSN(:),ANMK(:),ALPHA(:)
C
C OPENS, READS & CLOSES INPUT/PARAMETERS FILE
      FILE_IN='V4-10.LiEr.MCSCF.S_1.5.dat'
      OPEN(9,FILE=FILE_IN,STATUS='OLD',ERR=2017)
      WRITE(6,'(/" V^S=3/2_k=6 FROM FILE: ",A30)') FILE_IN
C . GLOBAL PARAMETERS FOR THE PEC
      READ(9,'(A100)') LINE
      WRITE(6,'(" VH6_INIT: ",A100)') LINE
      READ(9,*) SCL,NPT,NPTMX,ESHIFT
      WRITE(6,'(13X,F8.5,2(2X,I4),2X,F10.6)') SCL,NPT,NPTMX,ESHIFT
C . INPUT PARAMETERS FOR RKHS INTER-/EXTRAPOLATION 
      READ(9,'(A100)') LINE
      WRITE(6,'(" VH6_INIT: ",A100)') LINE
      READ(9,*) N,M,S,ASYMP,NCSN,CSN_PRINT
      WRITE(6,'(13X,3(I2,2X),L2,2X,I3,2X,L2)')N,M,S,ASYMP,NCSN,CSN_PRINT
      IF (ASYMP) THEN
         ALLOCATE (CSN(NCSN))
         ALLOCATE (RR(NPT+NCSN))
         READ(9,'(A100)') LINE
         WRITE(6,'(" VH6_INIT: ",A100)') LINE
         READ(9,*) RA,CSN,RR(NPT+1:NPT+NCSN)
         WRITE(6,'(13X,F5.2,2(2X,F14.2),2(2X,F6.2))') RA,CSN,
     &                                                RR(NPT+1:NPT+NCSN)
      ELSE
         read(9,*)
         read(9,*)
         ALLOCATE (RR(NPTMX))
         NCSN=0
         RA=0.D0
      END IF
      ALLOCATE (VV(NPTMX))
      read(9,*)
      DO I=1,NPTMX
         READ(9,*) RR(I),DUMMY,VV(I)
         RR(I)=RR(I)*A02AA ! CONVERTS FROM a0 TO AA
!         write(27,*) rr(i),vv(i)
      END DO
      CLOSE(9)
C COMPUTES RKHS COEFFICIENTS ANMK AND ALPHA
      NN=NPT+NCSN
      ALLOCATE (ANMK(N))
      ALLOCATE (ALPHA(NN))
      CALL RKHS_INIT(RR,VV,NPT,N,M,S,ASYMP,CSN,NCSN,RA,ANMK,ALPHA,
     & CSN_PRINT)
C
      IF (ASYMP) DEALLOCATE (CSN)
C
      RETURN
C------------------------------- VH6_EVAL ------------------------------
      ENTRY VH6_EVAL(X,VH6)
C-----------------------------------------------------------------------
      VH6=V_RKHS(X,RR,NN,N,M,S,ANMK,ALPHA)
C
      RETURN
 2017 PRINT*, " ERROR IN VH6_INIT: COULD NOT OPEN FILE ", FILE_IN
      STOP
      END SUBROUTINE VH6_INIT
C******************************* VL8_INIT ******************************
      SUBROUTINE VL8_INIT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C-----------------------------------------------------------------------
C  DATE (LAST UPDATE): 07/07/14       STATUS: FINISHED, UNDER TEST     |
C  AUTHOR: MAYKEL LEONARDO GONZALEZ MARTINEZ                           |
C-----------------------------------------------------------------------
C  THIS SUBROUTINE EVALUATES THE |
C-----------------------------------------------------------------------
C USES VMORSE, RKHS_INIT, V_RKHS
C PARAMETERS
      PARAMETER (A02AA=0.52917726D0)
C 'EXTRA' VARIABLES:
      INTEGER S
      CHARACTER FILE_IN*26,LINE*100
      LOGICAL ASYMP,CSN_PRINT
      ALLOCATABLE RR(:),VV(:),CSN(:),ANMK(:),ALPHA(:)
C
C OPENS, READS & CLOSES INPUT/PARAMETERS FILE
      FILE_IN='V4-10.LiEr.MCSCF.S_0.5.dat'
      OPEN(9,FILE=FILE_IN,STATUS='OLD',ERR=2018)
      WRITE(6,'(/" V^S=1/2_k=8 FROM FILE: ",A30)') FILE_IN
C . GLOBAL PARAMETERS FOR THE PEC
      READ(9,'(A100)') LINE
      WRITE(6,'(" VL8_INIT: ",A100)') LINE
      READ(9,*) SCL,NPT,NPTMX,ESHIFT
      WRITE(6,'(13X,F8.5,2(2X,I4),2X,F10.6)') SCL,NPT,NPTMX,ESHIFT
C . INPUT PARAMETERS FOR RKHS INTER-/EXTRAPOLATION 
      READ(9,'(A100)') LINE
      WRITE(6,'(" VL8_INIT: ",A100)') LINE
      READ(9,*) N,M,S,ASYMP,NCSN,CSN_PRINT
      WRITE(6,'(13X,3(I2,2X),L2,2X,I3,2X,L2)')N,M,S,ASYMP,NCSN,CSN_PRINT
      IF (ASYMP) THEN
         ALLOCATE (CSN(NCSN))
         ALLOCATE (RR(NPT+NCSN))
         READ(9,'(A100)') LINE
         WRITE(6,'(" VL8_INIT: ",A100)') LINE
         READ(9,*) RA,CSN,RR(NPT+1:NPT+NCSN)
         WRITE(6,'(13X,F5.2,2(2X,F14.2),2(2X,F6.2))') RA,CSN,
     &                                                RR(NPT+1:NPT+NCSN)
      ELSE
         read(9,*)
         read(9,*)
         ALLOCATE (RR(NPTMX))
         NCSN=0
         RA=0.D0
      END IF
      ALLOCATE (VV(NPTMX))
      read(9,*)
      DO I=1,NPTMX
         READ(9,*) RR(I),DUMMY1,DUMMY2,VV(I)
         RR(I)=RR(I)*A02AA ! CONVERTS FROM a0 TO AA
!         write(28,*) rr(i),vv(i)
      END DO
      CLOSE(9)
C COMPUTES RKHS COEFFICIENTS ANMK AND ALPHA
      NN=NPT+NCSN
      ALLOCATE (ANMK(N))
      ALLOCATE (ALPHA(NN))
      CALL RKHS_INIT(RR,VV,NPT,N,M,S,ASYMP,CSN,NCSN,RA,ANMK,ALPHA,
     & CSN_PRINT)
C
      IF (ASYMP) DEALLOCATE (CSN)
C
      RETURN
C------------------------------- VL8_EVAL ------------------------------
      ENTRY VL8_EVAL(X,VL8)
C-----------------------------------------------------------------------
      VL8=V_RKHS(X,RR,NN,N,M,S,ANMK,ALPHA)
C
      RETURN
 2018 PRINT*, " ERROR IN VL8_INIT: COULD NOT OPEN FILE ", FILE_IN
      STOP
      END SUBROUTINE VL8_INIT
C******************************* VH8_INIT ******************************
      SUBROUTINE VH8_INIT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C-----------------------------------------------------------------------
C  DATE (LAST UPDATE): 07/07/14       STATUS: FINISHED, UNDER TEST     |
C  AUTHOR: MAYKEL LEONARDO GONZALEZ MARTINEZ                           |
C-----------------------------------------------------------------------
C  THIS SUBROUTINE EVALUATES THE |
C-----------------------------------------------------------------------
C USES VMORSE, RKHS_INIT, V_RKHS
C PARAMETERS
      PARAMETER (A02AA=0.52917726D0)
C 'EXTRA' VARIABLES:
      INTEGER S
      CHARACTER FILE_IN*26,LINE*100
      LOGICAL ASYMP,CSN_PRINT
      ALLOCATABLE RR(:),VV(:),CSN(:),ANMK(:),ALPHA(:)
C
C OPENS, READS & CLOSES INPUT/PARAMETERS FILE
      FILE_IN='V4-10.LiEr.MCSCF.S_1.5.dat'
      OPEN(9,FILE=FILE_IN,STATUS='OLD',ERR=2019)
      WRITE(6,'(/" V^S=3/2_k=8 FROM FILE: ",A30)') FILE_IN
C . GLOBAL PARAMETERS FOR THE PEC
      READ(9,'(A100)') LINE
      WRITE(6,'(" VH8_INIT: ",A100)') LINE
      READ(9,*) SCL,NPT,NPTMX,ESHIFT
      WRITE(6,'(13X,F8.5,2(2X,I4),2X,F10.6)') SCL,NPT,NPTMX,ESHIFT
C . INPUT PARAMETERS FOR RKHS INTER-/EXTRAPOLATION 
      READ(9,'(A100)') LINE
      WRITE(6,'(" VH8_INIT: ",A100)') LINE
      READ(9,*) N,M,S,ASYMP,NCSN,CSN_PRINT
      WRITE(6,'(13X,3(I2,2X),L2,2X,I3,2X,L2)')N,M,S,ASYMP,NCSN,CSN_PRINT
      IF (ASYMP) THEN
         ALLOCATE (CSN(NCSN))
         ALLOCATE (RR(NPT+NCSN))
         READ(9,'(A100)') LINE
         WRITE(6,'(" VH8_INIT: ",A100)') LINE
         READ(9,*) RA,CSN,RR(NPT+1:NPT+NCSN)
         WRITE(6,'(13X,F5.2,2(2X,F14.2),2(2X,F6.2))') RA,CSN,
     &                                                RR(NPT+1:NPT+NCSN)
      ELSE
         read(9,*)
         read(9,*)
         ALLOCATE (RR(NPTMX))
         NCSN=0
         RA=0.D0
      END IF
      ALLOCATE (VV(NPTMX))
      read(9,*)
      DO I=1,NPTMX
         READ(9,*) RR(I),DUMMY1,DUMMY2,VV(I)
         RR(I)=RR(I)*A02AA ! CONVERTS FROM a0 TO AA
!         write(29,*) rr(i),vv(i)
      END DO
      CLOSE(9)
C COMPUTES RKHS COEFFICIENTS ANMK AND ALPHA
      NN=NPT+NCSN
      ALLOCATE (ANMK(N))
      ALLOCATE (ALPHA(NN))
      CALL RKHS_INIT(RR,VV,NPT,N,M,S,ASYMP,CSN,NCSN,RA,ANMK,ALPHA,
     & CSN_PRINT)
C
      IF (ASYMP) DEALLOCATE (CSN)
C
      RETURN
C------------------------------- VH8_EVAL ------------------------------
      ENTRY VH8_EVAL(X,VH8)
C-----------------------------------------------------------------------
      VH8=V_RKHS(X,RR,NN,N,M,S,ANMK,ALPHA)
C
      RETURN
 2019 PRINT*, " ERROR IN VH8_INIT: COULD NOT OPEN FILE ", FILE_IN
      STOP
      END SUBROUTINE VH8_INIT
C******************************* VL10_INIT ******************************
      SUBROUTINE VL10_INIT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C-----------------------------------------------------------------------
C  DATE (LAST UPDATE): 07/07/14       STATUS: FINISHED, UNDER TEST     |
C  AUTHOR: MAYKEL LEONARDO GONZALEZ MARTINEZ                           |
C-----------------------------------------------------------------------
C  THIS SUBROUTINE EVALUATES THE |
C-----------------------------------------------------------------------
C USES VMORSE, RKHS_INIT, V_RKHS
C PARAMETERS
      PARAMETER (A02AA=0.52917726D0)
C 'EXTRA' VARIABLES:
      INTEGER S
      CHARACTER FILE_IN*26,LINE*100
      LOGICAL ASYMP,CSN_PRINT
      ALLOCATABLE RR(:),VV(:),CSN(:),ANMK(:),ALPHA(:)
C
C OPENS, READS & CLOSES INPUT/PARAMETERS FILE
      FILE_IN='V4-10.LiEr.MCSCF.S_0.5.dat'
      OPEN(9,FILE=FILE_IN,STATUS='OLD',ERR=2020)
      WRITE(6,'(/" V^S=1/2_k=10 FROM FILE: ",A30)') FILE_IN
C . GLOBAL PARAMETERS FOR THE PEC
      READ(9,'(A100)') LINE
      WRITE(6,'(" VL10_INIT: ",A100)') LINE
      READ(9,*) SCL,NPT,NPTMX,ESHIFT
      WRITE(6,'(13X,F8.5,2(2X,I4),2X,F10.6)') SCL,NPT,NPTMX,ESHIFT
C . INPUT PARAMETERS FOR RKHS INTER-/EXTRAPOLATION 
      READ(9,'(A100)') LINE
      WRITE(6,'(" VL10_INIT: ",A100)') LINE
      READ(9,*) N,M,S,ASYMP,NCSN,CSN_PRINT
      WRITE(6,'(13X,3(I2,2X),L2,2X,I3,2X,L2)')N,M,S,ASYMP,NCSN,CSN_PRINT
      IF (ASYMP) THEN
         ALLOCATE (CSN(NCSN))
         ALLOCATE (RR(NPT+NCSN))
         READ(9,'(A100)') LINE
         WRITE(6,'(" VL10_INIT: ",A100)') LINE
         READ(9,*) RA,CSN,RR(NPT+1:NPT+NCSN)
         WRITE(6,'(13X,F5.2,2(2X,F14.2),2(2X,F6.2))') RA,CSN,
     &                                                RR(NPT+1:NPT+NCSN)
      ELSE
         read(9,*)
         read(9,*)
         ALLOCATE (RR(NPTMX))
         NCSN=0
         RA=0.D0
      END IF
      ALLOCATE (VV(NPTMX))
      read(9,*)
      DO I=1,NPTMX
         READ(9,*) RR(I),DUMMY1,DUMMY2,DUMMY3,VV(I)
         RR(I)=RR(I)*A02AA ! CONVERTS FROM a0 TO AA
!         write(30,*) rr(i),vv(i)
      END DO
      CLOSE(9)
C COMPUTES RKHS COEFFICIENTS ANMK AND ALPHA
      NN=NPT+NCSN
      ALLOCATE (ANMK(N))
      ALLOCATE (ALPHA(NN))
      CALL RKHS_INIT(RR,VV,NPT,N,M,S,ASYMP,CSN,NCSN,RA,ANMK,ALPHA,
     & CSN_PRINT)
C
      IF (ASYMP) DEALLOCATE (CSN)
C
      RETURN
C------------------------------- VL10_EVAL ------------------------------
      ENTRY VL10_EVAL(X,VL10)
C-----------------------------------------------------------------------
      VL10=V_RKHS(X,RR,NN,N,M,S,ANMK,ALPHA)
C
      RETURN
 2020 PRINT*, " ERROR IN VL10_INIT: COULD NOT OPEN FILE ", FILE_IN
      STOP
      END SUBROUTINE VL10_INIT
C******************************* VH10_INIT ******************************
      SUBROUTINE VH10_INIT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C-----------------------------------------------------------------------
C  DATE (LAST UPDATE): 07/07/14       STATUS: FINISHED, UNDER TEST     |
C  AUTHOR: MAYKEL LEONARDO GONZALEZ MARTINEZ                           |
C-----------------------------------------------------------------------
C  THIS SUBROUTINE EVALUATES THE |
C-----------------------------------------------------------------------
C USES VMORSE, RKHS_INIT, V_RKHS
C PARAMETERS
      PARAMETER (A02AA=0.52917726D0)
C 'EXTRA' VARIABLES:
      INTEGER S
      CHARACTER FILE_IN*26,LINE*100
      LOGICAL ASYMP,CSN_PRINT
      ALLOCATABLE RR(:),VV(:),CSN(:),ANMK(:),ALPHA(:)
C
C OPENS, READS & CLOSES INPUT/PARAMETERS FILE
      FILE_IN='V4-10.LiEr.MCSCF.S_1.5.dat'
      OPEN(9,FILE=FILE_IN,STATUS='OLD',ERR=2021)
      WRITE(6,'(/" V^S=3/2_k=10 FROM FILE: ",A30)') FILE_IN
C . GLOBAL PARAMETERS FOR THE PEC
      READ(9,'(A100)') LINE
      WRITE(6,'(" VH10_INIT: ",A100)') LINE
      READ(9,*) SCL,NPT,NPTMX,ESHIFT
      WRITE(6,'(13X,F8.5,2(2X,I4),2X,F10.6)') SCL,NPT,NPTMX,ESHIFT
C . INPUT PARAMETERS FOR RKHS INTER-/EXTRAPOLATION 
      READ(9,'(A100)') LINE
      WRITE(6,'(" VH10_INIT: ",A100)') LINE
      READ(9,*) N,M,S,ASYMP,NCSN,CSN_PRINT
      WRITE(6,'(13X,3(I2,2X),L2,2X,I3,2X,L2)')N,M,S,ASYMP,NCSN,CSN_PRINT
      IF (ASYMP) THEN
         ALLOCATE (CSN(NCSN))
         ALLOCATE (RR(NPT+NCSN))
         READ(9,'(A100)') LINE
         WRITE(6,'(" VH10_INIT: ",A100)') LINE
         READ(9,*) RA,CSN,RR(NPT+1:NPT+NCSN)
         WRITE(6,'(13X,F5.2,2(2X,F14.2),2(2X,F6.2))') RA,CSN,
     &                                                RR(NPT+1:NPT+NCSN)
      ELSE
         read(9,*)
         read(9,*)
         ALLOCATE (RR(NPTMX))
         NCSN=0
         RA=0.D0
      END IF
      ALLOCATE (VV(NPTMX))
      read(9,*)
      DO I=1,NPTMX
         READ(9,*) RR(I),DUMMY1,DUMMY2,DUMMY3,VV(I)
         RR(I)=RR(I)*A02AA ! CONVERTS FROM a0 TO AA
!         write(31,*) rr(i),vv(i)
      END DO
      CLOSE(9)
C COMPUTES RKHS COEFFICIENTS ANMK AND ALPHA
      NN=NPT+NCSN
      ALLOCATE (ANMK(N))
      ALLOCATE (ALPHA(NN))
      CALL RKHS_INIT(RR,VV,NPT,N,M,S,ASYMP,CSN,NCSN,RA,ANMK,ALPHA,
     & CSN_PRINT)
C
      IF (ASYMP) DEALLOCATE (CSN)
C
      RETURN
C------------------------------- VH10_EVAL ------------------------------
      ENTRY VH10_EVAL(X,VH10)
C-----------------------------------------------------------------------
      VH10=V_RKHS(X,RR,NN,N,M,S,ANMK,ALPHA)
C
      RETURN
 2021 PRINT*, " ERROR IN VH10_INIT: COULD NOT OPEN FILE ", FILE_IN
      STOP
      END SUBROUTINE VH10_INIT
C******************************** SWITCH *******************************
      FUNCTION SWITCH(BETA,X0,X)
      IMPLICIT NONE
C-----------------------------------------------------------------------
C  DATE (LAST UPDATE): 07/07/14       STATUS: FINISHED, UNDER TEST     |
C  AUTHOR: MAYKEL LEONARDO GONZALEZ MARTINEZ                           |
C-----------------------------------------------------------------------
C  THIS FUNCTION
C-----------------------------------------------------------------------
      REAL*8 SWITCH,BETA,X0,X
      SWITCH=1.D0/(1.D0+EXP(BETA*(X-X0)))
      RETURN
      END FUNCTION SWITCH
C*************************** RKHS_INIT *********************************
      SUBROUTINE RKHS_INIT(RR,VV,NPT,N,M,S,ASYMP,CSN,NCSN,RA,ANMK,ALPHA,
     & CSN_PRINT)
      IMPLICIT NONE
C-----------------------------------------------------------------------
C  DATE (LAST UPDATE): 20/06/12       STATUS: FINISHED, TESTED         |
C  AUTHOR: MAYKEL LEONARDO GONZALEZ MARTINEZ                           |
C-----------------------------------------------------------------------
C   THIS ROUTINE INITIALISES THE ANMK AND ALPHA ARRAYS USED BY V_RKHS  |
C            EQUATIONS FROM J. CHEM. PHYS. 113, 3960 (2000)            |
C-----------------------------------------------------------------------
C USES ANMK_EVAL (+ LAPACK'S DGESV & DCOPY)
      INTEGER NN,N,M,S,NCSN,NPT,I,J,K,IPIV,INFO
      REAL*8 RR(NPT+NCSN),VV(NPT),CSN(MAX(1,NCSN)),RA
      REAL*8 ANMK(N),ALPHA(NPT+NCSN),GAMMAS,BETA,RMIN,RMAX,CDISP
      LOGICAL ASYMP,CSN_PRINT
      ALLOCATABLE GAMMAS(:,:),BETA(:),IPIV(:)
C
      IF (S.LE.0) STOP " ERROR IN RKHS_INIT: S <= 0"
      IF (ASYMP .AND. NCSN.GT.N) STOP " ERROR IN RKHS_INIT: NCSN > N"
      NN=NPT+NCSN
      ALLOCATE (GAMMAS(NN,NN))
      ALLOCATE (BETA(NN))
      ALLOCATE (IPIV(NN))
C GETS A_{NMK=0,N) COEEFICIENTS USING EQ.(8)
      CALL ANMK_EVAL(N,M,ANMK)
      LI : DO I=1,NN ! LOOPS EVALUATE GAMMA^(S)_{IJ} FROM EQ.(9) & BETA_J FROM EQ.(10)
         LJ : DO J=1,NN
            IF (J.LE.NPT) THEN
               RMIN=MIN(RR(I),RR(J))
               RMAX=MAX(RR(I),RR(J))
               GAMMAS(J,I)=0.D0
               DO K=0,N-1
                  GAMMAS(J,I)=GAMMAS(J,I)+ANMK(K+1)*(RMIN/RMAX)**(S*K)
               END DO
               GAMMAS(J,I)=GAMMAS(J,I)*RMAX**(-S*(M+1)) ! Q^{N,M}_1(RI^S,RJ^S) FROM EQ.(1)
            ELSE
               K=J-NPT-1
               GAMMAS(J,I)=ANMK(K+1)*RR(I)**(S*K)/RA**(S*(K+M+1))
            END IF
         END DO LJ
         IF (I.LE.NPT) THEN
            BETA(I)=VV(I)
         ELSE
            K=I-NPT-1
            BETA(I)=-CSN(K+1)*RA**(-S*(K+M+1))
         END IF
      END DO LI
      CALL DGESV(NN,1,GAMMAS,NN,IPIV,BETA,NN,INFO) ! SOLVES EQ.(8) FOR ALPHA
      CALL DCOPY(NN,BETA,1,ALPHA,1) ! COPIES SOLUTION INTO ALPHA
      DEALLOCATE (GAMMAS,BETA,IPIV)
c (IF REQUESTED = CSN_PRINT) PRINTS OUT THE N DISPERSION COEFFICIENTS
      IF (CSN_PRINT) THEN
         WRITE(6,*) "RKHS_INIT: DISPERSION COEFFICIENTS"
         DO K=0,N-1
            CDISP=0.D0
            DO I=1,NN
               CDISP=CDISP+ALPHA(I)*RR(I)**(S*K)
            END DO
            CDISP=-CDISP*ANMK(K+1) ! C_{S(K+M+1)} FROM EQ.(6)
            IF (K+1.LE.NCSN) THEN
               WRITE(6,1) S*(K+M+1),CDISP,CSN(K+1),
     &                    ABS((CDISP-CSN(K+1))/(CDISP+CSN(K+1)))*200
            ELSE
               WRITE(6,2) S*(K+M+1),CDISP
            END IF
         END DO
      END IF
C
      RETURN
 1    FORMAT(2X,"C",I2," = ",F18.6," (out) ",F18.6," (in)",
     &       E15.2," (% error)")
 2    FORMAT(2X,"C",I2," = ",F18.6," (out) ")
      END SUBROUTINE RKHS_INIT
C**************************** V_RKHS ***********************************
      FUNCTION V_RKHS(R,RR,NN,N,M,S,ANMK,ALPHA)
      IMPLICIT NONE
C-----------------------------------------------------------------------
C  DATE (LAST UPDATE): 19/06/12       STATUS: FINISHED, TESTED         |
C  AUTHOR: MAYKEL LEONARDO GONZALEZ MARTINEZ                           |
C-----------------------------------------------------------------------
C     THIS FUNCTION INTERPOLATES/EXTRAPOLATES USING THE RKHS METHOD    |
C            EQUATIONS FROM J. CHEM. PHYS. 113, 3960 (2000)            |
C-----------------------------------------------------------------------
      INTEGER NN,N,M,S,I,K
      REAL*8 R,RR(NN),ANMK(N),ALPHA(NN),V_RKHS,RMIN,RMAX,QNM1
C
      V_RKHS=0.D0
      DO I=1,NN ! LOOPS EVALUATE Q^{N,M}_1(RI^S,R^S) FROM EQ.(1)
         RMIN=MIN(R,RR(I))
         RMAX=MAX(R,RR(I))
         QNM1=0.D0
         DO K=0,N-1
            QNM1=QNM1+ANMK(K+1)*(RMIN/RMAX)**(S*K)
         END DO
         QNM1=QNM1*RMAX**(-S*(M+1))
         V_RKHS=V_RKHS+QNM1*ALPHA(I)! EVALUATES V USING EQ.(4)
      END DO
C
      RETURN
      END FUNCTION V_RKHS
C****************************** ANMK_EVAL ******************************
      SUBROUTINE ANMK_EVAL(N,M,ANMK)
      IMPLICIT NONE
C-----------------------------------------------------------------------
C  DATE (LAST UPDATE): 19/06/12       STATUS: FINISHED, TESTED         |
C  AUTHOR: MAYKEL LEONARDO GONZALEZ MARTINEZ                           |
C-----------------------------------------------------------------------
C  THIS ROUTINE EVALUATES THE (ARRAY) OF ANMK_{N M 0:N-1) COEFFICIENTS |
C           USING EQ.(7) IN J. CHEM. PHYS. 113, 3960 (2000)            |
C [NOTICE: A_{NMK} == BETA^{N,M}_K IN J. CHEM. PHYS. 112, 4415 (2000)] |
C-----------------------------------------------------------------------
C USES FACTORIAL, POCHHAMMER (+ BETA)
      INTEGER I,J,N,M,K
      INTEGER*8 FACTORIAL,POCHHAMMER
      REAL*8 BETA,ANMK(N),FCT
C STATEMENT FUNCTIONS
C  BETA FUNCTION WITH INTEGER ARGS: BETA(I,J) = (I-1)!(J-1)!/(I+J-1)!
      BETA(I,J)=FACTORIAL(I-1)*FACTORIAL(J-1)/DFLOAT(FACTORIAL(I+J-1))
C
      FCT=N*N*BETA(N,M+1)
      DO K=0,N-1
         ANMK(K+1)=FCT*POCHHAMMER(-N+1,K)*POCHHAMMER(M+1,K)
     &             /DFLOAT(POCHHAMMER(N+M+1,K)*FACTORIAL(K))
      END DO
C
      RETURN
      END SUBROUTINE ANMK_EVAL
C****************************** FACTORIAL ******************************
      FUNCTION FACTORIAL(N)
      IMPLICIT NONE
C-----------------------------------------------------------------------
C  DATE (LAST UPDATE): 19/06/12       STATUS: FINISHED, TESTED         |
C  AUTHOR: MAYKEL LEONARDO GONZALEZ MARTINEZ                           |
C-----------------------------------------------------------------------
C             THIS FUNCTION EVALUATES THE FACTORIAL OF N               |
C    (IT USES A VERY SIMPLE ALGORITHM, IT SHOULD BE OK FOR SMALL N)    |
C-----------------------------------------------------------------------
      INTEGER N,FCT10(0:10),I
      INTEGER*8 FACTORIAL
      DATA FCT10/1,1,2,6,24,120,720,5040,40320,362880,3628800/
C
      SELECT CASE (N)
      CASE (:-1)
         STOP " ERROR IN FACTORIAL: N < 0"
      CASE (0:10)
         FACTORIAL=FCT10(N)
      CASE DEFAULT
         FACTORIAL=FCT10(10)
         DO I=11,N
            FACTORIAL=FACTORIAL*I
         END DO
      END SELECT
C
      RETURN
      END FUNCTION FACTORIAL
C****************************** POCHAMMER ******************************
      FUNCTION POCHHAMMER(I,N)
      IMPLICIT NONE
C-----------------------------------------------------------------------
C  DATE (LAST UPDATE): 19/06/12       STATUS: FINISHED, TESTED         |
C  AUTHOR: MAYKEL LEONARDO GONZALEZ MARTINEZ                           |
C-----------------------------------------------------------------------
C    THIS FUNCTION EVALUATES THE POCHHAMMER'S SYMBOL (I, N INTEGERS)   |
C            (I)N = {1 IF N = 0; I(I+1)...(I+N-1) IF N > 0             |
C  (IT USES A VERY SIMPLE ALGORITHM, IT SHOULD BE OK FOR SMALL I & N)  |
C-----------------------------------------------------------------------
      INTEGER I,J,N
      INTEGER*8 POCHHAMMER
C
      SELECT CASE (N)
      CASE (:-1)
         STOP " ERROR IN POCHHAMMER: N < 0"
      CASE (0)
         POCHHAMMER=1
      CASE DEFAULT
         POCHHAMMER=I
         DO J=1,N-1
            POCHHAMMER=POCHHAMMER*(I+J)
         END DO
      END SELECT
C
      RETURN
      END FUNCTION POCHHAMMER
