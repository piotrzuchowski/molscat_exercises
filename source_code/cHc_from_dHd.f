      FUNCTION cHc_from_dHd(XJ1C,XJ2C,XJC,XMC,XJ1R,XJ2R,XJR,XMR,dHd,IX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C  DATE (LAST UPDATE): 14/05/13        STATUS: FINISHED, TESTED        |
C  AUTHOR: MAYKEL LEONARDO GONZALEZ MARTINEZ                           |
C-----------------------------------------------------------------------
C    THIS FUNCTION COMPUTES A MATRIX ELEMENT IN A "COUPLED" BASIS      |
C                cHc = <(j1 j2)j mj|H|(j1' j2')j' mj'>                 |
C    USING THE MATRIX ELEMENTS IN THE ASSOCIATED "DECOUPLED" BASIS     |
C                dHd = <(j1 j2)j mj|H|(j1' j2')j' mj'>                 |
C-----------------------------------------------------------------------
C USES THRJ, PARITY (+ dHd, EXTERNAL)
C 'EXTRA' VARIABLES: INTEGER ARRAY TO PASS dHd OTHER INDEXES IT MAY NEED
      DIMENSION IX(*)
C GETS REAL QUANTUM NUMBERS
C . "COLUMN" QUANTUM NUMBERS - BRA
      J1C=2*XJ1C
      J2C=2*XJ2C
      JC =2*XJC
      MC =2*XMC
C . "ROW" QUANTUM NUMBERS - KET
      J1R=2*XJ1R
      J2R=2*XJ2R
      JR =2*XJR
      MR =2*XMR
C COMPUTES cHc FROM dHd AND COUPLING COEFFICIENTS
      cHc_from_dHd=0.D0
      CLOOP: DO M1C=-J1C,J1C,2
         M2C=MC-M1C
         IF (ABS(M2C).GT.J2C) CYCLE CLOOP
         XM1C=0.5D0*M1C
         XM2C=0.5D0*M2C
         RLOOP: DO M1R=-J1R,J1R,2
            M2R=MR-M1R
            IF (ABS(M2R).GT.J2R) CYCLE RLOOP
            XM1R=0.5D0*M1R
            XM2R=0.5D0*M2R
            cHc_from_dHd=cHc_from_dHd
     &       +THRJ(XJ1C,XJ2C,XJC, XM1C,XM2C,-XMC)
     &       *THRJ(XJ1R,XJ2R,XJR, XM1R,XM2R,-XMR)
     &       *dHd(XJ1C,XM1C,XJ2C,XM2C,XJ1R,XM1R,XJ2R,XM2R,IX)
            IF (ISNAN(cHc_from_dHd)) STOP "  ERROR: isNaN(cHc_from_dHd)"
         END DO RLOOP
      END DO CLOOP
      cHc_from_dHd=cHc_from_dHd
     &             *PARSGN((J1C-J2C+MC+J1R-J2R+MR)/2)
     &             *SQRT((JC+1.D0)*(JR+1.D0))
C
      RETURN
      END FUNCTION cHc_from_dHd
