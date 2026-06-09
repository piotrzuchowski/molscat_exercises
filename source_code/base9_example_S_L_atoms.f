      SUBROUTINE BAS9IN(PRTP,IBOUND)
      USE efvs
      USE potential
      USE basis_data
      USE physical_constants
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C-----------------------------------------------------------------------
C  DATE (LAST UPDATE): 22/08/14        STATUS: FINISHED, TESTED        |
C  AUTHOR: MAYKEL LEONARDO GONZALEZ MARTINEZ                           |
C-----------------------------------------------------------------------
C    THIS SUBROUTINE & ITS ENTRIES SET THE BASIS SET (I.E. QUANTUM     |
C    NUMBERS, COUPLING MATRIX ELEMENTS...) FOR THE COLLISION TYPES:    |
C                  ATOM(S) + ATOM(L) IN MAGNETIC FIELDS                |
C                          USING THE BASIS SET                         |
C          |i1 ms1> |s1 ms1> |i2 mi2> |(l2 s2)j2 mj2> |L ML>           |
C   HYPERFINE TERMS CAN BE INCLUDED/OMITTED WITH HFFLG = TRUE/FALSE    |
C    * ALL MONOMER QUANTUM NUMBERS & MTOT ARE TWICE THEIR VALUES *     |
C-----------------------------------------------------------------------
C PARAMETERS

      PARAMETER (L2MAX=7)
C EXTERNAL
      DIMENSION IXS1L(7),IXSK(6)
C 'EXTRA' VARIABLES:
      CHARACTER*8 PRTP(4),QNAME(10),CL2(L2MAX)
      LOGICAL LEVIN,EIN,LCNT
      DIMENSION ROTI(12),ELEVEL(1000),JLEVEL(4000),
     1   ISYM(10),ISYM2(10),WT(2)
      DIMENSION JLEV(1),VL(1),IV(1),J(1),L(1),CENT(1),LAM(1)
C
      CHARACTER A2SP1,B2SP1
      LOGICAL ODD,DIAG2,DIAG3,HFFLG,M2S1,D1C1R,D2C2R,DHF12
      DIMENSION AJ(3,3),BJ(3,3)
      DIMENSION IS(10)
C COMMONS
      COMMON /CMBASE/ ROTI,ELEVEL,EMAX,WT,SPNUC,NLEVEL,JLEVEL,
     1   IDUM1,IDUM2,JSTEP,ISYM,J2MIN,J2MAX,J2STEP,ISYM2,JHALF,IDENT,
     2   MXJL,MXEL
      COMMON /VLFLAG/ IVLFL
      COMMON /VLUSE/ NCONST,NRSQ,VCONST(VLCSIZ),IIREF,
     1               IFIELD,NNCONST,MCONST(VLCSIZ)


