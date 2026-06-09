
module file_dir
 implicit none
 integer, parameter :: len_ch=300
 character(len=*), parameter :: P_pattern='P_pattern.data', dispinp_pattern='pattern.dispinp'
 character(len=*), parameter :: file_A1_tmp='atomA1.tmp', file_B1_tmp='atomB1.tmp'
 character(len=*), parameter :: file_ele='Electrostatic.dat'
 character(len=*), parameter :: file_ind='Induction.dat'
 character(len=*), parameter :: file_dis='Dispersion.dat'
 character(len=*), parameter :: file_geoparm='geoparm.d'
 character(len=*), parameter :: file_dimer='dimer.cnf'
 character(len=*), parameter :: file_part_charges_A='Part_Charges_A.tmp', file_part_charges_B='Part_Charges_B.tmp'
 character(len=len_ch) :: P_file, dispinp_file
end module file_dir

module double
 implicit none
 integer, parameter :: dbl = kind(1.0d0)
end module double


!--------------------------------
!    MODULE : CONST
!--------------------------------
module const
  use double
  implicit none
  
  integer, parameter :: quadp = selected_real_kind(33, 4931)
  integer, parameter :: dblep = selected_real_kind(15, 307)
  integer, parameter :: STDOUT = 6
  integer, parameter :: MAX_ANUM = 36, MAX_NSITE=100, MAX_IDOF = 10, MAX_L = 6
  
  ! Energy component indices
  integer, parameter :: n_ener_inter = 8, n_ener = n_ener_inter+2
  integer, parameter :: E_tot = 1, E_els = 2, E_exh = 3, E_ind = 4, E_exhind = 5, E_dsp = 6, E_exhdsp = 7
  integer, parameter :: E_dhf = 8, E_mon(2) = [9,10]
  
  ! PES component indices
  integer, parameter :: n_FF_comp = 4
  integer, parameter :: F_els = 1, F_pol = 2, F_vdw = 3, F_exp = 4
  
  ! Combination rule indices (geometric, arithmetic)
  integer, parameter :: CR_geo = 1, CR_arith = 2
  
  ! Parameter class indices (global, site-specific, or site-pair-specific)
  integer, parameter :: PC_glb = 2, PC_s = 1, PC_sp = 0
  
  ! IDOF type indices (angular, linear, chiral)
  integer, parameter :: MC_ANG = 1, MC_LIN = 2, MC_CHI = 3
  
  real(kind=dbl), parameter :: PI = 2.0d0*dacos(0.0d0)
  real(kind=dbl), parameter :: twopi = 2.0d0*PI
  real(kind=dbl), parameter :: rad2d = 180.0d0/PI
  real(kind=dbl), parameter :: d2rad = PI/180.0d0
  real(kind=dbl), parameter :: a0 = 0.529177249d0
  real(kind=dbl), parameter :: bohr2A = a0
  real(kind=dbl), parameter :: A2bohr = 1/bohr2A
  real(kind=dbl), parameter :: J2cal = 0.239005736d0
  real(kind=dbl), parameter :: h2kcal = 627.510d0
  real(kind=dbl), parameter :: kcal2h = 1/h2kcal
  real(kind=dbl), parameter :: h2ev = 27.211396132d0
  real(kind=dbl), parameter :: atm2au = 3.443967d-9
  character(1), parameter :: orbSymbol(0:MAX_L) = ['s','p','d','f','g','h','i']
  real(8), parameter :: VDWRad(MAX_ANUM) = A2bohr*[ &
       1.20,                                                                                1.40, &
       1.82,0.00,                                                  0.00,1.70,1.55,1.52,1.47,1.54, &
       2.27,1.73,                                                  0.00,2.10,1.80,1.80,1.75,1.88, &
       2.65,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,1.63,1.40,1.39,1.87,0.00,1.85,1.90,1.85,2.02]
  real(8), parameter :: covRad(MAX_ANUM) = A2bohr*[ &
       0.37,                                                                                0.32, &
       1.34,0.90,                                                  0.82,0.77,0.75,0.73,0.71,0.69, &
       1.54,1.30,                                                  1.18,1.11,1.06,1.02,0.99,0.97, &
       1.96,1.74,1.44,1.36,1.25,1.27,1.39,1.25,1.26,1.21,1.38,1.31,1.26,1.22,1.19,1.16,1.14,1.10]
  real(8), parameter :: atomMass(0:MAX_ANUM) = [ &
  0.0000, &
  1.0079,                                                                                                                4.0026, &
  6.9410,9.0122,                                                                      10.811,12.011,14.001,15.999,18.998,20.180, &
  22.990,24.305,                                                                      26.981,28.085,30.974,32.065,35.453,39.948, &
  39.098,40.078,44.956,47.867,50.941,51.996,54.938,55.845,58.933,58.693,63.546,65.382,69.723,72.640,74.922,78.971,79.904,83.798]
  real(8), parameter :: electroNeg(MAX_ANUM) = [ &
       2.20,                                                                                0.00, &
       0.98,1.57,                                                  2.04,2.55,3.04,3.44,3.98,0.00, &
       0.93,1.31,                                                  1.61,1.90,2.19,2.58,3.16,0.00, &
       0.82,1.00,1.36,1.54,1.63,1.66,1.55,1.83,1.88,1.91,1.90,1.65,1.81,2.01,2.18,2.55,2.96,3.00]
  real(8), parameter :: atomPol(MAX_ANUM) = [ &
       4.50,                                                                                1.38, &
       164.,37.8,                                                  20.5,11.0,7.43,6.04,3.76,2.67, &
       162.,73.4,                                                  46.0,36.5,24.9,19.4,14.6,11.1, &
       293.,169.,130.,100.,85.0,80.0,65.0,60.0,55.0,50.0,45.0,38.8,50.0,40.0,29.5,26.2,21.9,17.1]
  character(2), parameter :: asymbol(0:MAX_ANUM) = [ &
       'Mb', &
       'H ',                                                                                'He', &
       'Li','Be',                                                  'B ','C ','N ','O ','F ','Ne', &
       'Na','Mg',                                                  'Al','Si','P ','S ','Cl','Ar', &
       'K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr']
  integer, parameter :: maxOccL(MAX_ANUM) = [ &
       0,                                0, &
       0,0,                    1,1,1,1,1,1, &
       1,1,                    1,1,1,1,1,1, &
       1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2]
end module const


!--------------------------------
!    MODULE : UTILITY
!--------------------------------
module utility
  use const
  implicit none
  
  contains
    
  ! Counts the number of lines in a file, returns 0 if file does not exist
  function nLines(filename)
    implicit none
    character(*) :: filename
    integer :: nLines
    
    integer, parameter :: nLinesUnit = 100
    
    nLines = 0
    open(nLinesUnit, file = filename, status = 'old', ERR=10)
    
    do
       read(nLinesUnit, *, END=10)
       nLines = nLines + 1
    end do
    
10  close(nLinesUnit)
  end function nLines
  
  
  ! Parses a string to an array of reals. All array elements after end of valid reals are set to zero.
  !    Ends after the first string which is not a valid real.
  subroutine splitLineReal(line, arr)
    implicit none
    character(*) :: line
    real(8) :: arr(:)
    
    integer :: i
    integer :: start, arrIndex
    logical :: blank, blankChar
    
    arr(:) = 0.0d0
    
    blank = .true.
    arrIndex = 1
    do i = 1,len(line)
       blankChar = index('0123456789.+-EeDd', line(i:i)) == 0
       
       if (blank) then
          if (.not. blankChar) then
             start = i
             blank = .false.
          endif
       else
          if (blankChar .or. i == len(line)) then
             read(line(start:i-1),*,err=10) arr(arrIndex)
             arrIndex = arrIndex + 1
             blank = .true.
          endif
       endif
    enddo
10  continue
  end subroutine splitLineReal
  
  function to_upper(strIn) result(strOut)
    ! Adapted from http://www.star.le.ac.uk/~cgp/fortran.html (25 May 2012) Original author: Clive Page    
    implicit none
    character(len=*), intent(in) :: strIn
    character(len=len(strIn)) :: strOut
    integer :: i,j
    
    do i = 1, len(strIn)
       j = iachar(strIn(i:i))
       if (j>= iachar("a") .and. j<=iachar("z") ) then
          strOut(i:i) = achar(iachar(strIn(i:i))-32)
       else
          strOut(i:i) = strIn(i:i)
       end if
    end do
  end function to_upper
  
  ! Convert string to atomic number
  function getAnum(str)
    integer :: getAnum
    character(*) :: str
    integer :: i
    
    do i = 1,MAX_ANUM
       if (trim(str) == trim(asymbol(i))) then
          getAnum = i
          return
       endif
    enddo
    
    read(str,*,err=10) getAnum
    return
    
10  continue
    getAnum = 0
  end function getAnum
  
  
  ! Gives set of 5 euler angles equivalent to the given 6 angles,
  !    but in the ranges ([0,pi],[0,2pi],[0,2pi],[0,pi],[0,2pi])
  subroutine fixDimerAngles(alphaA0, betaA, gammaA, alphaB, betaB, gammaB)
    real(8), intent(in) :: alphaA0
    real(8), intent(inout) :: betaA, gammaA, alphaB, betaB, gammaB
    real(8) :: alphaA
    
    alphaA = alphaA0
    call fixAngles(alphaA, betaA, gammaA)
    call fixAngles(alphaB, betaB, gammaB)
    alphaB = modulo(alphaB - alphaA, twopi)
  end subroutine fixDimerAngles
  
  ! Gives set of euler angles equivalent to the given one,
  !    but in the ranges ([0,2pi],[0,pi],[0,2pi])
  subroutine fixAngles(alpha, beta, gamma)
    real(8), intent(inout) :: alpha, beta, gamma
    
    beta = modulo(beta, twopi)
    if (beta > pi) then
       beta = twopi - beta
       alpha = alpha + pi
       gamma = gamma + pi
    endif
    alpha = modulo(alpha, twopi)
    gamma = modulo(gamma, twopi)
  end subroutine fixAngles
  
  
  ! Difference between two angles between 0 and 2pi, accounting for periodicity
  function angleDiff(angle1, angle2)
    real(8) :: angle1, angle2
    real(8) :: angleDiff
    
    angleDiff = angle2 - angle1
    if (angleDiff > twopi/2) angleDiff = angleDiff - twopi
    if (angleDiff < -1.0d0*twopi/2) angleDiff = angleDiff + twopi
  end function angleDiff
  
  
  ! Construct the rotation matrix coresponding to
  ! the ZYZ Euler angles ad,bd,gd (in radians).
  subroutine rotmat(ad,bd,gd,u)
    real(8), intent(in) :: ad, bd, gd
    real(8), intent(out) :: u(3,3)
    real(8) :: sad, sbd, sgd, cad, cbd, cgd
    
    sad = dsin(ad)
    sbd = dsin(bd)
    sgd = dsin(gd)
    cad = dcos(ad)
    cbd = dcos(bd)
    cgd = dcos(gd)
    u(1,1) = cad*cbd*cgd - sad*sgd
    u(1,2) = -cad*cbd*sgd - sad*cgd
    u(1,3) = cad*sbd
    u(2,1) = sad*cbd*cgd + cad*sgd
    u(2,2) = -sad*cbd*sgd + cad*cgd
    u(2,3) = sad*sbd
    u(3,1) = -sbd*cgd
    u(3,2) = sbd*sgd
    u(3,3) = cbd
    return
  end subroutine rotmat
  
  
  ! Construct matrix for rotation of a point about the line passing through p
  ! in the direction of unit vector q by an angle theta. Result is in homogeneous coordinates.
  subroutine rotmat2(p, q, theta, u)
    real(8), intent(out) :: u(4,4)
    real(8), intent(in) :: p(3), q(3), theta
    real(8) :: ct, st
    
    ct = cos(theta)
    st = sin(theta)
    u(1,1) = q(1)**2 + (q(2)**2 + q(3)**2)*ct
    u(1,2) = q(1)*q(2)*(1.0d0 - ct) - q(3)*st
    u(1,3) = q(1)*q(3)*(1.0d0 - ct) + q(2)*st
    u(1,4) = (p(1)*(q(2)**2 + q(3)**2) - q(1)*(p(2)*q(2) + p(3)*q(3)))*(1.0d0 - ct) + (p(2)*q(3) - p(3)*q(2))*st
    u(2,1) = q(1)*q(2)*(1.0d0 - ct) + q(3)*st
    u(2,2) = q(2)**2 + (q(1)**2 + q(3)**2)*ct
    u(2,3) = q(2)*q(3)*(1.0d0 - ct) - q(1)*st
    u(2,4) = (p(2)*(q(1)**2 + q(3)**2) - q(2)*(p(1)*q(1) + p(3)*q(3)))*(1.0d0 - ct) + (p(3)*q(1) - p(1)*q(3))*st
    u(3,1) = q(1)*q(3)*(1.0d0 - ct) - q(2)*st
    u(3,2) = q(2)*q(3)*(1.0d0 - ct) + q(1)*st
    u(3,3) = q(3)**2 + (q(1)**2 + q(2)**2)*ct
    u(3,4) = (p(3)*(q(1)**2 + q(2)**2) - q(3)*(p(1)*q(1) + p(2)*q(2)))*(1.0d0 - ct) + (p(1)*q(2) - p(2)*q(1))*st
    u(4,1) = 0.0d0
    u(4,2) = 0.0d0
    u(4,3) = 0.0d0
    u(4,4) = 1.0d0
  end subroutine rotmat2
  
  
  ! Parameter combination rules
  real(8) function combine(val1, val2, rule)
    real(8), intent(in) :: val1, val2
    integer, intent(in) :: rule
    
    combine = 0.0d0
    select case(rule)
    case(CR_geo)
       combine = sqrt(val1*val2)
    case(CR_arith)
       combine = 0.5d0*(val1 + val2)
    end select
  end function combine
  
  
  ! Sigmoid decode with min, max
  real(8) function decodeSig(x, min, max)
    real(8), intent(in) :: x, min, max
    real(8) :: avg, diff, const
    
    avg = (max + min)/2.0d0
    diff = max - min
    
    decodeSig = (diff/2.0d0)*(x/(1.0d0 + abs(x))) + avg
  end function decodeSig
  
  ! Sigmoid encode with min, max
  real(8) function encodeSig(f, min, max)
    real(8), intent(in) :: f, min, max
    real(8) :: f1, avg, diff, const
    
    if (f >= max) then
       f1 = 0.999d0*max + 0.001d0*min
    elseif (f <= min) then
       f1 = 0.999d0*min + 0.001d0*max
    else
       f1 = f
    endif
    
    avg = (max + min)/2.0d0
    diff = max - min
    const = (f1 - avg)*2.0d0/diff
    
    encodeSig = const/(1.0d0 - abs(const))
  end function encodeSig
  
  
  ! Encode parameter for fitting procedure
  real(8) function encodeLog(x)
    real(8), intent(in) :: x
    
    encodeLog = log(max(x,1.0d-8))
  end function encodeLog
  
  ! Decode parameter for fitting procedure
  real(8) function decodeLog(x)
    real(8), intent(in) :: x
    
    decodeLog = exp(x)
  end function decodeLog
  
  
  ! Encode parameter for fitting procedure
  real(8) function encodePow(x, pow)
    real(8), intent(in) :: x, pow
    
    if (x >= 0.0d0) then
       encodePow = x**(1.0d0/pow)
    else
       encodePow = -1.0d0*((-1.0d0*x)**(1.0d0/pow))
    end if
  end function encodePow
  
  ! Decode parameter for fitting procedure
  real(8) function decodePow(x, pow)
    real(8), intent(in) :: x, pow
    
    if (x >= 0.0d0) then
       decodePow = x**pow
    else
       decodePow = -1.0d0*((-1.0d0*x)**pow)
    end if
  end function decodePow
  
  ! Cartesian norm of 3-vector
  real(8) function cartnorm(vect)
    real(8), intent(in) :: vect(3)
    
    cartnorm = sqrt(vect(1)**2 + vect(2)**2 + vect(3)**2) 
  end function cartnorm
  
  ! Tang-Toennies damping funciton
  real(8) function TT(n, delta, r)
    integer, intent(in) :: n
    real(8), intent(in) :: delta, r
    
    real(8) :: rScaled, sum, fact, rtoi
    integer :: i
    
    rScaled = r*delta
    rtoi = 1.0
    fact = 1.0
    sum = 1.0
    do i = 1,n
       fact = fact*i
       rtoi = rtoi*rScaled
       sum = sum + rtoi/fact
    end do
    
    TT = 1.0 - exp(-rScaled)*sum
  end function TT

end module utility


!--------------------------------
!    MODULE : WEIGHT
!--------------------------------
module weight
  implicit none

contains
  
  real(8) function gweight(e, e0, weightScale)
    real(8), intent(in) :: e, e0, weightScale
    
    gweight = exp(-(e - e0)*weightScale/(4.0d0 + abs(e0)))
    return
  end function gweight
  
end module weight

