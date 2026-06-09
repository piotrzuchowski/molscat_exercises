
! This module reads in the input_long file, and performs geometric
!   manipulations of dimer configurations

module dimer
 use utility
 use const
 implicit none
 integer, parameter :: maxbond = 8
 character(1), parameter :: monlabelU(2) = ['A','B']
 character(1), parameter :: monlabelL(2) = ['a','b']
 
 integer :: nsite(2), nsitemax, nstype, natom(2), nsptype
 integer :: nmc(2), nmcmax, nmcPhys(2), nmcOptim(2), nmcgroupPhys(2), nmcgroupOptim(2)
 integer :: molchar(2)
 real(8) :: molsize(2), IP(2), DM(2), molPZ(2), molmass(2)!, molrad(2)
 logical :: homog
 integer :: stype(2,MAX_NSITE) ! site type by site
 integer :: sanum(2,MAX_NSITE) ! anum by site
 integer, allocatable :: typeMap(:,:) ! Map from pairs of site types to site pair types
 integer, allocatable :: typeMapInv(:,:) ! Map from site pair types to pairs of site types
 real(8), allocatable :: sit0(:,:,:) ! base geometry
 real(8), allocatable :: sitOA(:,:,:) ! base geometry with updated OA site positions
 integer, allocatable :: bonds(:,:,:) ! usused
 integer, allocatable :: stcount(:,:) ! Number of each site type in each monomer
 character(5), allocatable :: alabel(:,:)
 integer, allocatable :: fftype(:,:) ! site types for initial guiding potential
 integer, allocatable :: stanum(:) ! anum by site type    0 = OA
 logical, allocatable :: incInter(:,:) ! Which site types include which PES components
 integer, allocatable :: mcvect(:,:,:) ! Pairs of site indices defining axis of IDOFs
 integer, allocatable :: mcset(:,:,:) ! Sets of site indices defining what is included in IDOFs
 integer, allocatable :: nmcset(:,:) ! Sizes of sets in mcset
 integer, allocatable :: mctype(:,:) ! type of IDOFs
 integer, allocatable :: mcgroup(:,:) ! group of IDOF. IDOFs in the same group move together with a fixed ratio
 real(8), allocatable :: mcfact(:,:) ! IDOF unit multiplier, for groups with multiple IDOFs
 real(8), allocatable :: rAB0(:,:,:)
 
 real(8) :: sit(2,MAX_NSITE,3)
 !$OMP THREADPRIVATE(sit)
 
contains
  
  ! Reads in the main input file and initializes dimer data
  ! If this subrouine is called without providing an input file, then
  !    sit, nsite, stype, and sanum must first be initialized
  subroutine read_dimer(input_file)
    use utility
    implicit none
    character(*), intent(in), optional :: input_file
    character(1000) :: cbuf
    real(8) :: rbuf(3), mass, val
    integer :: i, j, k
    
    ! Read input file
    if (present(input_file)) then
       open(unit=101,file=trim(input_file))
       do i = 1,2
          read(101,*) nsite(i), nmcPhys(i), nmcOptim(i), molchar(i), IP(i), DM(i), molPZ(i)
          nmc(i) = nmcPhys(i) + nmcOptim(i)
       enddo
    endif
    
    nsitemax = max(nsite(1), nsite(2))
    nmcmax = max(nmc(1), nmc(2))
    allocate(alabel(2,nsitemax))
    allocate(sit0(2,nsitemax,3))
    allocate(sitOA(2,nsitemax,3))
    allocate(fftype(2,nsitemax))
    allocate(bonds(2,nsitemax,0:maxbond))
    allocate(mcvect(2,nmcmax,2))
    allocate(mcset(2,nmcmax,nsitemax))
    allocate(nmcset(2,nmcmax))
    allocate(mctype(2,nmcmax))
    allocate(mcgroup(2,nmcmax))
    allocate(mcfact(2,nmcmax))
    
    sit0(1:2,1:nsitemax,1:3) = sit(1:2,1:nsitemax,1:3)
    sitOA(1:2,1:nsitemax,1:3) = sit(1:2,1:nsitemax,1:3)
    nstype = 0
    natom(:) = 0
    mcgroup(:,:) = 0
    mcfact(:,:) = 1.0d0
    do i = 1,2
       
       bonds(i,:,:) = 0
       do j = 1, nsite(i)
          if (present(input_file)) then
             read(101,*) alabel(i,j),(sit0(i,j,k),k=1,3), &
                  sanum(i,j),stype(i,j),fftype(i,j)!,(bonds(i,j,k),k=1,maxbond)
             
             sit(i,j,1:3) = sit0(i,j,1:3)
             sitOA(i,j,1:3) = sit0(i,j,1:3)
             
             alabel(i,j) = to_upper(alabel(i,j))
             do k = 1,j-1
                if (alabel(i,j) == alabel(i,k)) then
                   write(*,*) 'Error: duplicate site label:',alabel(i,j)
                   stop
                endif
             enddo             
          endif
          
          nstype = max(nstype, stype(i,j))
          
          if (sanum(i,j) > 0) natom(i) = natom(i) + 1
          
          do k = 1, maxbond
             if(bonds(i,j,k).ne.0) bonds(i,j,0) = bonds(i,j,0) + 1
          end do
       enddo
       
       if (present(input_file)) then
          do j = 1,nmc(i)
             read(101,*) mctype(i,j), mcvect(i,j,1), mcvect(i,j,2), nmcset(i,j), mcgroup(i,j), mcfact(i,j)
             if (nmcset(i,j) > 0) read(101,*) mcset(i,j,1:nmcset(i,j))
          enddo
       endif
    enddo
    
    allocate(incInter(nstype,n_FF_comp))
    incInter(:,:) = .true.
    
    if (present(input_file)) then
       do i = 1,nstype
          read(101,*) incInter(i,1:n_FF_comp)
       enddo
       
       close(101)
    endif
    
    do i = 1,2
       nmcgroupPhys(i) = 0
       do j = 1,nmcPhys(i)
          nmcgroupPhys(i) = max(nmcgroupPhys(i), mcgroup(i,j))
       enddo
       
       nmcgroupOptim(i) = 0
       do j = nmcPhys(i)+1,nmc(i)
          nmcgroupOptim(i) = max(nmcgroupOptim(i), mcgroup(i,j))
          do k = 1,nmcset(i,j)
             if (sanum(i,mcset(i,j,k)) /= 0) then
                write(*,*) 'ERROR: A physical atom was selected for off-atomic site optimization'
                stop
             endif
          enddo
       enddo
    enddo
    
    do i = 1,2
       molsize(i) = 0.0d0
       molmass(i) = 0.0d0
       do j = 1,nsite(i)
          if (sanum(i,j) > 0) then
             mass = atomMass(sanum(i,j))
          else
             mass = 0.0d0
          endif
          molmass(i) = molmass(i) + mass
          val = sqrt(sit0(i,j,1)**2 + sit0(i,j,2)**2 + sit0(i,j,3)**2)
          molsize(i) = max(molsize(i), val)
       enddo
    enddo
    
    homog = nsite(1) == nsite(2)
    if (homog) then
       do i = 1,nsite(1)
          do j = 1,3
             if (sit0(1,i,j) /= sit0(2,i,j)) homog = .false.
          enddo
          if (stype(1,i) /= stype(2,i)) homog = .false.
       enddo
    endif
    
    allocate( rAB0(2,nsitemax,nsitemax) )
    do i = 1,2
       do j = 1,nsite(i)
          do k = 1,nsite(i)
             rAB0(i,j,k) = sqrt( &
                  (sit0(i,j,1) - sit0(i,k,1))**2 + &
                  (sit0(i,j,2) - sit0(i,k,2))**2 + &
                  (sit0(i,j,3) - sit0(i,k,3))**2)
          enddo
       enddo
    enddo
    
    allocate(stanum(nstype))
    stanum(:) = -1
    do i = 1,2
       do j = 1,nsite(i)
          if (stanum(stype(i,j)) == -1) then
             stanum(stype(i,j)) = sanum(i,j)
          elseif (stanum(stype(i,j)) /= sanum(i,j)) then
             write(*,*) 'ERROR: Conflicting atom types'
             stop
          endif
       enddo
    enddo
    
    allocate(stcount(2,nstype))
    stcount(:,:) = 0
    do i = 1,2
       do j = 1,nsite(i)
          stcount(i,stype(i,j)) = stcount(i,stype(i,j)) + 1
       enddo
    enddo
    
    nsptype = nstype*(nstype + 1)/2
    allocate( typeMap(nstype,nstype) )
    allocate( typeMapInv(nsptype,2) )
    
    typeMap(:,:) = 0
    typeMapInv(:,:) = 0
    nsptype = 0
    do i = 1,nstype
       do j = i,nstype
          nsptype = nsptype + 1
          typeMap(i,j) = nsptype
          typeMap(j,i) = nsptype
          typeMapInv(nsptype,1) = i
          typeMapInv(nsptype,2) = j
       enddo
    enddo
    
  end subroutine read_dimer
  
  
  ! Wrapper routine for moveDimer that updates sitOA with given coordinates
  subroutine setOACoor(OACoor)
    real(8), intent(in) :: OACoor(:,:)
    real(8) :: mc(2,MAX_IDOF)
    integer :: i, j, k
    
    mc(:,:) = 0.0d0
    do i = 1,2
       do j = 1,nmcgroupOptim(i)
          do k = nmcPhys(i)+1,nmc(i)
             if (mcgroup(i,k) == j) then
                mc(i,k) = mcfact(i,k)*OACoor(i,j)
             endif
          enddo
       enddo
    enddo
    
    sitOA(1:2,1:nsitemax,1:3) = sit0(1:2,1:nsitemax,1:3)
    call moveDimer(0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, mc)
    sitOA(1:2,1:nsitemax,1:3) = sit(1:2,1:nsitemax,1:3)
  end subroutine setOACoor
  
  
  ! Update sit by transforming sitOA according to the given coordinates
  subroutine moveDimer(RCOM, betaA, gammaA, alphaB, betaB, gammaB, mc)
    implicit none
    real(8), intent(in) :: RCOM, betaA, gammaA, alphaB, betaB, gammaB
    real(8), intent(in), optional :: mc(:,:)
    real(8) :: p(3), q(3), offset(3), v(4), u3(3,3), u4(4,4)
    integer :: i, j, k
    
    call rotmat(0.0d0, betaA, gammaA, u3)
    do i = 1,nsite(1)
       sit(1,i,:) = matmul(u3, sitOA(1,i,:))
    enddo
    
    call rotmat(alphaB, betaB, gammaB, u3)
    do i = 1,nsite(2)
       sit(2,i,:) = matmul(u3, sitOA(2,i,:))
       sit(2,i,3) = sit(2,i,3) + RCOM
    enddo
    
    if (present(mc)) then
       do i = 1,2
          do j = 1,nmc(i)
             
             if (abs(mc(i,j)) > 1.0d-20) then
                
                if (mctype(i,j) == MC_CHI) then ! chiral IDOF
                   offset(:) = 0.0d0
                   if (i == 2) offset(3) = RCOM
                   do k = 1,nsite(i)
                      sit(i,k,:) = sit(i,k,:) - 2.0d0*(sit(i,k,:) - offset(:))
                   enddo
                else
                   p(:) = sit(i,mcvect(i,j,1),:)
                   q(:) = sit(i,mcvect(i,j,2),1:3) - sit(i,mcvect(i,j,1),1:3)
                   q(:) = q(:)/sqrt(q(1)**2 + q(2)**2 + q(3)**2)
                   
                   if (mctype(i,j) == MC_ANG) then ! rotational IDOF
                      call rotmat2(p, q, mc(i,j), u4)
                      do k = 1,nmcset(i,j)
                         v(1:3) = sit(i,mcset(i,j,k),:)
                         v(4) = 1.0d0
                         v(:) = matmul(u4, v)
                         sit(i,mcset(i,j,k),:) = v(1:3)
                      enddo
                   elseif (mctype(i,j) == MC_LIN) then ! linear IDOF
                      do k = 1,nmcset(i,j)
                         sit(i,mcset(i,j,k),:) = sit(i,mcset(i,j,k),:) + mc(i,j)*q(:)
                      enddo
                   endif
                endif
                
             endif
             
          enddo
       enddo
    endif
    
  end subroutine moveDimer
  
 ! Returns cartesian distance between sites
 function cDist(mon1, at1, mon2, at2)
   real(8) :: cDist
   integer :: mon1, at1, mon2, at2
   
   cDist = sqrt( (sit(mon1,at1,1) - sit(mon2,at2,1))**2 + &
        (sit(mon1,at1,2) - sit(mon2,at2,2))**2 + &
        (sit(mon1,at1,3) - (sit(mon2,at2,3)))**2 )
 end function cDist
 
 ! Returns sum of VDW radii between given atoms
 function VDWsum(ata, atb)
   real(8) :: VDWSum
   integer, intent(in) :: ata, atb
   
   VDWsum = VDWRad(sanum(1,ata)) + VDWRad(sanum(2,atb))
 end function VDWsum
 
 ! Returns sum of covalent radii between given atoms
 function covSum(ata, atb)
   real(8) :: covSum
   integer, intent(in) :: ata, atb
   
   covSum = covRad(sanum(1,ata)) + covRad(sanum(2,atb))
 end function CovSum
 
 ! Returns distance between closest atoms
 function closestContact()
   real(8) :: closestContact
   integer :: i,j

   closestContact = 1.0d100
   do i = 1,nsite(1)
      do j = 1,nsite(2)
         if (sanum(1,i) > 0 .and. sanum(2,j) > 0) then
            closestContact = min(closestContact, cDist(1,i,2,j))
         endif
      enddo
   enddo
 end function closestContact
 
 ! Returns minimum fraction of covalent distance between any two atoms
 function minCovDist()
   real(8) :: minCovDist
   integer :: i,j
   
   minCovDist = 1.0d100
   do i = 1,nsite(1)
      do j = 1,nsite(2)
         if (sanum(1,i) > 0 .and. sanum(2,j) > 0) then
            minCovDist = min(minCovDist, cDist(1,i,2,j)/covSum(i,j))
         endif
      enddo
   enddo
 end function minCovDist

 ! Returns minimum fraction of VDW distance between any two atoms
 function minVDWDist()
   real(8) :: minVDWDist
   integer :: i,j
   
   minVDWDist = 1.0d100
   do i = 1,nsite(1)
      do j = 1,nsite(2)
         if (sanum(1,i) > 0 .and. sanum(2,j) > 0) then
            minVDWDist = min(minVDWDist, cDist(1,i,2,j)/VDWSum(i,j))
         endif
      enddo
   enddo
 end function minVDWDist
 
 ! Returns COM distance such that closest contact point is given ratio of sum of VDW radii
 function VDWDist(frac)
   real(8) :: VDWDist
   real(8), intent(in) :: frac
   real(8) :: distArr(nstype)
   integer :: i
   
   distArr(:) = 0.0d0
   do i = 1,nstype
      if (stanum(i) > 0) then
         distArr(i) = frac*VDWRad(stanum(i))
      endif
   enddo
   
   VDWDist = minDistSum(distArr)
 end function VDWDist
 
 ! Returns COM distance such that closest contact point is given ratio of sum of covalent radii
 function covDist(frac)
   real(8) :: covDist
   real(8), intent(in) :: frac
   real(8) :: distArr(nstype)
   integer :: i
   
   distArr(:) = 0.0d0
   do i = 1,nstype
      if (stanum(i) > 0) then
         distArr(i) = frac*covRad(stanum(i))
      endif
   enddo
   
   covDist = minDistSum(distArr)
 end function CovDist
 
 ! Returns COM distance such that closest contact point is given distance
 function minDist(inDist)
   real(8) :: minDist
   real(8), intent(in) :: inDist
   real(8) :: distArr(nstype)
   
   distArr(:) = 0.5d0*inDist
   minDist = minDistSum(distArr)
 end function MinDist
 
 ! Returns the difference in COM distance for given orientation such that
 ! closest atoms are separated by sum of corresponding distances in array (all units in bohrs)
 function minDistSum(indist)
   use double
   implicit none
   real(kind=dbl) :: minDistSum
   real(kind=dbl), intent(in) :: indist(nstype)
   real(kind=dbl) :: startDist
   real(kind=dbl) :: R, rab, step, minpairdist, rad1, rad2
   integer :: i, j
   logical :: stage2
   
   startDist = 0.0d0
   step = 1.0d100
   do i = 1,nstype
      startDist = max(startDist, indist(i))
      if (indist(i) > 0.0d0) step = min(step,indist(i))
   enddo
   startDist = 2.0d0*startDist + molsize(1) + molsize(2) + 10.0d0
   step = 0.2d0*step
   
   R = startDist
   stage2 = .false.
   do while (R > step)
      minpairdist = 1.0d100
      do i = 1,nsite(1)
         do j = 1,nsite(2)
            rad1 = indist(stype(1,i))
            rad2 = indist(stype(2,j))
            if (rad1 > 0.0d0 .and. rad2 > 0.0d0) then
               rab = sqrt( (sit(1,i,1) - sit(2,j,1))**2 + &
                    (sit(1,i,2) - sit(2,j,2))**2 + &
                    (sit(1,i,3) - (sit(2,j,3) + R))**2 ) 
               minpairdist = min(minpairdist, rab/(rad1 + rad2))
            endif
         enddo
      enddo
      
      if (minpairdist < 1.0d0) then
         if (stage2) then
            exit
         else
            stage2 = .true.
            R = R + step
            step = step*0.05d0
         endif
      else
         R = R - step
      endif
   enddo
   
   minDistSum = R + 0.5d0*step
 end function MinDistSum
 
 
 ! Symmetry-safe metric for difference between two dimers aka 'epsilon_mn'
 function epsMetric(sit1, sit2)
   real(8) :: epsMetric
   real(8), target :: sit1(:,:,:), sit2(:,:,:)
   
   integer :: i, j, idim, imon1, imon2, ityp1, ityp2, isit1, isit2
   real(8), pointer :: siti(:,:,:)
   real(8) :: distVect(2,2*MAX_NSITE*MAX_NSITE)
   real(8) :: tmp, sum
   integer :: nr, nr_new, minr
   
   do idim = 1,2
      if (idim == 1) then
         siti => sit1
      else
         siti => sit2
      endif
      nr = 0
      
      ! Loop over monomers and dimer
      do imon1 = 1,2
         do imon2 = imon1,2
            if (imon1 /= imon2 .or. nmcPhys(imon1) > 0) then
               
               ! Loop over site pair types
               do ityp1 = 1,nstype
                  do ityp2 = ityp1,nstype
                     nr_new = nr
                     
                     ! Append all distances of given type pair
                     do isit1 = 1,nsite(imon1)
                        do isit2 = 1,nsite(imon2)
                           if ((stype(imon1,isit1) == ityp1 .and. stype(imon2,isit2) == ityp2) .or. &
                                (stype(imon1,isit1) == ityp2 .and. stype(imon2,isit2) == ityp1)) then
                              nr_new = nr_new + 1
                              sum = 0.0d0
                              do i = 1,3
                                 sum = sum + (siti(imon1,isit1,i) - siti(imon2,isit2,i))**2
                              enddo
                              distVect(idim,nr_new) = sqrt(sum)
                           endif
                        enddo
                     enddo
                     
                     ! Sort added distances
                     do i = nr+1,nr_new-1
                        minr = i
                        do j = i+1,nr_new
                           if (distVect(idim,j) < distVect(idim,minr)) minr = j 
                        enddo
                        
                        tmp = distVect(idim,i)
                        distVect(idim,i) = distVect(idim,minr)
                        distVect(idim,minr) = tmp
                     enddo
                     
                     nr = nr_new
                  enddo
               enddo
               
            endif
         enddo
      enddo
      
   enddo
   
   ! Compute RMSD of the two distance arrays
   epsMetric = 0.0d0
   do i = 1,nr
      epsMetric = epsMetric + &
           (distVect(1,i) - distVect(2,i))**2/((1.0d0 + distVect(1,i))*(1.0d0 + distVect(2,i)))
   enddo
   epsMetric = sqrt(epsMetric/nr)
   
 end function epsMetric
 
 subroutine printXYZ()
   integer :: i,j
   
   do i = 1,2
      do j = 1,nsite(i)
         write(*,'(I3,A,3F15.10)') &
              sanum(i,j),' ',sit(i,j,1:3)*bohr2A
      enddo
   enddo
 end subroutine printXYZ
 
end module dimer
