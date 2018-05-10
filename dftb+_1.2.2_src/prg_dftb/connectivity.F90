!!* Tools to create and update connectivity matrices and sprint coordinates
module Connectivity
#include "assert.h"
#include "allocate.h"
  use Accuracy
  use FileId
  
  implicit none
  
  private

  public :: OConnect
  public :: destroy, initialise, update, setRefBonds, isDisociation, printCoords, getCOM, reinitialise, getDistance

  character(*), parameter :: transOut = "trans.out"
  character(*), parameter :: formatGeoOut = "(I5,F16.8,F16.8,F16.8)"
  character(*), parameter :: title = "(A8,I8,I8,I8,A8,I8)" 


  !!* data type to store connectivity vectors and sprint coordinates to track
  !!* reactive events 
  type OConnect
  private
    integer :: nAtom
    integer :: COMindex
    real(dp), pointer :: C(:,:)
    integer :: oldC     !* compared with C for dissociation
    real(dp), pointer :: d(:,:)
    real(dp), pointer :: dRatio(:,:)
    real(dp), pointer :: dRef(:,:)
    real(dp), pointer :: sprintMat(:,:)
    real(dp), pointer :: sprintCoords(:)
    real(dp), pointer :: transPoint(:,:) 
    real(dp), pointer :: transPoint2(:,:)
    real(dp), pointer :: transPoint3(:,:)
    integer, pointer :: species(:)

    logical :: trans = .false.
    logical :: seccondTrans = .false. 
    logical :: thirdTrans = .false.
    logical :: disociation = .false. 

    integer :: firstTransIdx
    integer :: seccondTransIdx
    integer :: thirdTransIdx
 
   !!* This could be the start of more complicated constraints in the future 
   !!* integer :: pointer :: breakBond(2)
   
  end type OConnect

  interface destroy
    module procedure Connectivity_destroy
  end interface

  interface initialise
    module procedure Connectivity_initialise
  end interface

  interface getCOM
    module procedure Connectivity_getCOM
  end interface

  interface getDistance
    module procedure Connectivity_getDistance
  end interface

  interface update
    module procedure Connectivity_update
  end interface

  !!*  interface sprint
  !!*  module procedure sprint_coords
  !!*end interface

  interface setRefBonds
    module procedure Connectivity_setRefBonds
  end interface

  interface isDisociation
    module procedure Connectivity_isDisociation
  end interface

  interface printCoords
    module procedure Connectivity_printCoords
  end interface
 
interface reinitialise
module procedure Connectivity_reinitialise
end interface


contains

  subroutine Connectivity_destroy(self)
    type(OConnect), pointer    :: self

    DEALLOCATE(self%C)
    DEALLOCATE(self%d)
    DEALLOCATE(self%dRatio)
    DEALLOCATE(self%sprintMat)
    DEALLOCATE(self%sprintCoords)
    DEALLOCATE(self%transPoint)
    DEALLOCATE(self%transPoint2)
    DEALLOCATE(self%transPoint3)
    DEALLOCATE(self%species)
    DEALLOCATE(self%dRef)     
    DEALLOCATE(self)
    
    end subroutine Connectivity_destroy

  !!* Populates C, d and dRef vectors from DOI: 10.1002/jcc.23790
  !!* Also sets up connectivity matrix specifically for SPRINT COORDS

  subroutine Connectivity_getCOM(self, Coord, masses, COM1, COM2, distance, constraint)
    type(OConnect), pointer    :: self
    real(dp), intent(in)      :: Coord(:,:)
    real(dp), intent(in)      :: masses(:,:)
    real(dp), intent(inout)      :: COM1(:)
    real(dp), intent(inout)      :: COM2(:)
    real(dp), intent(out)      :: distance
    real(dp), intent(out)      :: constraint(:,:)
    
    integer :: i
    integer :: ii
    integer :: jj
    real(dp) :: totalMass1
    real(dp) :: totalMass2
    real(dp) :: COMdistance

    do i = 1, self%COMindex
    totalMass1 = totalMass1 + masses(1,i)
    end do

    do i = (self%COMindex + 1), self%nAtom
          totalMass2 = totalMass2 + masses(1,i)
    end do 

    do ii = 1, self%COMindex
       do jj = 1, (self%nAtom - self%COMindex)
          COM1(jj) = COM1(jj)+(masses(1,ii)*Coord(jj,ii))
       end do
    end do
    COM1(:) = COM1(:)*(1/totalMass1)

    do ii = (self%COMindex + 1), self%nAtom
       do jj = 1, (self%nAtom - self%COMindex)
          COM2(jj) = COM2(jj)+(masses(1,ii)*Coord(jj,ii))
       end do
    end do
    COM2(:) = COM2(:)*(1/totalMass2)

    COMdistance =(COM1(1) - COM2(1))**2
    COMdistance = COMdistance + (COM1(2) - COM2(2))**2
    COMdistance = COMdistance + (COM1(3) - COM2(3))**2
    distance = sqrt(COMdistance)


    do ii = 1, self%nAtom
       do jj = 1, 3
          constraint(jj,ii) = 1 / (2 * distance)
          constraint(jj,ii) = constraint(jj,ii) * 2 * (COM1(jj) - COM2(jj))
          if ( ii .gt. self%COMindex ) then
          constraint(jj,ii) = constraint(jj,ii) * (-1.0_dp *masses(1,ii))/totalMass2 
          else
          constraint(jj,ii) = constraint(jj,ii) * masses(1,ii)/totalMass1
          end if           
       end do
    end do
  end subroutine Connectivity_getCOM

  subroutine Connectivity_getDistance(self, Coord, masses, COM1, COM2, distance, constraint)
    type(OConnect), pointer    :: self
    real(dp), intent(in)      :: Coord(:,:)
    real(dp), intent(in)      :: masses(:,:)
    real(dp), intent(inout)      :: COM1(:)
    real(dp), intent(inout)      :: COM2(:)
    real(dp), intent(out)      :: distance
    real(dp), intent(out)      :: constraint(:,:)

    integer :: i
    integer :: ii
    integer :: jj
    real(dp) :: totalMass1
    real(dp) :: totalMass2
    real(dp) :: COMdistance


    COMdistance =(Coord(1,self%COMindex) - Coord(1,self%COMindex + 1))**2
    COMdistance = COMdistance + (Coord(2,self%COMindex) - Coord(2,self%COMindex + 1))**2
    COMdistance = COMdistance + (Coord(3,self%COMindex) - Coord(3,self%COMindex + 1))**2
    distance = sqrt(COMdistance)

    do ii = 1, self%nAtom
        do jj = 1, 3
            write(*,*) "ii , jj ", ii, jj
            if ( ii == self%COMindex ) then
                constraint(jj,ii) = ( Coord(jj,self%COMindex) - Coord(jj,self%COMindex + 1 ) ) / distance
            else if ( ii == self%COMindex + 1) then
                constraint(jj,ii) = ( Coord(jj,self%COMindex + 1) - Coord(jj,self%COMindex) ) / distance
            else
                constraint(jj,ii) = 0.0_dp
            end if
        end do
    end do





  end subroutine Connectivity_getDistance

  !!* Populates C, d and dRef vectors from DOI: 10.1002/jcc.23790
  !!* Also sets up connectivity matrix specifically for SPRINT COORDS

  subroutine Connectivity_initialise(self, Coord)
    type(OConnect), pointer    :: self
    real(dp), intent(in)      :: Coord(:,:)

    integer :: ii
    integer :: jj

    INITALLOCATE_PARR(self%C,(self%nAtom,self%nAtom))
    INITALLOCATE_PARR(self%d,(self%nAtom,self%nAtom))
    INITALLOCATE_PARR(self%dRatio,(self%nAtom,self%nAtom))
    INITALLOCATE_PARR(self%transPoint,(3,self%nAtom))
    INITALLOCATE_PARR(self%transPoint2,(3,self%nAtom))
    INITALLOCATE_PARR(self%transPoint3,(3,self%nAtom))
   
    
    do ii = 1, self%nAtom
       do jj = 1, self%nAtom
          self%d(ii,jj) = (Coord(1,ii) - Coord(1,jj))**2 + &
  (Coord(2,ii) - Coord(2,jj))**2 + (Coord(3,ii) - Coord(3,jj))**2
          self%d(ii,jj) = sqrt(self%d(ii,jj))
           self%d(ii,jj) = self%d(ii,jj) * 0.529177_dp
          !! Maybe dRef can be assigned from starting Geom 
          !! No atempted Currently
       end do
    end do

   self%dRatio = self%d / self%dRef 
    
    where(self%d(:,:) < self%dRef(:,:))
         self%C(:,:) = 1.0_dp
    elsewhere
         self%C(:,:) = 0.0_dp
    end where
    

!    self%sprintMat(:,:) = (1.0_dp - (self%d(:,:)/self%dRef(:,:))**6.0_dp)  &
!          / (1.0_dp - (self%d(:,:)/self%dRef(:,:))**12.0_dp)

  end subroutine Connectivity_initialise

!!* Populates C, d and dRef vectors from DOI: 10.1002/jcc.23790
!!* Also sets up connectivity matrix specifically for SPRINT COORDS

subroutine Connectivity_reinitialise(self, Coord)
type(OConnect), pointer    :: self
real(dp), intent(in)      :: Coord(:,:)

integer :: ii
integer :: jj

do ii = 1, self%nAtom
do jj = 1, self%nAtom
self%d(ii,jj) = (Coord(1,ii) - Coord(1,jj))**2 + &
(Coord(2,ii) - Coord(2,jj))**2 + (Coord(3,ii) - Coord(3,jj))**2
self%d(ii,jj) = sqrt(self%d(ii,jj))
self%d(ii,jj) = self%d(ii,jj) * 0.529177_dp
!! Maybe dRef can be assigned from starting Geom
!! No atempted Currently
end do
end do

self%dRatio = self%d / self%dRef

where(self%d(:,:) < self%dRef(:,:))
self%C(:,:) = 1.0_dp
elsewhere
self%C(:,:) = 0.0_dp
end where


!    self%sprintMat(:,:) = (1.0_dp - (self%d(:,:)/self%dRef(:,:))**6.0_dp)  &
!          / (1.0_dp - (self%d(:,:)/self%dRef(:,:))**12.0_dp)

end subroutine Connectivity_reinitialise



  !!* Populates C, d and dRef vectors from DOI: 10.1002/jcc.23790
  !!* Also sets up connectivity matrix specifically for SPRINT COORDS

  subroutine Connectivity_update(self, Coord, reaxCount,  transition,  RevIdx )
    type(OConnect), pointer    :: self
    real(dp), intent(in)      :: Coord(:,:)
    integer, intent(in)      :: reaxCount
    logical, intent(out)     :: transition
    integer, intent(out)      :: RevIdx(:)

    real(dp) :: bond    = 0.0_dp
    real(dp) :: nonBond = 0.0_dp

    integer :: ii
    integer :: jj
    integer :: dd
    logical :: mask(self%nAtom, self%nAtom)
    integer :: atom1
    integer :: atom2
    integer :: atom3
    
    do ii = 1, self%nAtom
       do jj = 1, self%nAtom
          self%d(ii,jj) = (Coord(1,ii) - Coord(1,jj))**2 + &
  (Coord(2,ii) - Coord(2,jj))**2 + (Coord(3,ii) - Coord(3,jj))**2
          self%d(ii,jj) = sqrt(self%d(ii,jj))
           self%d(ii,jj) = self%d(ii,jj) * 0.529177_dp
       end do
    end do
    
    self%dRatio = self%d / self%dRef
    
    !!* If the end of transition window is reached, reset everything
    if(reaxCount == 1) then
       self%trans = .false.
       self%seccondTrans = .false.
       self%thirdTrans = .false.
       self%oldC = 0
       revIdx(:) = 0
    end if
    if (reaxCount == 1) then
         where(self%d(:,:) < self%dRef(:,:))
             self%C(:,:) = 1.0_dp
         elsewhere
             self%C(:,:) = 0.0_dp
         end where
    end if
   
    do ii = 1, self%nAtom
    bond    = 0.0_dp
    nonBond    = 100.0_dp
    do jj = 1, self%nAtom
       if(self%C(ii,jj) == 1.0_dp .AND. self%dRatio(ii,jj) > bond) then
            bond = self%dRatio(ii,jj)
            atom1 = ii 
            atom2 = jj
       else if(self%C(ii,jj) == 0.0_dp .AND. self%dRatio(ii,jj) < nonBond) then
            nonBond = self%dRatio(ii,jj)
            atom3 = jj
    end if
    end do
    if (bond > nonBond) then
          self%trans = .true.
          transition = .true.
          self%transPoint(:,:) = Coord(:,:)
          !count bonds prior to transition
          mask = self%C == 1
          self%oldC = count(mask)
    end if
    end do

   
    revIdx(atom1) = 1
    revIdx(atom2) = 1
    revIdx(atom3) = 1

    

!    self%sprintMat(:,:) = (1.0_dp - (self%d(:,:)/self%dRef(:,:))**6.0_dp)  &
!          / (1.0_dp - (self%d(:,:)/self%dRef(:,:))**12.0_dp)

  end subroutine Connectivity_update

  !!* Populates vector of reference bonds from DOI: 10.1002/jcc.23790
  !!* These values are currently hardwired and this needs to be altered in the
  !!* future
  subroutine Connectivity_setRefBonds(self, specie0, nAtom, COMindex)
  type(OConnect), pointer    :: self
  integer, intent(in)      :: specie0(:)
  integer, intent(in)      :: nAtom
  integer, intent(in)      :: COMindex
  integer :: ii
  integer :: jj
  
  INITALLOCATE_P(self)
  INITALLOCATE_PARR(self%dRef,(nAtom,nAtom))
  INITALLOCATE_PARR(self%species,(nAtom))

  self%nAtom = size(specie0)
  self%COMindex = COMindex

  self%species = specie0
  
 

  !!* Coded quickly and poorly clean up in future
  !!* For now this iterates through all pairs of atoms and assigns 
  !!* an ideal bond length based upon type (MASS)
  !!* Only H C and O so far, all values hardcoded

  do ii=1,self%nAtom
     do jj=1,self%nAtom
     if (self%species(ii) == 1 ) then
        if (self%species(jj) == 1) then
           self%dRef(ii,jj) = 0.8_dp
        else if (self%species(jj) == 2) then
           self%dRef(ii,jj) = 1.2_dp
        else if (self%species(jj) == 3) then
           self%dRef(ii,jj) = 1.2_dp
        end if
     else if (self%species(ii) == 2 ) then
        if (self%species(jj) == 1) then
           self%dRef(ii,jj) = 1.2_dp
        else if (self%species(jj) == 2) then
           self%dRef(ii,jj) = 1.6_dp
        else if (self%species(jj) == 3) then
           self%dRef(ii,jj) = 1.6_dp
        end if
     else if (self%species(ii) == 3_dp ) then
        if (self%species(jj) == 1) then
           self%dRef(ii,jj) = 1.2_dp
        else if (self%species(jj) == 2) then
           self%dRef(ii,jj) = 1.6_dp
        else if (self%species(jj) == 3) then
           self%dRef(ii,jj) = 1.6_dp
        end if
     end if
     end do
  end do

 end subroutine Connectivity_setRefBonds

  !!* Compares number of bonds pre and post reaction to check
  !!* whether dissociation event occured


  subroutine Connectivity_isDisociation(self, dissociation)
  type(OConnect), pointer    :: self
  logical, intent(out)      :: dissociation

  logical :: mask(self%nAtom, self%nAtom)

  mask = self%C == 1
  if (count(mask) < self%oldC) then
     dissociation = .true.
  end if

end subroutine Connectivity_isDisociation

 !!* Populates C, d and dRef vectors from DOI: 10.1002/jcc.23790
 !!* Also sets up connectivity matrix specifically for SPRINT COORDS

  subroutine Connectivity_printCoords(self, geomStep, fdTrans)
  type(OConnect), pointer    :: self
  integer, intent(in)      :: geomStep
  integer, intent(in)      :: fdTrans
  integer :: ii 

 close(fdTrans, status='keep')
 open(fdTrans, file="trans.out", position="append", status="old")

 write (fdTrans,title) "Atoms", self%firstTransIdx, self%seccondTransIdx, &
        &self%thirdTransIdx, "GeomStep", geomStep

  write (fdTrans,'(A20)') "First transition"    

  do ii = 1, self%nAtom
      write(fdTrans,formatGeoOut) self%species(ii), &
      &self%transPoint(1,ii), self%transPoint(2,ii), self%transPoint(3,ii)
  end do

  if (self%seccondTrans) then
  write (fdTrans,'(A20)') "Seccond transition"
  do ii = 1, self%nAtom
      write(fdTrans,formatGeoOut) self%species(ii), &
      &self%transPoint2(1,ii), self%transPoint2(2,ii), self%transPoint2(3,ii)
  end do
  end if

  if (self%thirdTrans) then
  write (fdTrans,'(A20)') "Third transition"
  do ii = 1, self%nAtom
      write(fdTrans,formatGeoOut) self%species(ii), &
      &self%transPoint3(1,ii), self%transPoint3(2,ii), self%transPoint3(3,ii)
  end do
  end if

  close(fdTrans, status='keep')

end subroutine Connectivity_printCoords



end module Connectivity

