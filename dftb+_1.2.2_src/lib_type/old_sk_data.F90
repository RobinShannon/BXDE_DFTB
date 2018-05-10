!!* Contains type for representing the data stored in the old SK-file format and
!!* subroutines to read that data from file.
module OldSKData
#include "allocate.h"
#include "assert.h"  
  use Accuracy
  use Constants
  use RepSpline, only : TRepSplineIn
  use RepPoly, only : TRepPolyIn
  use FileId
  use Message
  implicit none
  private

  public :: TOldSKData, readFromFile
  
  !!* Represents the Slater-Koster data in an SK file.
  type TOldSKData
    real(dp) :: dist                          !* Grid separation
    integer :: nGrid                          !* Nr. of grid points
    real(dp) :: skSelf(4)                     !* Atomic eigenvalues.
    real(dp) :: skHubbU(4)                    !* Hubbard Us
    real(dp) :: skOcc(4)                      !* Occupations
    real(dp) :: mass                          !* Mass of the atom
    real(dp), pointer :: skHam(:,:)           !* Table for H
    real(dp), pointer :: skOver(:,:)          !* Table for S
  end type TOldSKData


  !!* Reads the data from an SK-file.
  interface readFromFile
    module procedure OldSKData_readFromFile
  end interface


  integer, parameter :: nSKInter = 20     !* nr. of interactions in the SK-file
  integer, parameter :: nSKInterOld = 10  !* nr. of ints. in old SK-file

  !! Mapping between old (spd) and new (spdf) interactions in the SK-table
  integer, parameter :: iSKInterOld(nSKInterOld) &
      & = (/8, 9, 10, 13, 14, 15, 16, 18, 19, 20/)


contains

  !!* Reads the data from an SK-file.
  !!* @param skData Contains the content of the SK-file on exit
  !!* @param fileName Name of the file to read the data from
  !!* @param homo Is it a homonuclear SK-file?
  !!* @param iSp1 Index of 1st interacting specie (for error messages only)
  !!* @param iSp1 Index of 2nd interacting specie (for error messages only)  
  !!* @param repSplineIn Repulsive spline part of the SK-file.
  !!* @param repPolyIn Repulsive polynomial part of the SK-file.
  subroutine OldSKData_readFromFile(skData, fileName, homo, iSp1, iSp2, &
      &repSplineIn, repPolyIn)
    type(TOldSKData), intent(out) :: skData
    character(len=*), intent(in) :: fileName
    logical, intent(in) :: homo
    integer, intent(in), optional :: iSp1, iSp2
    type(TRepSplineIn), intent(out), optional :: repSplineIn
    type(TRepPolyIn), intent(out), optional :: repPolyIn

    integer, save :: file = -1
    character(lc) :: chDummy
    logical :: tExtended, noRep
    integer :: nShell, nInt
    integer :: ii, jj, iGrid
    real(dp) :: rTmp, rTmp2, rDummy
    real(dp) :: coeffs(2:9), polyCutoff
    integer :: iostat

    ASSERT(present(repSplineIn) .eqv. present(iSp1))
    ASSERT(present(iSp1) .eqv. present(iSp2))

    if (file == -1) then
      file = getFileId()
    end if
    open(file, file=fileName, status="old", action="read", iostat=iostat)
    call checkIoError(iostat, fileName, "Unable to open file")
    rewind(file)
    
    read (file, '(A1)', iostat=iostat) chDummy
    call checkIoError(iostat, fileName, "Unable to read 1st line")
    if (chDummy == '@') then
      tExtended = .true.
      nShell = 4
    else
      tExtended = .false.
      nShell = 3
      rewind(file)
    end if

    read (file,*, iostat=iostat) skData%dist, skData%nGrid
    call checkIoError(iostat, fileName, "Unable to read 1st data line")
    skData%nGrid = skData%nGrid - 1
    if (homo) then
      skData%skSelf(nShell+1:) = 0.0_dp;
      skData%skHubbU(nShell+1:) = 0.0_dp;
      skData%skOcc(nShell+1:) = 0.0_dp;
      read (file,*, iostat=iostat) (skData%skSelf(ii), ii = nShell, 1, -1), &
          &rDummy, &
          &(skData%skHubbU(ii), ii = nShell, 1, -1),&
          &(skData%skOcc(ii), ii = nShell, 1, -1)
      call checkIoError(iostat, fileName, "Unable to read 2nd data line")
      read (file,*, iostat=iostat) skData%mass, (coeffs(ii), ii = 2, 9), &
          &polyCutoff, rDummy, (rDummy, ii = 12, 20)
      call checkIoError(iostat, fileName, "Unable to read 3rd data line")
      skData%mass = skData%mass * amu__au  !convert to atomic units
    else
      read (file,*, iostat=iostat) rDummy, (coeffs(ii), ii = 2, 9), polyCutoff, &
          &(rDummy, ii = 11, 20)
      call checkIoError(iostat, fileName, "Unable to read 1st data line")
    end if

    if (present(repPolyIn)) then
      repPolyIn%polyCoeffs(:) = coeffs(:)
      repPolyIn%cutoff = polyCutoff
    end if

    INITALLOCATE_PARR(skData%skHam, (skData%nGrid, nSKInter))
    INITALLOCATE_PARR(skData%skOver, (skData%nGrid, nSKInter))
    skData%skHam(:,:) = 0.0_dp
    skData%skOver(:,:) = 0.0_dp
    do iGrid = 1, skData%nGrid
      if (tExtended) then
        read (file, *, iostat=iostat) &
            &(skData%skHam(iGrid, ii), ii = 1, nSKInter), &
            &(skData%skOver(iGrid, ii), ii = 1, nSKInter)
        call checkIoError(iostat, fileName, "Reading error for integrals")
      else
        read(file,*, iostat=iostat) &
            &(skData%skHam(iGrid, iSKInterOld(ii)), ii = 1, nSKInterOld), &
            &(skData%skOver(iGrid, iSKInterOld(ii)), ii = 1, nSKInterOld)
        call checkIoError(iostat, fileName, "Reading error for integrals")
      end if
    end do

    if (.not. present(repSplineIn)) then
      close(file)
      return
    end if
    
    do
      read (file, '(A)', iostat=ii) chDummy
      if (ii /= 0) then
        noRep = .true.
        exit
      elseif (chDummy == "Spline") then
        noRep = .false.
        exit
      end if
    end do
    if (noRep) then
      write (chDummy, "(A,I2,A,I2,A)") "Spline repulsive for species pair ", &
          &iSp1, "-", iSp2, "not found."
      call error(chDummy)
    end if

    read (file, *, iostat=iostat) nInt, repSplineIn%cutoff
    call checkIoError(iostat, fileName, "Error in reading repulsive")
    read (file, *, iostat=iostat) (repSplineIn%expCoeffs(ii), ii = 1, 3)
    call checkIoError(iostat, fileName, "Error in reading repulsive")
    
    INITALLOCATE_PARR(repSplineIn%xStart, (nInt))
    INITALLOCATE_PARR(repSplineIn%spCoeffs, (4, nInt - 1))

    rTmp = 0.0_dp
    do jj = 1, nInt
      if (jj == nInt) then
        read (file, *, iostat=iostat) repSplineIn%xStart(jj), rTmp2,  &
            &(repSplineIn%spLastCoeffs(ii), ii = 1, 6)
        call checkIoError(iostat, fileName, "Error in reading repulsive")
      else
        read (file, *, iostat=iostat) repSplineIn%xStart(jj), rTmp2, &
            &(repSplineIn%spCoeffs(ii,jj), ii = 1, 4)
        call checkIoError(iostat, fileName, "Error in reading repulsive")
      end if
      if (jj /= 1 .and. abs(repSplineIn%xStart(jj) - rTmp) > 1e-8_dp) then
        write (chDummy, "(A,I2,A,I2,A)") "Repulsive not continuous for specie &
            &pair ", iSp1, "-", iSp2, "."
        call error(chDummy)
      end if
      rTmp = rTmp2
    end do
    repSplineIn%cutoff = rTmp2
    close(file)

  contains

    subroutine checkIoError(iostat, fileName, msg)
      integer, intent(in) :: iostat
      character(*), intent(in) :: fileName, msg

      if (iostat /= 0) then
        call error("IO error in SK-file '" // trim(fileName) &
            &// "' (" // trim(msg) // ")")
      end if
    
    end subroutine checkIoError

  end subroutine OldSKData_readFromFile
  

end module OldSKData
