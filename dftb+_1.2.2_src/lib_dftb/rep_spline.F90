!!* Implements a repulsive potential between two atoms represented by cubic 
!!* splines.
module RepSpline
#include "assert.h"
#include "allocate.h"  
  use Accuracy
  use Bisect
  implicit none
  private

  public :: TRepSplineIn, ORepSpline, init, destruct
  public :: getCutoff, getEnergy, getEnergyDeriv

  !!* Initialisation type for ORepSpline
  type TRepSplineIn
    real(dp), pointer :: xStart(:) => null() !* Starting pos. for each spline
    real(dp), pointer :: spCoeffs(:,:) => null() !* Spline coeffs (3, nSpline-1)
    real(dp) :: spLastCoeffs(6)         !* Coeffs. of the last spline
    real(dp) :: expCoeffs(3)            !* Coeffs for the exponential head
    real(dp) :: cutoff                  !* Cutoff for the last spline
  end type TRepSplineIn


  !!* Contains the spline representation of a repulsive.
  type ORepSpline
    private
    integer :: nSpline                  ! Nr. of splines.
    real(dp), pointer :: xStart(:) => null()  ! Starting point for each spline
    real(dp), pointer :: spCoeffs(:,:) => null()  ! Spline coeffs (3, nSpline-1)
    real(dp) :: spLastCoeffs(6)         ! Coeffs of the last spline
    real(dp) :: expCoeffs(3)            ! Exponential head
    real(dp) :: cutoff                  ! Cutoff of the last spline
    logical :: tInit = .false.          ! Initialisation status
  end type ORepSpline


  !!* Initialises spline repulsive.
  interface init
    module procedure RepSpline_init
  end interface

  !!* Frees spline repulsive.
  interface destruct
    module procedure RepSpline_destruct
  end interface

  !!* Returns cutoff of the repulsive.
  interface getCutoff
    module procedure RepSpline_getCutoff
  end interface

  !!* Returns energy of the repulsive for a given distance.
  interface getEnergy
    module procedure RepSpline_getEnergy
  end interface

  !!* Returns gradient of the repulsive for a given distance.
  interface getEnergyDeriv
    module procedure RepSpline_getEnergyDeriv
  end interface


contains

  !!* Initialises spline repulsive.
  !!* @param self Spline repulsive.
  !!* @param inp Input parameters for the spline repulsive.
  subroutine RepSpline_init(self, inp)
    type(ORepSpline), intent(out) :: self
    type(TRepSplineIn), intent(in) :: inp


    ASSERT(.not. self%tInit)
    ASSERT(size(inp%xStart) > 0)
    ASSERT(size(inp%spCoeffs, dim=1) == 4)
    ASSERT(size(inp%spCoeffs, dim=2) == size(inp%xStart) - 1)
    ASSERT(inp%cutoff >= 0.0_dp)

    self%nSpline = size(inp%xStart)
    ALLOCATE_PARR(self%xStart, (self%nSpline))
    ALLOCATE_PARR(self%spCoeffs, (4, self%nSpline - 1))
    self%xStart(:) = inp%xStart(:)
    self%spCoeffs(:,:) = inp%spCoeffs(:,:)
    self%spLastCoeffs(:) = inp%spLastCoeffs(:)
    self%expCoeffs(:) = inp%expCoeffs
    self%cutoff = inp%cutoff
    self%tInit = .true.
    
  end subroutine RepSpline_init


  
  !!* Frees spline repulsive.
  !!* @param self Spline repulsive.
  subroutine RepSpline_destruct(self)
    type(ORepSpline), intent(inout) :: self

    DEALLOCATE_PARR(self%xStart)
    DEALLOCATE_PARR(self%spCoeffs)
    self%tInit = .false.
    
  end subroutine RepSpline_destruct

  

  !!* Returns cutoff of the repulsive.
  !!* @param self Spline repulsive.
  !!* @return Cutoff.
  function RepSpline_getCutoff(self) result(cutoff)
    type(ORepSpline), intent(in) :: self
    real(dp) :: cutoff

    cutoff = self%cutoff

  end function RepSpline_getCutoff



  !!* Returns energy of the repulsive for a given distance.
  !!* @param self Spline repulsive.
  !!* @param rr Distance between interacting atoms.
  subroutine RepSpline_getEnergy(self, res, rr)
    type(ORepSpline), intent(in) :: self
    real(dp), intent(out) :: res
    real(dp), intent(in) :: rr
    
    integer :: imatch, ii
    real(dp) :: xh, xv

    !* below this distance, contributions are meaningless, as SK parameters and
    !* repulsive paramers are not defined
    if (rr < minNeighDist) then
      res = 0.0_dp
    elseif (rr < self%xStart(1)) then
      res = exp(-self%expCoeffs(1)*rr + self%expCoeffs(2)) + self%expCoeffs(3)
    elseif (rr > self%cutoff) then
      res = 0.0_dp
    else
      !* find the point in the table to use
      call bisection(imatch, self%xStart, rr)
      
      xv = rr - self%xStart(imatch)
      xh = xv
      if (imatch < self%nSpline) then
        res = self%spCoeffs(1, imatch)
        do ii = 2, 4
          res = res + self%spCoeffs(ii, imatch) * xh
          xh = xh * xv
        end do
      else
        res = self%spLastCoeffs(1)
        do ii = 2, 6
          res = res + self%spLastCoeffs(ii) * xh
          xh = xh * xv
        end do
      end if
    end if
    
  end subroutine RepSpline_getEnergy



  !!* Returns gradient of the repulsive for a given distance.
  !!* @param self Spline repulsive.
  !!* @param res Resulting contribution
  !!* @param x Actual vector between atoms
  subroutine RepSpline_getEnergyDeriv(self, res, xx)
    type(ORepSpline), intent(in) :: self
    real(dp), intent(out) :: res(3)
    real(dp), intent(in) :: xx(3)
    
    integer :: imatch, ii
    real(dp) :: rr, xh, xv

    rr = sqrt(sum(xx**2))
    !* below this distance, contributions are meaningless, as SK parameters and
    !* repulsive paramers are not defined
    if (rr < minNeighDist) then
      res(:) = 0.0_dp
    elseif (rr < self%xStart(1)) then
      res(:) = -self%expCoeffs(1) &
          &* exp(-self%expCoeffs(1) * rr + self%expCoeffs(2))
      res(:) = res(:) * xx(:) / rr
    elseif (rr > self%cutoff) then
      res(:) = 0.0_dp
    else
      !* find the point in the table to use
      call bisection(imatch, self%xStart, rr)
      
      xv = rr - self%xStart(imatch)
      xh = 1.0_dp
      res(:) = 0.0_dp
      if (imatch < self%nSpline) then
        do ii = 2, 4
          res(:) = res(:) + real(ii-1, dp) * self%spCoeffs(ii, imatch) * xh
          xh = xh * xv
        end do
      else
        do ii = 2, 6
          res(:) = res(:) + real(ii-1, dp) * self%spLastCoeffs(ii) * xh
          xh = xh * xv
        end do
      end if
      res(:) = res(:) * xx(:) / rr
    end if

  end subroutine RepSpline_getEnergyDeriv

  
end module RepSpline
