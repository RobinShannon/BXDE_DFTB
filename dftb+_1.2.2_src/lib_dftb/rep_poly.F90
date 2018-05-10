!!* Implements a repulsive potential between two atoms represented by
!!* a polynomial of 9th degree
module RepPoly
#include "assert.h"
#include "allocate.h"  
  use Accuracy
  use Bisect
  implicit none
  private

  public :: powMin, powMax
  public :: TRepPolyIn, ORepPoly, init, destruct
  public :: getCutoff, getEnergy, getEnergyDeriv
  
  !! Minimal and maximal power appearing in the polynomial
  integer, parameter :: powMin = 2
  integer, parameter :: powMax = 9
  integer, parameter :: powMin1 = max(powMin, 1)


  !!* Initialisation type for ORepPoly
  type TRepPolyIn
    real(dp) :: polyCoeffs(powMin:powMax)         !* Polynomial coefficients
    real(dp) :: cutoff                            !* Cutoff distance
  end type TRepPolyIn


  !!* Contains the polynomial representation of a repulsive.
  type ORepPoly
    private
    real(dp) :: polyCoeffs(powMin:powMax)
    real(dp) :: cutoff
    logical :: tInit = .false.
  end type ORepPoly


  !!* Initialises polynomial repulsive.
  interface init
    module procedure RepPoly_init
  end interface

  !!* Frees polynomial repulsive.
  interface destruct
    module procedure RepPoly_destruct
  end interface

  !!* Returns cutoff of the repulsive.
  interface getCutoff
    module procedure RepPoly_getCutoff
  end interface

  !!* Returns energy of the repulsive for a given distance.
  interface getEnergy
    module procedure RepPoly_getEnergy
  end interface

  !!* Returns gradient of the repulsive for a given distance.
  interface getEnergyDeriv
    module procedure RepPoly_getEnergyDeriv
  end interface

  
contains

  !!* Initialises polynomial repulsive.
  !!* @param self Polynomial repulsive.
  !!* @param inp Input parameters for the polynomial repulsive.
  subroutine RepPoly_init(self, inp)
    type(ORepPoly), intent(out) :: self
    type(TRepPolyIn), intent(in) :: inp

    
    ASSERT(.not. self%tInit)
    ASSERT(inp%cutoff >= 0.0_dp)

    self%polyCoeffs(:) = inp%polyCoeffs(:)
    self%cutoff = inp%cutoff
    self%tInit = .true.
    
  end subroutine RepPoly_init


  
  !!* Frees polynomial repulsive.
  !!* @param self polynomial repulsive.
  subroutine RepPoly_destruct(self)
    type(ORepPoly), intent(inout) :: self

    ASSERT(self%tInit)
    self%tInit = .false.
    
  end subroutine RepPoly_destruct

  

  !!* Returns cutoff of the repulsive.
  !!* @param self Polynomial repulsive.
  !!* @return Cutoff.
  function RepPoly_getCutoff(self) result(cutoff)
    type(ORepPoly), intent(in) :: self
    real(dp) :: cutoff

    ASSERT(self%tInit)
    cutoff = self%cutoff

  end function RepPoly_getCutoff



  !!* Returns energy of the repulsive for a given distance.
  !!* @param self Polynomial repulsive.
  !!* @param rr Distance between interacting atoms.
  subroutine RepPoly_getEnergy(self, res, rr)
    type(ORepPoly), intent(in) :: self
    real(dp), intent(out) :: res
    real(dp), intent(in) :: rr

    real(dp) :: rrr
    integer :: ii

    ASSERT(self%tInit)

    if (rr >= self%cutoff) then
      res = 0.0_dp
    else
      rrr = self%cutoff - rr
      res = self%polyCoeffs(powMax)
      do ii = powMax-1, powMin, -1
        res = res * rrr + self%polyCoeffs(ii)
      end do
      do ii = powMin, 1, -1
        res = res * rrr
      end do
    end if
    
  end subroutine RepPoly_getEnergy



  !!* Returns gradient of the repulsive for a given distance.
  !!* @param self Polynomial repulsive.
  !!* @param res Resulting contribution
  !!* @param x Actual vector between atoms
  subroutine RepPoly_getEnergyDeriv(self, res, xx)
    type(ORepPoly), intent(in) :: self
    real(dp), intent(out) :: res(3)
    real(dp), intent(in) :: xx(3)

    integer :: ii
    real(dp) :: rr, rrr, xh

    ASSERT(self%tInit)
    
    rr = sqrt(sum(xx**2))
    if (rr >= self%cutoff) then
      res(:) = 0.0_dp
    else
      rrr = self%cutoff - rr
      xh = real(powMax, dp) * self%polyCoeffs(powMax)
      do ii = powMax-1, powMin1, -1
        xh = xh * rrr + real(ii, dp) * self%polyCoeffs(ii)
      end do
      do ii = powMin1, 2, -1
        xh = xh * rrr
      end do
      res(:) = -xh * xx(:) / rr
    end if

  end subroutine RepPoly_getEnergyDeriv

  
end module RepPoly
