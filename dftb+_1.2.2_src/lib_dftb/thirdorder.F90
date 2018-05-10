!> Routines implementing the full 3rd order DFTB.
module ThirdOrder
#include "assert.h"
#include "allocate.h"
  use Accuracy
  use Short_gamma, only : expGammaCutoff
  use Periodic, only : TNeighborList, getNrOfNeighbors
  implicit none
  private

  public :: TThirdOrderInp, OThirdOrder
  public :: init, getCutoff, updateCoords, updateCharges
  public :: getShiftAtom, getEnergyPerAtom, addGradientDC

  !> Input for the 3rd order module
  type TThirdOrderInp
    integer :: nAtom                           ! nr. of atoms
    real(dp), allocatable :: hubbU(:)          ! hubb U for each atom
    real(dp), allocatable :: hubbUDeriv(:)     ! dU/dq for each atom
    logical, allocatable :: damped(:)          ! interaction damping? (nAtom)
    real(dp) :: dampExp                        ! Exponential for damping
  end type TThirdOrderInp

  !> Internal status of third order.
  type OThirdOrder
    integer :: nSpecie, nAtom
    real(dp), allocatable :: UU(:)
    real(dp), allocatable :: dUdq(:)
    real(dp), allocatable :: shift1(:), shift2(:)
    real(dp), allocatable :: charges(:)
    real(dp), allocatable :: cutoffs(:,:)
    integer, allocatable :: nNeigh(:,:), nNeighMax(:)
    real(dp), allocatable :: gamma3ab(:,:), gamma3ba(:,:)
    logical, allocatable :: damped(:)
    real(dp) :: xi  ! damping exponent
    real(dp) :: mCutoff
  end type OThirdOrder

  interface init
    module procedure ThirdOrder_init
  end interface

  interface getCutoff
    module procedure ThirdOrder_getCutoff
  end interface

  interface updateCoords
    module procedure ThirdOrder_updateCoords
  end interface

  interface updateCharges
    module procedure ThirdOrder_updateCharges
  end interface

  interface getShiftAtom
    module procedure ThirdOrder_getShiftPerAtom
  end interface

  interface getEnergyPerAtom
    module procedure ThirdOrder_getEnergyPerAtom
  end interface

  interface addGradientDC
    module procedure ThirdOrder_addGradientDC
  end interface


contains

  !> Initializes instance.
  !! \param self Instance.
  !! \param inp Input data.
  subroutine ThirdOrder_init(self, inp)
    type(OThirdOrder), intent(inout) :: self
    type(TThirdOrderInp), intent(in) :: inp

    integer :: iSp1, iSp2

    self%nAtom = inp%nAtom
    self%nSpecie = size(inp%hubbU)
    ALLOCATE_(self%UU, (self%nSpecie))
    ALLOCATE_(self%dUdq, (self%nSpecie))
    self%UU(:) = inp%hubbU
    self%dUdq(:) = inp%hubbUDeriv

    ! Get cutoff (normal SCC cutoff)
    ALLOCATE_(self%cutoffs, (self%nSpecie, self%nSpecie))
    do iSp1 = 1, self%nSpecie
      do iSp2 = iSp1, self%nSpecie
        self%cutoffs(iSp2, iSp1) = &
            & expGammaCutoff(self%UU(iSp2), self%UU(iSp1))
        self%cutoffs(iSp1, iSp2) = self%cutoffs(iSp2, iSp1)
      end do
    end do
    self%mCutoff = maxval(self%cutoffs)
    
    ALLOCATE_(self%nNeigh, (self%nSpecie, self%nAtom))
    ALLOCATE_(self%nNeighMax, (self%nAtom))
    ALLOCATE_(self%shift1, (self%nAtom))
    ALLOCATE_(self%shift2, (self%nAtom))
    ALLOCATE_(self%charges, (self%nAtom))
    ALLOCATE_(self%gamma3ab, (0, self%nAtom))
    ALLOCATE_(self%gamma3ba, (0, self%nAtom))
    ALLOCATE_(self%damped, (self%nSpecie))
    self%damped = inp%damped
    self%xi = inp%dampExp
    
  end subroutine ThirdOrder_init


  !> Returns real space cutoff.
  !! \param self Instance.
  !! \return Cutoff.
  function ThirdOrder_getCutoff(self) result(cutoff)
    type(OThirdOrder), intent(inout) :: self
    real(dp) :: cutoff

    cutoff = self%mCutoff
    
  end function ThirdOrder_getCutoff


  !> Updates data structures if there are changed coordinates for the instance.
  !! \param self Instance.
  !! \param neighList Neighbor list.
  !! \param species Species for all atoms (nAllAtom).
  subroutine ThirdOrder_updateCoords(self, neighList, species)
    type(OThirdOrder), intent(inout) :: self
    type(TNeighborList), intent(in) :: neighList
    integer, intent(in) :: species(:)

    integer :: iNeigh, iAt1, iAt2, iSp1, iSp2
    logical :: damping
    real(dp) :: rr

    self%nNeigh(:,:) = 0.0_dp
    do iAt1 = 1, self%nAtom
      iSp1 = species(iAt1)
      do iSp2 = 1, self%nSpecie
        self%nNeigh(iSp2, iAt1) = &
            &getNrOfNeighbors(neighList, self%cutoffs(iSp2, iSp1), iAt1)
      end do
    end do
    self%nNeighMax = maxval(self%nNeigh, dim=1)
    
    if (size(self%gamma3ab, dim=1) < maxval(self%nNeighMax) + 1) then
      DEALLOCATE_(self%gamma3ab)
      DEALLOCATE_(self%gamma3ba)
      ALLOCATE_(self%gamma3ab, (0:maxval(self%nNeighMax), self%nAtom))
      ALLOCATE_(self%gamma3ba, (0:maxval(self%nNeighMax), self%nAtom))
    end if

    self%gamma3ab = 0.0_dp
    self%gamma3ba = 0.0_dp
    do iAt1 = 1, self%nAtom
      iSp1 = species(iAt1)
      do iNeigh = 0, self%nNeighMax(iAt1)
        iAt2 = neighList%iNeighbor(iNeigh, iAt1)
        iSp2 = species(iAt2)
        if (iNeigh <= self%nNeigh(iSp2, iAt1)) then
          rr = sqrt(neighList%neighDist2(iNeigh, iAt1))
          damping = self%damped(iSp1) .or. self%damped(iSp2)
          self%gamma3ab(iNeigh, iAt1) = gamma3(self%UU(iSp1), &
              &self%UU(iSp2), self%dUdq(iSp1), rr, damping, self%xi)
          self%gamma3ba(iNeigh, iAt1) = gamma3(self%UU(iSp2), &
              &self%UU(iSp1), self%dUdq(iSp2), rr, damping, self%xi)
        end if
      end do
    end do
    
  end subroutine ThirdOrder_updateCoords


  !> Signalizes changed charges for the instance.
  !! \param self Instance.
  !! \param species Species (nAtom).
  !! \param neighList Neighbor list.
  !! \param qq Charges.
  !! \param q0 Reference charges
  !! \param img2CentCell Mapping on atoms in central cell.
  subroutine ThirdOrder_updateCharges(self, species, neighList, qq, q0, &
      &img2CentCell)
    type(OThirdOrder), intent(inout) :: self
    integer, intent(in) :: species(:)
    type(TNeighborList), intent(in) :: neighList
    real(dp), intent(in) :: qq(:,:,:), q0(:,:,:)
    integer, intent(in) :: img2CentCell(:)

    integer :: iAt1, iAt2f, iNeigh

    ASSERT(size(species) == self%nAtom)
    ASSERT(size(qq, dim=2) == self%nAtom)
    ASSERT(size(q0, dim=2) == self%nAtom)

    self%charges = sum(qq(:,:,1) - q0(:,:,1), dim=1)

    self%shift1(:) = 0.0_dp
    self%shift2(:) = 0.0_dp
    do iAt1 = 1, self%nAtom
      do iNeigh = 0, self%nNeighMax(iAt1)
        iAt2f = img2CentCell(neighList%iNeighbor(iNeigh, iAt1))
        self%shift1(iAt1) = self%shift1(iAt1) &
            &+ self%gamma3ab(iNeigh, iAt1) * self%charges(iAt2f)
        self%shift2(iAt1) = self%shift2(iAt1) &
            &+ self%gamma3ba(iNeigh, iAt1) * self%charges(iAt2f)**2
        if (iAt2f /= iAt1) then
          self%shift1(iAt2f) = self%shift1(iAt2f) &
              &+ self%gamma3ba(iNeigh, iAt1) * self%charges(iAt1)
          self%shift2(iAt2f) = self%shift2(iAt2f) &
              &+ self%gamma3ab(iNeigh, iAt1) * self%charges(iAt1)**2
        end if
      end do
      self%shift1(iAt1) = self%shift1(iAt1) * self%charges(iAt1)
    end do
    self%shift1 = 2.0_dp / 3.0_dp * self%shift1
    self%shift2 = 1.0_dp / 3.0_dp * self%shift2

  end subroutine ThirdOrder_updateCharges


  !> Returns shifts per atom.
  !! \param self Instance.
  !! \param shift Shift per atom.
  subroutine ThirdOrder_getShiftPerAtom(self, shift)
    type(OThirdOrder), intent(inout) :: self
    real(dp), intent(out) :: shift(:)

    ASSERT(size(shift) == self%nAtom)

    shift(:) = self%shift1 + self%shift2
    
  end subroutine ThirdOrder_getShiftPerAtom


  !> Returns energy per atom.
  !! \param self Instance.
  !! \param energyPerAtom Energy per atom.
  subroutine ThirdOrder_getEnergyPerAtom(self, energyPerAtom)
    type(OThirdOrder), intent(inout) :: self
    real(dp), intent(out) :: energyPerAtom(:)

    ASSERT(size(energyPerAtom) == self%nAtom)

    energyPerAtom = self%shift2 * self%charges
    
  end subroutine ThirdOrder_getEnergyPerAtom


  !> Add gradient component resulting from the derivative of the potential.
  !! \param self Instance.
  !! \param neighList Neighbor list.
  !! \param species Specie for each atom.
  !! \param coords Coordinate of each atom.
  !! \param img2CentCell Mapping of atoms to cetnral cell.
  !! \param derivs Gradient on exit.
  subroutine ThirdOrder_addGradientDC(self, neighList, species, coords, &
      &img2CentCell, derivs)
    type(OThirdOrder), intent(inout) :: self
    type(TNeighborList), intent(in) :: neighList
    integer, intent(in) :: species(:)
    real(dp), intent(in) :: coords(:,:)
    integer, intent(in) :: img2CentCell(:)
    real(dp), intent(inout) :: derivs(:,:)

    integer :: iAt1, iAt2, iAt2f, iSp1, iSp2, iNeigh
    real(dp) :: rab, tmp, tmp3(3)
    logical :: damping
   
    do iAt1 = 1, self%nAtom
      iSp1 = species(iAt1)
      do iNeigh = 1, self%nNeighMax(iAt1)
        iAt2 = neighList%iNeighbor(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)
        iSp2 = species(iAt2f)
        if (iAt1 == iAt2f) then
          continue
        end if
        if (iNeigh <= self%nNeigh(iSp2, iAt1)) then
          rab = sqrt(neighList%neighDist2(iNeigh, iAt1))
          damping = self%damped(iSp1) .or. self%damped(iSp2)
          tmp = self%charges(iAt1) / 3.0_dp &
              &* (self%charges(iAt1) * self%charges(iAt2f) &
              &* gamma3pR(self%UU(iSp1), self%UU(iSp2), self%dUdQ(iSp1), rab, &
              &damping, self%xi)&
              &+ self%charges(iAt2f)**2 &
              &* gamma3pR(self%UU(iSp2), self%UU(iSp1), self%dUdQ(iSp2), rab, &
              &damping, self%xi))
          tmp3 = tmp * (coords(:, iAt1) - coords(:, iAt2)) / rab
          derivs(:, iAt1) = derivs(:, iAt1) + tmp3
          derivs(:, iAt2f) = derivs(:, iAt2f) - tmp3
        end if
      end do
    end do
   
  end subroutine ThirdOrder_addGradientDC

          

!!! Helper functions

  ! Gamma_AB = dgamma_AB/dUa * (dUa/dQa)
  function gamma3(Ua, Ub, dUa, rab, damping, xi) result(res)
    real(dp), intent(in) :: Ua, Ub, dUa, rab
    logical, intent(in) :: damping
    real(dp), intent(in) :: xi
    real(dp) :: res

    res = gamma2pU(Ua, Ub, rab, damping, xi) * dUa

  end function gamma3


  ! dGamma_AB/dr
  function gamma3pR(Ua, Ub, dUa, rab, damping, xi) result(res)
    real(dp), intent(in) :: Ua, Ub, dUa, rab
    logical, intent(in) :: damping
    real(dp), intent(in) :: xi
    real(dp) :: res

    res = gamma2pUpR(Ua, Ub, rab, damping, xi) * dUa

  end function gamma3pR
    
    
  ! dgamma_AB/dUa
  ! Sign convention: routine delivers dgamma_AB/dUa with the right sign.
  ! Energy contributions must be therefore summed with *positive* sign.
  function gamma2pU(Ua, Ub, rab, damping, xi) result(res)
    real(dp), intent(in) :: Ua, Ub, rab
    logical, intent(in) :: damping
    real(dp), intent(in) :: xi
    real(dp) :: res

    real(dp) :: tauA, tauB, tau, uu

    tauA = 3.2_dp * Ua      ! 16/5 * Ua
    tauB = 3.2_dp * Ub

    if (rab < tolSameDist) then
      ! NOTE: 0.5 to take care about the factor 1/2 in Gamma3_{aa}.
      res = 0.5_dp
    else if (abs(Ua - Ub) < minHubDiff) then
      tau = 0.5_dp * (tauA + tauB)
      res = -3.2_dp * shortpT_2(tau, rab)
      if (damping) then
        uu = 0.5_dp * (Ua + Ub)
        res = res * hh(uu, uu, rab, xi) &
            &- short_2(tau, rab) * hpU(uu, uu, rab, xi)
      end if
    else
      res = -3.2_dp * shortpT_1(tauA, tauB, rab)
      if (damping) then
        res = res * hh(Ua, Ub, rab, xi) &
            &- short_1(tauA, tauB, rab) * hpU(Ua, Ub, rab, xi)
      end if
    end if

  end function gamma2pU


  ! d^2gamma_AB/dUa*dr
  ! Sign convention: routine delivers d^2gamma_AB/dUa*dr with the right sign.
  ! Gradient contributions must be therefore summed with *positive* sign.
  function gamma2pUpR(Ua, Ub, rab, damping, xi) result(res)
    real(dp), intent(in) :: Ua, Ub, rab
    logical, intent(in) :: damping
    real(dp), intent(in) :: xi
    real(dp) :: res

    real(dp) :: tauA, tauB, tau, uu

    tauA = 3.2_dp * Ua
    tauB = 3.2_dp * Ub
    
    if (rab < tolSameDist) then
      res = 0.0_dp
    else if (abs(Ua - Ub) < minHubDiff) then
      tau = 0.5_dp * (tauA + tauB)
      res = -3.2_dp * shortpTpR_2(tau, rab)
      if (damping) then
        uu = 0.5_dp * (Ua + Ub)
        res = res * hh(uu, uu, rab, xi) &
            &- 3.2_dp * shortpT_2(tau, rab) * hpR(uu, uu, rab, xi) &
            &- shortpR_2(tau, rab) * hpU(uu, uu, rab, xi)&
            &- short_2(tau, rab) * hpUpR(uu, uu, rab, xi)
      end if
    else
      res = -3.2_dp *shortpTpR_1(tauA, tauB, rab)
      if (damping) then
        res = res * hh(Ua, Ub, rab, xi) &
            &- 3.2_dp * shortpT_1(tauA, tauB, rab) * hpR(Ua, Ub, rab, xi)&
            &- shortpR_1(tauA, tauB, rab) * hpU(Ua, Ub, rab, xi)&
            &- short_1(tauA, tauB, rab) * hpUpR(Ua, Ub, rab, xi)
      end if
    end if

  end function gamma2pUpR
  


!!! Helper function


  ! S1(tauA,tauB,r): Short range SCC when tauA <> tauB and r <> 0
  function short_1(tauA, tauB, rab) result(res)
    real(dp), intent(in) :: tauA, tauB, rab
    real(dp) :: res

    res = exp(-tauA * rab) * ff(tauA, tauB, rab) &
        &+ exp(-tauB * rab) * ff(tauB, tauA, rab)

  end function short_1

  
  ! S2(tau,r), short range SCC when tauA = tauB = tau and r <> 0.
  function short_2(tau, rab) result(res)
    real(dp), intent(in) :: tau, rab
    real(dp) :: res

    res = exp(-tau * rab) * gg(tau, rab)

  end function short_2


  ! dS1(tauA,tauB,r)/dtauA
  function shortpT_1(tauA, tauB, rab) result(res)
    real(dp), intent(in) :: tauA, tauB, rab
    real(dp) :: res

    res = exp(-tauA * rab) * (fpT1(tauA, tauB, rab) &
          &- rab * ff(tauA, tauB, rab)) &
          &+ exp(-tauB * rab) * fpT2(tauB, tauA, rab)

  end function shortpT_1


  ! dS2(tauA,tauB,r)/dtauA
  function shortpT_2(tau, rab) result(res)
    real(dp), intent(in) :: tau, rab
    real(dp) :: res

    res = exp(-tau * rab) * (gpT(tau, rab) - rab * gg(tau, rab))

  end function shortpT_2


  ! dS1(tauA,tauB,r)/dr
  function shortpR_1(tauA, tauB, rab) result(res)
    real(dp), intent(in) :: tauA, tauB, rab
    real(dp) :: res

    res = exp(-tauA * rab) * (fpR(tauA, tauB, rab) &
        &- tauA *  ff(tauA, tauB, rab)) &
        &+ exp(-tauB * rab) * (fpR(tauB, tauA, rab) &
        &- tauB * ff(tauB, tauA, rab))

  end function shortpR_1
    

  ! dS2(tauA,tauB,r)/dr
  function shortpR_2(tau, rab) result(res)
    real(dp), intent(in) :: tau, rab
    real(dp) :: res

    res = exp(-tau * rab) * (gpR(tau, rab) - tau * gg(tau, rab))

  end function shortpR_2


  ! d^2S1(tauA,tauB,r)/dtauA*dr
  function shortpTpR_1(tauA, tauB, rab) result(res)
    real(dp), intent(in) :: tauA, tauB, rab
    real(dp) :: res

    res = exp(-tauA * rab) * (ff(tauA, tauB, rab) * (tauA * rab - 1.0_dp) &
        &- tauA * fpT1(tauA, tauB, rab) + fpT1pR(tauA, tauB, rab) &
        &- rab * fpR(tauA, tauB, rab)) &
        &+ exp(-tauB * rab) * (fpT2pR(tauB, tauA, rab) &
        &- tauB * fpT2(tauB, tauA, rab))

  end function shortpTpR_1


  ! d^2S2(tau,r)/dtau*dr
  function shortpTpR_2(tau, rab) result(res)
    real(dp), intent(in) :: tau, rab
    real(dp) :: res

    res = exp(-tau * rab) * ((tau * rab - 1.0_dp) * gg(tau, rab) &
        &- tau * gpT(tau, rab) + gpTpR(tau, rab) - rab * gpR(tau, rab))

  end function shortpTpR_2

  
  ! f(tauA,tauB,r)
  function ff(tauA, tauB, rab) result(res)
    real(dp), intent(in) :: tauA, tauB, rab
    real(dp) :: res

    res = 0.5_dp * tauA * tauB**4 / (tauA**2 - tauB**2)**2 &
        &- (tauB**6 - 3.0_dp * tauA**2 * tauB**4) &
        &/ ((tauA**2 - tauB**2)**3 * rab)

  end function ff


  ! df(tauA,tauB,r)/dtauA
  function fpT1(tauA, tauB, rab) result(res)
    real(dp), intent(in) :: tauA, tauB, rab
    real(dp) :: res

    res = -0.5_dp * (tauB**6 + 3.0_dp * tauA**2 * tauB**4) &
        &/ (tauA**2 - tauB**2)**3 &
        &- 12.0_dp * tauA**3 * tauB**4 &
        &/ ((tauA**2 - tauB**2)**4 * rab)

  end function fpT1


  ! df(tauB,tauA,rab)/dtauA
  function fpT2(tauB, tauA , rab) result(res)
    real(dp), intent(in) :: tauA, tauB, rab
    real(dp) :: res

    res = 2.0_dp * tauB**3 * tauA**3 / (tauB**2 - tauA**2)**3 &
        &+ 12.0_dp * tauB**4 * tauA**3 / ((tauB**2 - tauA**2)**4 * rab)

  end function fpT2
  

  ! df(tauA, tauB,r)/dr
  function fpR(tauA, tauB, rab) result(res)
    real(dp), intent(in) :: tauA, tauB, rab
    real(dp) :: res

    res = (tauB**6 - 3.0_dp * tauB**4 * tauA**2) &
        &/ (rab**2 * (tauA**2 - tauB**2)**3)
    
  end function fpR


  ! d^2f(tauA,tauB,r)/dtauA*dr
  function fpT1pR(tauA, tauB, rab) result(res)
    real(dp), intent(in) :: tauA, tauB, rab
    real(dp) :: res

    res = 12.0_dp * tauA**3 * tauB**4 / (rab**2 * (tauA**2 - tauB**2)**4)

  end function fpT1pR


  ! d^2f(tauB,tauA,r)/dtauA*dr
  function fpT2pR(tauB, tauA, rab) result(res)
    real(dp), intent(in) :: tauB, tauA, rab
    real(dp) :: res

    res = -12.0_dp * tauA**3 * tauB**4 / (rab**2 * (tauA**2 - tauB**2)**4)

  end function fpT2pR
  

  ! g(tau,r)
  function gg(tau, rab) result(res)
    real(dp), intent(in) :: tau, rab
    real(dp) :: res

    res = 1.0_dp / (48.0_dp * rab) * (48.0_dp + 33.0_dp * tau * rab &
        &+ 9.0_dp * tau**2 * rab**2 + tau**3 * rab**3)

  end function gg


  ! dg(tau,rab)/dtau
  function gpT(tau, rab) result(res)
    real(dp), intent(in) :: tau, rab
    real(dp) :: res

    res = 1.0_dp / 48.0_dp * (33.0_dp + 18.0_dp * tau * rab &
        &+ 3.0_dp * tau**2 * rab**2)

  end function gpT

  
  ! dg(tau,r)/dr
  function gpR(tau, rab) result(res)
    real(dp), intent(in) :: tau, rab
    real(dp) :: res

    res = (-48.0_dp + 9.0_dp * (tau * rab)**2 + 2.0_dp * (tau * rab)**3) &
        &/ (48.0_dp * rab**2)
    
  end function gpR

  
  ! d^2g(tau,r)/dtau*dr
  function gpTpR(tau, rab) result(res)
    real(dp), intent(in) :: tau, rab
    real(dp) :: res

    res = (3.0_dp * tau + tau**2 * rab) / 8.0_dp
    
  end function gpTpR
  

  ! Damping: h(Ua,Ub)
  function hh(Ua, Ub, rab, xi) result(res)
    real(dp), intent(in) :: Ua, Ub, rab, xi
    real(dp) :: res

    res = exp(-(0.5_dp * (Ua + Ub))**xi * rab**2)

  end function hh


  ! dh(Ua,Ub)/dUa
  function hpU(Ua, Ub, rab, xi) result(res)
    real(dp), intent(in) :: Ua, Ub, rab, xi
    real(dp) :: res

    res = -0.5_dp * xi * rab**2 * (0.5_dp * (Ua + Ub))**(xi - 1.0_dp) &
        &* hh(Ua, Ub, rab, xi)

  end function hpU


  ! dh(Ua,Ub)/dr
  function hpR(Ua, Ub, rab, xi) result(res)
    real(dp), intent(in) :: Ua, Ub, rab, xi
    real(dp) :: res

    res = -2.0_dp * rab * (0.5_dp * (Ua + Ub))**xi * hh(Ua, Ub, rab, xi)

  end function hpR


  ! dh(Ua,Ub)/dUa*dr
  function hpUpR(Ua, Ub, rab, xi) result(res)
    real(dp), intent(in) :: Ua, Ub, rab, xi
    real(dp) :: res

    res = xi * rab * (0.5_dp * (Ua + Ub))**(xi - 1.0_dp) &
        &* (rab**2 * (0.5_dp * (Ua + Ub))**xi - 1.0_dp) * hh(Ua, Ub, rab, xi)
    
  end function hpUpR
  
  
end module ThirdOrder
