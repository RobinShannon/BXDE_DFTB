!!* Dispersion a la Slater-Kirkwood as implemented by M. Elstner in old DFTB.
!!* @desc The expression as found in the old DFTB had been reimplemented.
!!*   The periodic case had been completely rewritten using a correct
!!*   Ewald summation (instead of the pure real space summation, which converges
!!*   very poorly).
!!* @note The expression for C6(iAt1,iAt2) is not the same as in the reference 
!!*   paper by M. Elstner, but the one found in the old code (implemented by
!!*   him as well). Furthermore, the expression for the dispersion energy in the
!!*   paper (eq. 9) misses a factor of 1/2.
!!* @ref Elstner et al., J. Chem. Phys., 114, 5149 (2001)
!!* @todo The generation of the reciprocal lattice vectors should not be done
!!*   localy, but somewhere outside, since the Coulomb module does the same.
module DispSlaterKirkwood
#include "allocate.h"
#include "assert.h"
  use Accuracy
  use Periodic, only: TNeighborList, getNrOfNeighborsForAll, getLatticePoints
  use Constants, only : pi
  use DispCommon
  implicit none
  private

  public :: TDispSlaKirkInp, ODispSlaKirk
  public :: init, destruct, updateCoords, updateLatVecs
  public :: getEnergies, addGradients
  public :: getRCutoff


  !!* Contains the initialisation data for the Slater-Kirkwood module
  type TDispSlaKirkInp
    real(dp), pointer :: polar(:) => null()  !* Atomic polarisabilities (nAtom)
    real(dp), pointer :: rWaals(:) => null()  !* Van der Waals radii (nAtom)
    real(dp), pointer :: charges(:) => null()  !* Effective charges (nAtom)
  end type TDispSlaKirkInp


  !!* Data for the Slater-Kirkwood type dispersion
  type ODispSlaKirk
    private
    real(dp), pointer :: c6(:,:)      !* Atomic polarisabilities (nAtom)
    real(dp), pointer :: rVdW2(:,:)   !* Van der Waals radii (nAtom)
    integer :: nAtom                  !* Nr. of atoms (without images)
    real(dp), pointer :: energies(:)  !* Energies
    real(dp), pointer :: gradients(:,:) !* Gradients (3, nAtom)
    real(dp) :: stress(3,3)           !* stress tensor components
    logical :: tPeriodic              !* If system is periodic
    real(dp) :: rCutoff               !* Real space cutoff
    real(dp) :: gCutoff               !* Reciprocal space cutoff
    real(dp) :: dampCutoff            !* Cutoff, where damping function = 1
    real(dp) :: eta                   !* Periodic summation parameter
    real(dp) :: vol                   !* Volume of the unit cell
    real(dp) :: maxR                  
    real(dp), pointer :: gLatPoint(:,:)  !* Temporary dirty solution
    logical :: coordsUpdated          !* If first coordinate update done
    logical :: tInit = .false.
  end type ODispSlaKirk


  !!* Consctructors
  interface init
    module procedure SlaterKirkwood_init
  end interface

  !!* Destructors
  interface destruct
    module procedure SlaterKirkwood_destruct
    module procedure SlaterKirkwoodInp_destruct
  end interface

  !!* Updates the coordinates of the atoms
  interface updateCoords
    module procedure SlaterKirkwood_updateCoords
  end interface

  !!* Update the lattice vectors
  interface updateLatVecs
    module procedure SlaterKirkwood_updateLatVecs
  end interface

  !!* Gets the atomic resolved energies due to the dispersion
  interface getEnergies
    module procedure SlaterKirkwood_getEnergies
  end interface

  !!* Adds the atomic resoved gradients due to the dispersion
  interface addGradients
    module procedure SlaterKirkwood_addGradients
  end interface

  !!* Returns the real space cutoff for the dispersion interaction
  interface getRCutoff
    module procedure SlaterKirkwood_getRCutoff
  end interface

  !!* Returns the reciprocal space cutoff for the dispersion interaction
  !interface getGCutoff
  !  module procedure SlaterKirkwood_getGCutoff
  !end interface

  !! Some magic constants for the damping function (see paper)
  integer, parameter :: nn_ = 7         ! N
  integer, parameter :: mm_ = 4         ! M
  real(dp), parameter :: dd_ = 3.0_dp   ! d

  
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Public routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!* Initializes a SlaterKirkwood instance.
  !!* @param sf Initialised instance at return.
  !!* @param inp Input parameters for Slater-Kirkwood.
  !!* @param latVecs Lattice vectors, if the system is periodic
  !!* @param recVecs Reciprocal vectors, if the system is periodic
  !!* @param vol Volume of the unit cell, if the system is periodic
  subroutine SlaterKirkwood_init(sf, inp, latVecs, recVecs, vol)
    type(ODispSlaKirk), intent(inout) :: sf
    type(TDispSlaKirkInp), intent(in) :: inp
    real(dp), intent(in), optional :: latVecs(:,:)
    real(dp), intent(in), optional :: recVecs(:,:)
    real(dp), intent(in), optional :: vol

    integer :: iAt1, iAt2
    real(dp) :: tol, rTmp, c6sum

    ASSERT(.not. sf%tInit)
    ASSERT(size(inp%polar) > 0)
    ASSERT(size(inp%polar) == size(inp%rWaals))
    ASSERT(size(inp%polar) == size(inp%charges))
    ASSERT(all(inp%polar >= 0.0_dp))
    ASSERT(all(inp%rWaals >= 0.0_dp))
    ASSERT(present(latVecs) .eqv. present(recVecs))
    ASSERT(present(latVecs) .eqv. present(vol))
    ASSERT_ENV(if (present(latVecs)) then)
    ASSERT_ENV(ASSERT(all(shape(latVecs) == (/ 3, 3 /))))
    ASSERT_ENV(ASSERT(all(shape(recVecs) == (/ 3, 3 /))))
    ASSERT_ENV(ASSERT(vol > 0.0_dp))
    ASSERT_ENV(end if)


    sf%nAtom = size(inp%polar)
    INITALLOCATE_PARR(sf%c6, (sf%nAtom, sf%nAtom))
    INITALLOCATE_PARR(sf%rVdW2, (sf%nAtom, sf%nAtom))
    sf%rCutoff = 0.0_dp
    sf%c6 = 0.0_dp
    sf%rVdW2 = 0.0_dp
    sf%maxR = 0.0_dp
    tol = epsilon(1.0_dp)
    do iAt1 = 1, sf%nAtom
      if (inp%polar(iAt1) < tol .or. inp%rWaals(iAt1) < tol) then
        cycle
      end if
      do iAt2 = 1, iAt1
        sf%c6(iAt2, iAt1) = 1.5_dp * inp%polar(iAt1) * inp%polar(iAt2)&
            &/ (sqrt(inp%polar(iAt1)/inp%charges(iAt1)) &
            &+ sqrt(inp%polar(iAt2)/inp%charges(iAt2)))
        rTmp = (inp%rWaals(iAt1)**3 + inp%rWaals(iAt2)**3) &
            &/ (inp%rWaals(iAt1)**2 + inp%rWaals(iAt2)**2)
        sf%rVdW2(iAt2, iAt1) = dd_ / rTmp**nn_
        sf%maxR = max(sf%maxR, rTmp)
        if (iAt1 /= iAt2) then
          sf%c6(iAt1, iAt2) = sf%c6(iAt2, iAt1)
          sf%rVdW2(iAt1, iAt2) = sf%rVdW2(iAt2, iAt1)
        end if
      end do
    end do
    sf%rCutoff = (maxval(sf%c6)/tolDispersion)**(1.0_dp/6.0_dp)

    sf%tPeriodic = present(latVecs)
    if (sf%tPeriodic) then
      sf%vol = vol
      !! Scaling down optimal eta (as suggested in the literature) is purely
      !! empirical, it reduces the real space summation, and seems to yield
      !! shorter execution times. (It doesn't influence the result.)
      sf%eta =  getOptimalEta(latVecs, sf%vol) / sqrt(2.0_dp)
      c6sum = sum(abs(sf%c6))
      sf%rCutoff = getMaxRDispersion(sf%eta, c6sum, sf%vol, &
          &tolDispersion)
      !! Cutoff, beyond which dispersion is purely 1/r^6 without damping
      sf%dampCutoff = getDampCutoff_(sf%maxR, tolDispDamp)
      sf%rCutoff = max(sf%rCutoff, sf%dampCutoff)
      sf%gCutoff = getMaxGDispersion(sf%eta, c6sum, tolDispersion)
      INIT_PARR(sf%gLatPoint)
      call getLatticePoints(sf%gLatPoint, recVecs, latVecs/(2.0_dp*pi), &
          &sf%gCutoff, onlyInside=.true., reduceByInversion=.true., &
          &withoutOrigin=.true.)
      sf%gLatPoint = matmul(recVecs, sf%gLatPoint)
    end if

    INITALLOCATE_PARR(sf%energies, (sf%nAtom))
    INITALLOCATE_PARR(sf%gradients, (3, sf%nAtom))
    sf%coordsUpdated = .false.
    sf%tInit = .true.

  end subroutine SlaterKirkwood_init


  !!* Destructs SlaterKirkwood instance.
  !!* @param sf SK-instance.
  subroutine SlaterKirkwood_destruct(sf)
    type(ODispSlaKirk), pointer :: sf

    ASSERT(sf%tInit)
    DEALLOCATE_PARR(sf%c6)
    DEALLOCATE_PARR(sf%rVdW2)
    DEALLOCATE_PARR(sf%energies)
    DEALLOCATE_PARR(sf%gradients)
    sf%tInit = .false.

  end subroutine SlaterKirkwood_destruct


  !!* Notifies the objects about changed coordinates
  !!* @param sf Object instance.
  !!* @param neigh Current neighbor list mapping
  !!* @param img2CentCell Current mapping of periodic images to the central cell
  !!* @param coords Current coordinates of the atoms
  subroutine SlaterKirkwood_updateCoords(sf, neigh, img2CentCell, coords)
    type(ODispSlaKirk), intent(inout) :: sf
    type(TNeighborList), intent(in) :: neigh
    integer, intent(in) :: img2CentCell(:)
    real(dp), intent(in) :: coords(:,:)

    integer, allocatable :: nNeighReal(:) ! Neighbors for real space summation
    integer, allocatable :: nNeighDamp(:) ! Nr. of neighbors with damping

    ASSERT(sf%tInit)
    
    ALLOCATE_(nNeighReal, (sf%nAtom))
    call getNrOfNeighborsForAll(nNeighReal, neigh, sf%rCutoff)
    sf%energies(:) = 0.0_dp
    sf%gradients(:,:) = 0.0_dp
    if (sf%tPeriodic) then
      !! Make Ewald summation for a pure 1/r^6 interaction
      call addDispEGr_per_atom(sf%nAtom, coords, nNeighReal, &
          &neigh%iNeighbor, neigh%neighDist2, img2CentCell, sf%c6, sf%eta, &
          &sf%vol, sf%gLatPoint, sf%energies, sf%gradients,sf%stress)
      !! Correct those terms, where damping is important
      ALLOCATE_(nNeighDamp, (sf%nAtom))
      call getNrOfNeighborsForAll(nNeighDamp, neigh, sf%dampCutoff)
      call addDispEnergyAndGrad_cluster_(sf%nAtom, coords, nNeighDamp, &
          &neigh%iNeighbor, neigh%neighDist2, img2CentCell, sf%c6, sf%rVdW2, &
          &sf%energies, sf%gradients, dampCorrection=-1.0_dp)
    else
      call addDispEnergyAndGrad_cluster_(sf%nAtom, coords, nNeighReal, &
          &neigh%iNeighbor, neigh%neighDist2, img2CentCell, sf%c6, sf%rVdW2, &
          &sf%energies, sf%gradients)
    end if
    sf%coordsUpdated = .true.

    DEALLOCATE_(nNeighReal)
    DEALLOCATE_(nNeighDamp)
    
  end subroutine SlaterKirkwood_updateCoords



  !!* Notifies object about updated lattice vectors
  !!* @param latVecs  New lattice vectors
  !!* @param recVecs  New reciprocal vectors
  !!* @param vol  New unit cell volume
  subroutine SlaterKirkwood_updateLatVecs(sf, latVecs, recVecs, vol)
    type(ODispSlaKirk), pointer :: sf
    real(dp), intent(in) :: latVecs(:,:), recVecs(:,:)
    real(dp), intent(in) :: vol
    
    real(dp) :: c6sum

    sf%vol = vol
    sf%eta =  getOptimalEta(latVecs, sf%vol) / sqrt(2.0_dp)
    c6sum = sum(abs(sf%c6))
    sf%rCutoff = getMaxRDispersion(sf%eta, c6sum, sf%vol, tolDispersion)
    !! Cutoff, beyond which dispersion is purely 1/r^6 without damping
    sf%dampCutoff = getDampCutoff_(sf%maxR, tolDispDamp)
    sf%rCutoff = max(sf%rCutoff, sf%dampCutoff)
    sf%gCutoff = getMaxGDispersion(sf%eta, c6sum, tolDispersion)
    DEALLOCATE_PARR(sf%gLatPoint)
    call getLatticePoints(sf%gLatPoint, recVecs, latVecs/(2.0_dp*pi), &
        &sf%gCutoff, onlyInside=.true., reduceByInversion=.true., &
        &withoutOrigin=.true.)
    sf%gLatPoint = matmul(recVecs, sf%gLatPoint)
    sf%coordsUpdated = .false.

  end subroutine SlaterKirkwood_updateLatVecs
  

  !!* Returns the atomic resolved energies due to the dispersion.
  !!* @param sf Object instance
  !!* @param energies Contains the atomic energy contributions on exit.
  subroutine SlaterKirkwood_getEnergies(sf, energies)
    type(ODispSlaKirk), intent(in) :: sf
    real(dp), intent(out) :: energies(:)

    ASSERT(sf%tInit)
    ASSERT(sf%coordsUpdated)
    ASSERT(size(energies) == sf%nAtom)
    
    energies(:) = sf%energies(:)
    
  end subroutine SlaterKirkwood_getEnergies
  

  !!* Adds the atomic gradients to the provided vector.
  !!* @param sf Object instance.
  !!* @param gradients The vector to increase by the gradients.
  subroutine SlaterKirkwood_addGradients(sf, gradients)
    type(ODispSlaKirk), intent(in) :: sf
    real(dp), intent(inout) :: gradients(:,:)

    ASSERT(sf%tInit)
    ASSERT(sf%coordsUpdated)
    ASSERT(all(shape(gradients) == (/ 3, sf%nAtom /)))
    
    gradients(:,:) = gradients(:,:) + sf%gradients(:,:)
    
  end subroutine SlaterKirkwood_addGradients


  !!* Estimates the real space cutoff of the dispersion interaction.
  !!* @param sf Object instance
  !!* @return Cutoff
  function SlaterKirkwood_getRCutoff(sf) result(cutoff)
    type(ODispSlaKirk), intent(in) :: sf
    real(dp) :: cutoff

    cutoff = sf%rCutoff

  end function SlaterKirkwood_getRCutoff

  
  !!* Estimates the reciprocal space cutoff of the dispersion interaction.
  !!* @param sf Object instance
  !!* @return Cutoff
  !function SlaterKirkwood_getGCutoff(sf) result(cutoff)
  !  type(ODispSlaKirk), pointer :: sf
  !  real(dp) :: cutoff
  !
  !  cutoff = sf%gCutoff
  !
  !end function SlaterKirkwood_getGCutoff

  
  !!* Destructs the object instance.
  !!* @param sf Object instance.
  subroutine SlaterKirkwoodInp_destruct(sf)
    type(TDispSlaKirkInp), intent(inout) :: sf

    DEALLOCATE_PARR(sf%polar)
    DEALLOCATE_PARR(sf%rWaals)
    DEALLOCATE_PARR(sf%charges)
    
  end subroutine SlaterKirkwoodInp_destruct

  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Private routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !!* Adds the energy per atom and the gradients for the cluster case
  !!* @param nAtom Nr. of atoms (without periodic images)
  !!* @param coords Coordinates of the atoms (including images)
  !!* @param nNeighbors Nr. of neighbors for each atom
  !!* @param iNeighbor Neighborlist.
  !!* @param neighDist2 Square distances of the neighbours.
  !!* @param img2CentCell Mapping into the central cell.
  !!* @param c6 Van der Waals coefficients (nAtom, nAtom)
  !!* @param rVdW2 Scaled inverse van der Waals radii (nAtom, nAtom)
  !!* @param energies Updated energy vector at return
  !!* @param gradients Updated gradient vector at return
  !!* @param dampCorrection Adds the provided value to the damping function
  !!*   (use -1.0 to sum up damped 1/r^6 terms and subtract pure 1/r^6 ones, in
  !!*   order to correct periodic Ewald sum for the short range damped terms.)
  subroutine addDispEnergyAndGrad_cluster_(nAtom, coords, nNeighbors, &
      &iNeighbor, neighDist2, img2CentCell, c6, rVdW2, energies, gradients, &
      &dampCorrection)
    integer, intent(in) :: nAtom
    real(dp), intent(in) :: coords(:,:)
    integer, intent(in) :: nNeighbors(:)
    integer, intent(in) :: iNeighbor(0:,:)
    real(dp), intent(in) :: neighDist2(0:,:)
    integer, intent(in) :: img2CentCell(:)
    real(dp), intent(in) :: c6(:,:)
    real(dp), intent(in) :: rVdW2(:,:)
    real(dp), intent(inout) :: energies(:)
    real(dp), intent(inout) :: gradients(:,:)
    real(dp), intent(in), optional :: dampCorrection

    integer :: iAt1, iNeigh, iAt2, iAt2f
    real(dp) :: dist2, dist, h0, h1, h2, rTmp
    real(dp) :: diff(3), gr(3)
    real(dp) :: corr

    if (present(dampCorrection)) then
      corr = dampCorrection
    else
      corr = 0.0_dp
    end if

    !! Cluster case => explicit sum of the contributions
    !! NOTE: the cluster summation also (ab)used in the periodic case, neighbors
    !! may go over the cell boundary -> img2CentCell needed for folding back.
    do iAt1 = 1, nAtom
      do iNeigh = 1, nNeighbors(iAt1)
        iAt2 = iNeighbor(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)
        if (c6(iAt2f, iAt1) == 0.0_dp) then
          cycle
        end if
        dist2 = neighDist2(iNeigh, iAt1)
        if (dist2 > minNeighDist2) then
          dist = sqrt(dist2)
          h0 = rVdW2(iAt2f, iAt1)
          h1 = exp(-1.0_dp * h0 * dist**nn_)
          h2 = 1.0_dp - h1
          !! Energy
          rTmp = -0.5_dp * c6(iAt2f, iAt1) * (h2**mm_ + corr) / dist**6
          energies(iAt1) = energies(iAt1) + rTmp
          if (iAt1 /= iAt2f) then
            energies(iAt2f) = energies(iAt2f) + rTmp
          end if
          !! Gradients
          diff(:) = (coords(:,iAt1) - coords(:,iAt2))
          gr(:) = -c6(iAt2f, iAt1) * diff(:) &
              &* (mm_*h2**(mm_-1)*h1*h0*nn_*dist**(nn_-8) &
              &- 6.0_dp * (h2**mm_ + corr) * dist**(-8))
          gradients(:,iAt1) = gradients(:,iAt1) + gr(:)
          gradients(:,iAt2f) = gradients(:,iAt2f) - gr(:)
        end if
      end do
    end do
    
  end subroutine addDispEnergyAndGrad_cluster_


  !!* Returns the distance, beyond that the damping function equals approx. 1.
  !!* @param r0 Length scaling parameter
  !!* @param tol  Tolerance value.
  !!* @return cutoff
  function getDampCutoff_(r0, tol) result(xx)
    real(dp), intent(in) :: r0, tol
    real(dp) :: xx

    !! solve: 1 - tol < (1-exp(-d*(r/r0)^N))^M  for r and hope that the 
    !! logarithm is not blowing up your computer.
    xx = r0 * (-1.0_dp/dd_ * log(1.0_dp &
        &- (1.0_dp - tol)**(1.0_dp/real(mm_,dp))))**(1.0_dp/real(nn_, dp))

  end function getDampCutoff_


end module DispSlaterKirkwood
