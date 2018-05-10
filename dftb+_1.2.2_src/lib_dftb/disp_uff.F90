!!* Dispersion a la UFF, similar to Thomas Heines approach in the deMon code.
!!* @ref L. Zheckov et al., JCTC 1, 841-847 (2005)
!!* @note Periodic case could be inaccurate, if two atoms are very close to
!!*   each other.
!!* @todo Take the reciprocal lattice vectors from outside.
module DispUff
#include "allocate.h"
#include "assert.h"
  use Accuracy
  use Periodic, only: TNeighborList, getNrOfNeighborsForAll, getLatticePoints
  use Constants, only: pi
  use DispCommon
  implicit none
  private

  public :: TDispUFFInp, ODispUFF
  public :: init, destruct, updateCoords, updateLatVecs
  public :: getEnergies, addGradients, getStress
  public :: getRCutoff

  !!* Input structure for the van der Waals initialization.
  type TDispUFFInp
    real(dp), pointer :: energies(:) => null()  !* potential depths (nSpecie)
    real(dp), pointer :: distances(:) => null() !* van der Waals radii (nSpecie)
  end type TDispUFFInp

  !!* Internal state of the van der Waals dispersion module.
  type ODispUFF
    private
    integer :: nAtom, nSpecie                     ! Nr. of atoms, species
    real(dp), pointer :: c6(:,:) => null()        ! Prefactors for r^-6
    real(dp), pointer :: c12(:,:) => null()       ! Prefactors for r^-12
    real(dp), pointer :: cPoly(:,:,:) => null()   ! Prefactors for polynomial
    real(dp), pointer :: r0(:,:) => null()        ! Switching radius
    real(dp) :: rCutoff                           ! Real space cutoff
    logical :: tPeriodic                          ! Periodic system?
    real(dp) :: vol                               ! Volume of the unit cell
    real(dp) :: eta                               ! Ewald summation parameter
    real(dp) :: ewaldRCut, ewaldGCut              ! Ewald cutoff radii
    real(dp), pointer :: gLatPoints(:,:) => null() ! Reciprocal lattice vectors
    real(dp), pointer :: energies(:) => null()
    real(dp), pointer :: gradients(:,:) => null()
    real(dp) :: stress(3,3) = 0.0_dp              ! stress tensor component
    logical :: coordsUpdated                      ! If first coordinate update
                                                  ! done
    logical :: tInit = .false.
  end type ODispUFF
  
  
  !!* Consctructors.
  interface init
    module procedure DispUFF_init
  end interface

  !!* Destructors.
  interface destruct
    module procedure DispUFF_destruct
    module procedure DispUFFInp_destruct
  end interface

  !!* Updates the coordinates of the atoms.
  interface updateCoords
    module procedure DispUFF_updateCoords
  end interface

  !!* Updates lattice vectors.
  interface updateLatVecs
    module procedure DispUFF_updateLatVecs
  end interface

  !!* Gets the atomic resolved energies due to the dispersion.
  interface getEnergies
    module procedure DispUFF_getEnergies
  end interface

  !!* Adds the atomic resoved gradients due to the dispersion.
  interface addGradients
    module procedure DispUFF_addGradients
  end interface
  
  !!* Returns the stress tensor due to dispersion.
  interface getStress
    module procedure DispUFF_getStress
  end interface
  
  !!* Returns the real space cutoff for the dispersion interaction.
  interface getRCutoff
    module procedure DispUFF_getRCutoff
  end interface

  !!!* Returns the reciprocal space cutoff for the dispersion interaction
  !interface getGCutoff
  !  module procedure DispUFF_getGCutoff
  !end interface

  
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Public routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!* Inits a DispUFF instance.
  !!* @param sf Initialised instance at return.
  !!* @param inp Specific input parameters for Slater-Kirkwood.
  !!* @param nAtom Nr. of atoms in the system.
  !!* @param species0 Species of every atom in the unit cell.
  !!* @param latVecs Lattice vectors, if system is periodic.
  !!* @param recVecs Reciprocal lattice vectors, if system is periodic.
  !!* @param vol Unit cell volume, if system is periodic.
  subroutine DispUFF_init(sf, inp, nAtom, species0, latVecs, recVecs, vol)
    type(ODispUFF), intent(inout) :: sf
    type(TDispUFFInp), intent(in) :: inp
    integer, intent(in) :: nAtom
    integer, intent(in), optional :: species0(:)
    real(dp), intent(in), optional :: latVecs(:,:), recVecs(:,:)
    real(dp), intent(in), optional :: vol

    integer :: iSp1, iSp2, iAt1
    real(dp), allocatable :: dij(:,:), rij(:,:)
    real(dp) :: preU0, preU5, preU10, c6sum

    ASSERT(.not. sf%tInit)
    ASSERT(size(inp%energies) > 0)
    ASSERT(size(inp%distances) == size(inp%energies))
    ASSERT(all(inp%energies >= 0.0_dp))
    ASSERT(all(inp%distances >= 0.0_dp))
    ASSERT(present(latVecs) .eqv. present(species0))
    ASSERT(present(latVecs) .eqv. present(recVecs))
    ASSERT(present(latVecs) .eqv. present(vol))
    ASSERT_ENV(if (present(latVecs)) then)
    ASSERT_ENV(ASSERT(all(shape(latVecs) == (/ 3, 3 /))))
    ASSERT_ENV(ASSERT(all(shape(recVecs) == (/ 3, 3 /))))
    ASSERT_ENV(ASSERT(vol > 0.0_dp))
    ASSERT_ENV(end if)

    sf%nSpecie = size(inp%energies)
    sf%nAtom = nAtom
    INITALLOCATE_PARR(sf%c6, (sf%nSpecie, sf%nSpecie))
    INITALLOCATE_PARR(sf%c12, (sf%nSpecie, sf%nSpecie))
    INITALLOCATE_PARR(sf%cPoly, (3, sf%nSpecie, sf%nSpecie))
    INITALLOCATE_PARR(sf%r0, (sf%nSpecie, sf%nSpecie))
    
    ALLOCATE_(dij, (sf%nSpecie, sf%nSpecie))
    ALLOCATE_(rij, (sf%nSpecie, sf%nSpecie))
    forall(iSp1=1:sf%nSpecie, iSp2=1:sf%nSpecie)
      dij(iSp1,iSp2) = sqrt(inp%energies(iSp1) * inp%energies(iSp2))
      rij(iSp1,iSp2) = sqrt(inp%distances(iSp1) * inp%distances(iSp2))
    end forall
    
    sf%c6 = 2.0_dp * dij * rij**6
    sf%c12 = dij * rij**12

    preU0 = 396.0_dp / 25.0_dp
    preU5 = 2.0_dp**(5.0_dp/6.0_dp) * 672.0_dp / 25.0_dp
    preU10 = -(2.0_dp**(2.0_dp/3.0_dp)) * 552.0_dp / 25.0_dp
    sf%cPoly(1,:,:) = preU0 * dij
    sf%cPoly(2,:,:) = preU5 * dij / rij**5
    sf%cPoly(3,:,:) = preU10 * dij / rij**10
    sf%r0 = 2.0_dp**(-1.0_dp/6.0_dp) * rij

    DEALLOCATE_(rij)
    DEALLOCATE_(dij)
    
    sf%tPeriodic = present(latVecs)
    if (sf%tPeriodic) then
      ! Cutoff for the direct summation of r^(-12) terms. To be sure, it is
      ! delivering the required accuracy, dispTol is strengthened by two orders
      ! of magnitude more.
      sf%rCutoff = (maxval(sf%c12)/(tolDispersion*1.0e-2_dp))**(1.0_dp/12.0_dp)
      ! Summing with loop to avoid creation of (nAtom, nAtom) tmp array.
      c6sum = 0.0_dp
      do iAt1 = 1, nAtom
        c6sum = c6sum + sum(abs(sf%c6(species0,species0(iAt1))))
      end do
      sf%vol = vol
      sf%eta =  getOptimalEta(latVecs, sf%vol) / sqrt(2.0_dp)
      sf%ewaldRCut = getMaxRDispersion(sf%eta, c6sum, sf%vol, tolDispersion)
      sf%ewaldGCut = getMaxGDispersion(sf%eta, c6sum, tolDispersion)
      INIT_PARR(sf%gLatPoints)
      call getLatticePoints(sf%gLatPoints, recVecs, latVecs/(2.0_dp*pi), &
          &sf%ewaldGCut, onlyInside=.true., reduceByInversion=.true., &
          &withoutOrigin=.true.)
      sf%gLatPoints = matmul(recVecs, sf%gLatPoints)
    else
      ! Cutoff for the direct real space summation of r^(-6) terms.
      sf%rCutoff = (maxval(sf%c6) / tolDispersion)**(1.0_dp/6.0_dp)
    end if

    INITALLOCATE_PARR(sf%energies, (sf%nAtom))
    INITALLOCATE_PARR(sf%gradients, (3, sf%nAtom))
    sf%tInit = .true.
    
  end subroutine DispUFF_init


  !!* Destorys DispUFF instance.
  !!* @param sf SK-instance.
  subroutine DispUFF_destruct(sf)
    type(ODispUFF), pointer :: sf

    ASSERT(sf%tInit)
    DEALLOCATE_PARR(sf%c6)
    DEALLOCATE_PARR(sf%c12)
    DEALLOCATE_PARR(sf%cPoly)
    DEALLOCATE_PARR(sf%r0)
    DEALLOCATE_PARR(sf%energies)
    DEALLOCATE_PARR(sf%gradients)
    DEALLOCATE_PARR(sf%gLatPoints)

  end subroutine DispUFF_destruct


  !!* Notifies the objects about changed coordinates
  !!* @param sf Object instance.
  !!* @param neigh Current neighbor list mapping
  !!* @param img2CentCell Images of every atom in the central cell.
  !!* @param coords Current coordinates of the atoms
  !!* @param species0 Species of the atoms in the unit cell.
  subroutine DispUFF_updateCoords(sf, neigh, img2CentCell, coords, species0)
    type(ODispUFF), pointer :: sf
    type(TNeighborList), intent(in) :: neigh
    integer, intent(in) :: img2CentCell(:)
    real(dp), intent(in) :: coords(:,:)
    integer, intent(in) :: species0(:)

    integer, allocatable :: nNeigh(:)

    ASSERT(sf%tInit)

    ALLOCATE_(nNeigh, (sf%nAtom))
    call getNrOfNeighborsForAll(nNeigh, neigh, sf%rCutoff)
    if (sf%tPeriodic) then
      call getDispEnergyAndGrad_cluster_(sf%nAtom, coords, species0, nNeigh, &
          &neigh%iNeighbor, neigh%neighDist2, img2CentCell, sf%c6, sf%c12, &
          &sf%cPoly, sf%r0, sf%energies, sf%gradients, removeR6=.true., &
          &stress=sf%stress, vol=sf%vol)
      call getNrOfNeighborsForAll(nNeigh, neigh, sf%ewaldRCut)
      call addDispEGr_per_species(sf%nAtom, coords, species0, nNeigh, &
          &neigh%iNeighbor, neigh%neighDist2, img2CentCell, &
          &sf%c6, sf%eta, sf%vol, sf%gLatPoints, sf%energies, sf%gradients, &
          &sf%stress)
    else
      call getDispEnergyAndGrad_cluster_(sf%nAtom, coords, species0, nNeigh, &
          &neigh%iNeighbor, neigh%neighDist2, img2CentCell, sf%c6, sf%c12, &
          &sf%cPoly, sf%r0, sf%energies, sf%gradients)
    end if
        
    sf%coordsUpdated = .true.
    DEALLOCATE_(nNeigh)
    
  end subroutine DispUFF_updateCoords


  !!* Notifies the object about updated lattice vectors.
  !!* @param latVecs  New lattice vectors
  !!* @param recVecs  New reciprocal vectors
  !!* @param vol  New unit cell volume
  !!* @param species0  Species in the unit cell.
  subroutine DispUFF_updateLatVecs(sf, latVecs, recVecs, vol, species0)
    type(ODispUFF), pointer :: sf
    real(dp), intent(in) :: latVecs(:,:), recVecs(:,:)
    real(dp), intent(in) :: vol
    integer, intent(in) :: species0(:)

    real(dp) :: c6sum
    integer :: iAt1

    ASSERT(sf%tInit)
    ASSERT(sf%tPeriodic)
    ASSERT(all(shape(latVecs) == (/ 3, 3 /)))
    ASSERT(all(shape(recVecs) == (/ 3, 3 /)))
    ASSERT(vol > 0)
    
    sf%rCutoff = (maxval(sf%c12)/(tolDispersion*1.0e-2_dp))**(1.0_dp/12.0_dp)
    ! Summing with loop to avoid creation of (nAtom, nAtom) tmp array.
    c6sum = 0.0_dp
    do iAt1 = 1, sf%nAtom
      c6sum = c6sum + sum(abs(sf%c6(species0,species0(iAt1))))
    end do
    sf%vol = vol
    sf%eta =  getOptimalEta(latVecs, sf%vol) / sqrt(2.0_dp)
    sf%ewaldRCut = getMaxRDispersion(sf%eta, c6sum, sf%vol, tolDispersion)
    sf%ewaldGCut = getMaxGDispersion(sf%eta, c6sum, tolDispersion)
    DEALLOCATE_PARR(sf%gLatPoints)
    call getLatticePoints(sf%gLatPoints, recVecs, latVecs/(2.0_dp*pi), &
        &sf%ewaldGCut, onlyInside=.true., reduceByInversion=.true., &
        &withoutOrigin=.true.)
    sf%gLatPoints = matmul(recVecs, sf%gLatPoints)

    sf%coordsUpdated = .false.
    
  end subroutine DispUFF_updateLatVecs
  

  !!* Returns the atomic resolved energies due to the dispersion.
  !!* @param sf Object instance
  !!* @param energies Contains the atomic energy contributions on exit.
  subroutine DispUFF_getEnergies(sf, energies)
    type(ODispUFF), pointer :: sf
    real(dp), intent(out) :: energies(:)

    ASSERT(sf%tInit)
    ASSERT(sf%coordsUpdated)
    ASSERT(size(energies) == sf%nAtom)
    
    energies(:) = sf%energies(:)
    
  end subroutine DispUFF_getEnergies
  

  !!* Adds the atomic gradients to the provided vector.
  !!* @param sf Object instance.
  !!* @param gradients The vector to increase by the gradients.
  subroutine DispUFF_addGradients(sf, gradients)
    type(ODispUFF), pointer :: sf
    real(dp), intent(inout) :: gradients(:,:)

    ASSERT(sf%tInit)
    ASSERT(sf%coordsUpdated)
    ASSERT(all(shape(gradients) == (/ 3, sf%nAtom /)))
    
    gradients(:,:) = gradients(:,:) + sf%gradients(:,:)
    
  end subroutine DispUFF_addGradients


  !!* Returns the stress tensor
  !!* @param sf Object instance.
  !!* @param stress tensor from the dispersion
  subroutine DispUFF_getStress(sf, stress)
    type(ODispUFF), pointer :: sf
    real(dp), intent(out) :: stress(3,3)

    ASSERT(sf%tInit)
    ASSERT(sf%coordsUpdated)

    stress = sf%stress

  end subroutine DispUFF_getStress


  !!* Estimates the real space cutoff of the dispersion interaction.
  !!* @param sf Object instance
  !!* @return Cutoff
  function DispUFF_getRCutoff(sf) result(cutoff)
    type(ODispUFF), pointer :: sf
    real(dp) :: cutoff

    ASSERT(sf%tInit)
    cutoff = max(sf%rCutoff, sf%ewaldRCut)

  end function DispUFF_getRCutoff

  
  !!!* Estimates the reciprocal space cutoff of the dispersion interaction.
  !!!* @param sf Object instance
  !!!* @return Cutoff
  !function DispUFF_getGCutoff(sf) result(cutoff)
  !  type(ODispUFF), pointer :: sf
  !  real(dp) :: cutoff
  !
  !  ASSERT(sf%tInit)
  !  cutoff = sf%gCutoff
  !
  !end function DispUFF_getGCutoff

  
  !!* Destructs the object instance.
  !!* @param sf Object instance.
  subroutine DispUFFInp_destruct(sf)
    type(TDispUFFInp), intent(inout) :: sf

    DEALLOCATE_PARR(sf%energies)
    DEALLOCATE_PARR(sf%distances)
    
  end subroutine DispUFFInp_destruct

  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Private routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !!* Returns the energy per atom and the gradients for the cluster case
  !!* @param nAtom Nr. of atoms (without periodic images)
  !!* @param coords Coordinates of the atoms (including images)
  !!* @param species Species of every atom.
  !!* @param nNeighbors Nr. of neighbors for each atom
  !!* @param iNeighbor Neighborlist.
  !!* @param neighDist2 Square distances of the neighbours.
  !!* @param img2CentCell Mapping into the central cell.
  !!* @param c6 Prefactors for the r^-6 potential
  !!* @param c12 Prefactors for the r^-12 potential
  !!* @param cPoly Prefactors for the polynomial part.
  !!* @param r0 Distances where polynomial repr. should change to LJ.
  !!* @param energies Updated energy vector at return
  !!* @param gradients Updated gradient vector at return
  !!* @param removeR6 If yes, the 1/r^6 term is substracted from every 
  !!*   interaction.
  subroutine getDispEnergyAndGrad_cluster_(nAtom, coords, species, nNeighbors, &
      &iNeighbor, neighDist2, img2CentCell, c6, c12, cPoly, r0, energies, &
      &gradients, removeR6, stress, vol)
    integer, intent(in) :: nAtom
    real(dp), intent(in) :: coords(:,:)
    integer, intent(in) :: species(:)
    integer, intent(in) :: nNeighbors(:), iNeighbor(0:,:)
    real(dp), intent(in) :: neighDist2(0:,:)
    integer, intent(in) :: img2CentCell(:)
    real(dp), intent(in) :: c6(:,:), c12(:,:), cPoly(:,:,:), r0(:,:)
    real(dp), intent(out) :: energies(:), gradients(:,:)
    logical , intent(in), optional :: removeR6
    real(dp), intent(out), optional :: stress(:,:)
    real(dp), intent(in), optional :: vol

    integer :: iAt1, iAt2, iAt2f, iSp1, iSp2, iNeigh, ii
    real(dp) :: rr, r2, r5, r6, r10, r12, k1, k2, dE, dGr, u0, u1, u2, f6
    real(dp) :: gr(3), vec(3)

    ASSERT_ENV(if (present(stress)) then)
    ASSERT_ENV(  ASSERT(all(shape(stress) == (/ 3, 3 /))))
    ASSERT_ENV(endif)

    !! Cluster case => explicit sum of the contributions
    if (present(removeR6)) then
      if (removeR6) then
        f6 = 0.0_dp
      else
        f6 = 1.0_dp
      end if
    else
      f6 = 1.0_dp
    end if
    energies = 0.0_dp
    gradients = 0.0_dp
    if (present(stress)) then
      stress = 0.0_dp
    end if
    do iAt1 = 1, nAtom
      iSp1 = species(iAt1)
      do iNeigh = 1, nNeighbors(iAt1)
        iAt2 = iNeighbor(iNeigh, iAt1)
        vec = coords(:,iAt1)-coords(:,iAt2)
        iAt2f = img2CentCell(iAt2)
        iSp2 = species(iAt2f)
        r2 = neighDist2(iNeigh, iAt1)
        rr = sqrt(r2)
        k1 = c6(iSp2, iSp1)
        r6 = r2 * r2 * r2
        if (rr > r0(iSp2, iSp1)) then
          ! Two atoms far enough: Lennard-Jones potential
          r12 = r6 * r6
          k2 = c12(iSp2, iSp1)
          dE = 0.5_dp * (-(k1 * f6) / r6 + k2 / r12)
          dGr = (6.0_dp * k1 * f6 / r6 - 12.0_dp * k2 / r12) / rr
        elseif (rr > minNeighDist) then
          ! Two atoms close: polynomial potential
          r10 = r2**5
          r5 = sqrt(r10)
          u0 = cPoly(1, iSp2, iSp1)
          u1 = cPoly(2, iSp2, iSp1)
          u2 = cPoly(3, iSp2, iSp1)
          dE = 0.5_dp * ((u0 - u1 * r5 - u2 * r10) + (1.0_dp - f6) * k1 / r6)
          dGr = (-5.0_dp * u1 * r5 - 10.0_dp * u2 * r10 &
              &- 6.0_dp * k1 * (1.0_dp - f6) / r6) / rr 
        else
          ! Two atoms at the same position -> forget it
          dE = 0.0_dp
          dGr = 0.0_dp
        end if
        energies(iAt1) = energies(iAt1) + dE
        if (iAt1 /= iAt2f) then
          energies(iAt2f) = energies(iAt2f) + dE
        end if
        gr(:) =  dGr * (coords(:,iAt1) - coords(:,iAt2)) / rr
        gradients(:,iAt1) = gradients(:,iAt1) + gr
        gradients(:,iAt2f) = gradients(:,iAt2f) - gr
        if (present(stress)) then
          if (iAt1 /= iAt2f) then
            do ii = 1, 3
              stress(:,ii) = stress(:,ii) - gr * vec(ii) / vol
            end do
          else
            do ii = 1, 3
              stress(:,ii) = stress(:,ii) - 0.5_dp * gr * vec(ii) / vol
            end do
          end if
        end if
      end do
    end do

  end subroutine getDispEnergyAndGrad_cluster_
  
  !!* Returns the real space part of the R6 stress tensor and all of the R12
  !!* part
  !!* @param nAtom Nr. of atoms (without periodic images)
  !!* @param coords Coordinates of the atoms (including any periodic images)
  !!* @param species Species of every atom.
  !!* @param nNeighbors Nr. of neighbors for each atom
  !!* @param iNeighbor Neighborlist.
  !!* @param neighDist2 Square distances of the neighbours.
  !!* @param img2CentCell Mapping into the central cell.
  !!* @param c12 Prefactors for the r^-12 potential
  !!* @param cPoly Prefactors for the polynomial part.
  !!* @param r0 Distances where polynomial repr. should change to LJ.
  !!* @param vol Cell volume
  !!* @param st Stress tensor
  subroutine addLocalStress_per_species_(nAtom, coords, species, nNeighbors, &
      &iNeighbor, neighDist2, img2CentCell, c6, c12, cPoly, r0, vol, st)
    integer, intent(in) :: nAtom
    real(dp), intent(in) :: coords(:,:)
    integer, intent(in) :: species(:)
    integer, intent(in) :: nNeighbors(:), iNeighbor(0:,:)
    real(dp), intent(in) :: neighDist2(0:,:)
    integer, intent(in) :: img2CentCell(:)
    real(dp), intent(in) :: c6(:,:), c12(:,:), cPoly(:,:,:), r0(:,:)
    real(dp), intent(in) :: vol
    real(dp), intent(inout) :: st(3,3)
    
    integer :: ii, jj, iAt1, iAt2, iAt2f, iSp1, iSp2, iNeigh
    real(dp) :: rr, r2, r5, r6, r10, r12, k1, dGr, u0, u1, u2
    real(dp) :: gr(3), vect(3)
    
    do iAt1 = 1, nAtom
      iSp1 = species(iAt1)
      do iNeigh = 1, nNeighbors(iAt1)
        iAt2 = iNeighbor(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)
        iSp2 = species(iAt2f)
        r2 = neighDist2(iNeigh, iAt1)
        rr = sqrt(r2)
        k1 = c6(iSp2, iSp1)
        r6 = r2 * r2 * r2
        if (rr > r0(iSp2, iSp1)) then
          ! Two atoms far enough: Lennard-Jones potential
          r12 = r6 * r6
          dGr = - 12.0_dp * c12(iSp2, iSp1) / (r12 * rr)
        elseif (rr > minNeighDist) then
          ! Two atoms close: polynomial potential
          r10 = r2**5
          r5 = sqrt(r10)
          u0 = cPoly(1, iSp2, iSp1)
          u1 = cPoly(2, iSp2, iSp1)
          u2 = cPoly(3, iSp2, iSp1)
          dGr = (-5.0_dp * u1 * r5 - 10.0_dp * u2 * r10 &
              &- 6.0_dp * k1 * (1.0_dp) / r6) / rr 
        else
          ! Two atoms at the same position -> forget it
          dGr = 0.0_dp
        end if
        vect = coords(:,iAt1) - coords(:,iAt2)
        gr(:) =  vect * dGr / rr
        
        if (iAt2f /= iAt1) then
          do ii = 1, 3
            do jj = 1, 3
              st(jj,ii) = st(jj,ii) &
                  & + vect(jj)*gr(ii) + gr(jj)*vect(ii)
            end do
          end do
        else
          do ii = 1, 3
            do jj = 1, 3
              st(jj,ii) = st(jj,ii) + 0.5_dp * &
                  & (vect(jj)*gr(ii) + gr(jj)*vect(ii))
            end do
          end do
        end if
      end do
    end do
    
    st = -0.5_dp * st / vol
        
  end subroutine addLocalStress_per_species_

end module DispUff
