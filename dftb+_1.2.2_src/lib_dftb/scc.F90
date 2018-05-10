!!* Functions and local variables for the SCC calculation.
module SCC
#include "allocate.h"
#include "assert.h"
  use accuracy
  use message
  use coulomb
  use short_gamma
  use FileID
  use constants
  use periodic
  use ExternalCharges
  use blasRoutines
  use CommonTypes
  use ChargeConstraints
  use shift
  implicit none

  private

  public :: TSCCInit, init_SCC, destruct_SCC
  public :: updateCoords_SCC, updateLatVecs_SCC, updateCharges_SCC
  public :: getSCCCutoff, getEnergyPerAtom_SCC
  public :: addForceDCSCC, getShiftPerAtom, getShiftPerL, getSCCEwaldPar
  public :: SCC_getOrbitalEquiv, addStressDCSCC

  !!* Data necessary for initialize the SCC module
  type TSCCInit
    type(TOrbitals), pointer :: orb
    real(dp), pointer :: hubbU(:,:) => null()
    logical, pointer :: tDampedShort(:) => null()
    real(dp) :: dampExp = 0.0_dp
    real(dp), pointer :: latVecs(:,:) => null()
    real(dp), pointer :: recVecs(:,:) => null()
    real(dp) :: volume = 0.0_dp
    real(dp), pointer :: extCharges(:,:) => null()
    real(dp), pointer :: blurWidths(:) => null()
    real(dp), pointer :: chrgConstraints(:,:) => null()
    real(dp), pointer :: thirdOrderOn(:,:) => null()
    real(dp) :: ewaldAlpha = 0.0_dp    ! if > 0 -> manual setting for alpha
  end type TSCCInit
  
  real(dp), parameter :: tolEwald = 1.0e-9_dp !* tolerance for Ewald

  
  !! Private module variables (suffixed with "_" for clarity)
  logical :: tInitialised_ = .false.      ! If module is initialised
  logical :: tCoordUp_                    ! If coordinates updated at least once
  logical :: tChrgUp_                     ! If charges updated at least once   
  integer :: nAtom_, nSpecie_             ! Nr. of atoms and species
  real(dp), allocatable :: invRMat_(:,:)  ! Stores 1/r between atom pairs
  real(dp), allocatable :: shiftPerAtom_(:)      ! Shift vector per atom
  real(dp), allocatable :: shiftPerL_(:,:)       ! Shift vector per l-shell
  real(dp), allocatable :: shortGamma_(:,:,:,:)  ! Short range interaction
  real(dp), allocatable :: shortCutoff_(:,:,:,:) ! Cutoff for short range int.
  real(dp) :: cutoff_                            ! Maximal cutoff
  real(dp), pointer :: gLatPoint_(:,:)     ! Lattice points for reciprocal Ewald
  integer, allocatable :: nNeighShort_(:,:,:,:) ! Nr. of neighbors for short
                                                ! range interaction
  integer, allocatable :: nNeighEwald_(:)       ! Nr. of neigh for real Ewald
  real(dp) :: maxREwald_                        ! Cutoff for real Ewald
  real(dp) :: alpha_                            ! Parameter for Ewald
  real(dp) :: volume_                           ! Cell volume
  real(dp), allocatable :: uniqHubbU_(:,:)      ! Uniq Us per species
  integer, allocatable :: nHubbU_(:)            ! Nr. of uniq Us per species
  integer, allocatable :: iHubbU_(:,:)          ! Mapping L-shell -> uniq U
  logical :: tExtChrg_                          ! Are external charges present?

  real(dp), allocatable :: deltaQ_(:,:)          ! Negative net charge
  real(dp), allocatable :: deltaQPerLShell_(:,:) ! Negative net charge per shell
  real(dp), allocatable :: deltaQAtom_(:)        ! Negative net charge per atom
  real(dp), allocatable :: deltaQUniqU_(:,:)     ! Negative net charge per U

  logical, allocatable :: tDampedShort_(:)      ! Damped short range? (nSpecie)
  real(dp) :: dampExp_                          ! Damping exponent

  logical :: tPeriodic_                          ! Is the system periodic?

  logical :: tChrgConstr_                       ! Shifts due charge constrains?
  type(OChrgConstr), save :: chrgConstr_        ! Object for charge constraints
  logical :: tThirdOrder_                       ! Shifts due to 3rd order
  type(OChrgConstr), save :: thirdOrder_
  logical :: tAutoEwald_
  

contains

  !!* Initialises the SCC module
  !!* @param nAtom Number of atoms in the system
  !!* @param nSpecie Number of species in the system
  !!* @param nShell Nr. of shells for each species
  !!* @param mShell Maximal nr. of shells for any species
  !!* @param mOrb   Maximal nr. of orbitals for any species
  !!* @param hubbU The Hubbard Us read from the SK files
  !!* @param tDampedShort True for species which needs damped interaction
  !!* @param latVec Lattice vectors (if the system is periodic)
  !!* @param recVec Reciprocal vectors (if the system is periodic)
  !!* @param volume Volume of the supercell (if the system is periodic)
  !!* @param extCharges Coordinate and charge of external point charges.
  !!* @param blurWidths widths of Gaussian blurring on charges (non-periodic
  !!* only)  
  subroutine init_SCC(inp)
    type(TSCCInit), intent(inout) :: inp

    integer :: mShell
    integer :: iSp1, iSp2, iU1, iU2, iL
    real(dp) :: maxGEwald

    nSpecie_ = size(inp%orb%nOrbSpecie)
    nAtom_ = size(inp%orb%nOrbAtom)
    mShell = inp%orb%mShell

    ASSERT(.not. tInitialised_)
    ASSERT(associated(inp%latVecs) .eqv. associated(inp%recVecs))
    ASSERT(associated(inp%latVecs) .eqv. (inp%volume > 0.0_dp))
    ASSERT(size(inp%hubbU, dim=1) == mShell)
    ASSERT(size(inp%hubbU, dim=2) == nSpecie_)
    ASSERT(size(inp%tDampedShort) == nSpecie_)
    ASSERT(associated(inp%extCharges) .or. .not. associated(inp%blurWidths))
    ASSERT_ENV(if (associated(inp%extCharges)) then)
    ASSERT_ENV(  ASSERT(size(inp%extCharges, dim=1) == 4))
    ASSERT_ENV(  ASSERT(size(inp%extCharges, dim=2) > 0))
    ASSERT_ENV(end if)
    
    ALLOCATE_(invRMat_, (nAtom_, nAtom_))
    ALLOCATE_(shiftPerAtom_, (nAtom_))
    ALLOCATE_(shiftPerL_, (mShell, nAtom_))
    ALLOCATE_(shortGamma_, (0, 0, 0, 0))
    INIT_PARR(gLatPoint_)
    tPeriodic_ = associated(inp%latVecs)
    tExtChrg_ = associated(inp%extCharges)
    
    !! Initialize Hubbard U's
    ALLOCATE_(uniqHubbU_, (mShell, nSpecie_))
    ALLOCATE_(nHubbU_, (nSpecie_))
    ALLOCATE_(iHubbU_, (mShell, nSpecie_))
    iHubbU_(:,:) = 0
    iHubbU_(1,:) = 1
    nHubbU_(:) = 1
    uniqHubbU_(:,:) = 0.0_dp
    uniqHubbU_(1,:) = inp%hubbU(1,:)
    do iSp1 = 1, nSpecie_
      do iL = 2, inp%orb%nShell(iSp1)
        do iU1 = 1, nHubbU_(iSp1)
          if (abs(inp%hubbU(iL,iSp1) - uniqHubbU_(iU1,iSp1)) < MinHubDiff) then
            iHubbU_(iL,iSp1) = iU1
            exit
          end if
        end do
        if (iHubbU_(iL,iSp1) == 0) then
          nHubbU_(iSp1) = nHubbU_(iSp1) + 1
          uniqHubbU_(nHubbU_(iSp1),iSp1) = inp%hubbU(iL,iSp1)
          iHubbU_(iL,iSp1) = nHubbU_(iSp1)
        end if
      end do
    end do
    
    !! Get cutoff for short range coulomb
    ALLOCATE_(shortCutoff_, (maxval(nHubbU_),maxval(nHubbU_),nSpecie_,nSpecie_))
    shortCutoff_(:,:,:,:) = 0.0_dp
    do iSp1 = 1, nSpecie_
      do iSp2 = iSp1, nSpecie_
        do iU1 = 1, nHubbU_(iSp1)
          do iU2 = 1, nHubbU_(iSp2)
            shortCutoff_(iU2, iU1, iSp2, iSp1) = &
                & expGammaCutoff(uniqHubbU_(iU2, iSp2), uniqHubbU_(iU1, iSp1))
            shortCutoff_(iU1, iU2, iSp1, iSp2) = &
                & shortCutoff_(iU2, iU1, iSp2, iSp1)
          end do
        end do
      end do
    end do
    cutoff_ = maxval(shortCutoff_)

    !! Initialize Ewald summation for the periodic case
    if (tPeriodic_) then
      volume_ = inp%volume
      tAutoEwald_ = inp%ewaldAlpha <= 0.0_dp
      if (tAutoEwald_) then
        alpha_ = getOptimalAlphaEwald(inp%latVecs, inp%recVecs, volume_, &
            &tolEwald)
      else
        alpha_ = inp%ewaldAlpha
      end if
      maxREwald_ = getMaxREwald(alpha_, tolEwald)
      maxGEwald = getMaxGEwald(alpha_, volume_, tolEwald)
      call getLatticePoints(gLatPoint_, inp%recVecs, inp%latVecs/(2.0_dp*pi), &
          &maxGEwald, onlyInside=.true., reduceByInversion=.true., &
          &withoutOrigin=.true.)
      gLatPoint_ = matmul(inp%recVecs, gLatPoint_)
      cutoff_ = max(cutoff_, maxREwald_)
    end if

    !! Number of neighbors for short range cutoff and real part of Ewald
    ALLOCATE_(nNeighShort_, (maxval(nHubbU_), maxval(nHubbU_),nSpecie_,nAtom_))
    if (tPeriodic_) then
      ALLOCATE_(nNeighEwald_, (nAtom_))
    end if

    !! Initialise external charges
    if (tExtChrg_) then
      if (tPeriodic_) then
        call init_ExtChrg(inp%extCharges, nAtom_, inp%latVecs, inp%recVecs, &
            &maxREwald_)
      else
        ASSERT(associated(inp%blurWidths))
        call init_ExtChrg(inp%extCharges, nAtom_, blurWidths=inp%blurWidths)
      end if
    end if

    tChrgConstr_ = associated(inp%chrgConstraints)
    if (tChrgConstr_) then
      call init(chrgConstr_, inp%chrgConstraints, 2)
    end if
    tThirdOrder_ = associated(inp%thirdOrderOn)
    if (tThirdOrder_) then
      ! Factor 1/6 in the energy is put into the Hubbard derivatives
      call init(thirdOrder_, inp%thirdOrderOn / 6.0_dp, 3)
    end if
      
    !! Initialise arrays for charge differences
    ALLOCATE_(deltaQ_, (inp%orb%mOrb, nAtom_))
    ALLOCATE_(deltaQPerLShell_, (mShell, nAtom_))
    ALLOCATE_(deltaQAtom_, (nAtom_))
    ALLOCATE_(deltaQUniqU_, (maxval(nHubbU_), nAtom_))
    
    !! Initialise short range damping
    ALLOCATE_(tDampedShort_, (nSpecie_))
    tDampedShort_(:) = inp%tDampedShort(:)
    dampExp_ = inp%dampExp
    
    tInitialised_ = .true.
    tCoordUp_ = .false.
    tChrgUp_ = .false.

  end subroutine init_SCC

  
  
  !!* Deallocates the arrays in the SCC module
  subroutine destruct_SCC()

    ASSERT(tInitialised_)
    
    DEALLOCATE_(invRMat_)
    DEALLOCATE_(shiftPerAtom_)
    DEALLOCATE_(shiftPerL_)
    DEALLOCATE_(shortGamma_)
    DEALLOCATE_(shortCutoff_)
    DEALLOCATE_PARR(gLatPoint_)
    DEALLOCATE_(nNeighShort_)
    DEALLOCATE_(nNeighEwald_)
    DEALLOCATE_(uniqHubbU_)
    DEALLOCATE_(nHubbU_)
    DEALLOCATE_(iHubbU_)
    DEALLOCATE_(deltaQ_)
    DEALLOCATE_(deltaQPerLShell_)
    DEALLOCATE_(deltaQAtom_)
    DEALLOCATE_(deltaQUniqU_)
    DEALLOCATE_(tDampedShort_)
    if (tChrgConstr_) then
      call destruct(chrgConstr_)
    end if
    if (tThirdOrder_) then
      call destruct(thirdOrder_)
    end if
  
  end subroutine destruct_SCC



  !!* Returns a minimal cutoff for the neighborlist, which must be passed
  !!* to various functions in this module.
  !!* @return cutoff The neighborlists, passed to scc routines, should contain
  !!*   neigbhor information at least up to that cutoff.
  function getSCCCutoff() result(cutoff)
    real(dp) :: cutoff

    ASSERT(tInitialised_)
    cutoff = cutoff_
    
  end function getSCCCutoff


  !!* Returns the currenty used alpha parameter of the Ewald-summation
  !!* @return Parameter in the Ewald summation.
  function getSCCEwaldPar() result(alpha)
    real(dp) :: alpha

    ASSERT(tInitialised_)
    alpha = alpha_

  end function getSCCEwaldPar
  

  !!* Updates the number of neighbors for the SCC module (local).
  !!* @param specie Species for each atom
  !!* @param neighList Neighbor list for the atoms in the system.
  subroutine updateNNeigh_(specie, neighList)
    integer, intent(in) :: specie(:)
    type(TNeighborList), intent(in) :: neighList

    integer :: iAt1, iSp2, iU1, iU2
    
    nNeighShort_(:,:,:,:) = 0
    do iAt1 = 1, nAtom_
      do iSp2 = 1, nSpecie_
        do iU1 = 1, nHubbU_(specie(iAt1))
          do iU2 = 1, nHubbU_(iSp2)
            nNeighShort_(iU2, iU1, iSp2, iAt1) = &
                & getNrOfNeighbors(neighList, &
                & shortCutoff_(iU2, iU1, iSp2, specie(iAt1)), iAt1)
          end do
        end do
      end do
    end do
    if (tPeriodic_) then
      call getNrOfNeighborsForAll(nNeighEwald_, neighList, maxREwald_)
    end if
    
  end subroutine updateNNeigh_
  


  !!* Updates the atom coordinates for the SCC module.
  !!* @param coord New coordinates of the atoms
  !!* @param specie Specie of the atoms (should not change during run)
  !!* @param mAngSpecie Maximal angular momentum of the species (should not
  !!* change during run)
  !!* @param neighList Neighbor list for the atoms.
  !!* @param img2CentCell Mapping to the central cell for the atoms
  subroutine updateCoords_SCC(coord, specie, neighList, img2CentCell)
    real(dp), intent(in) :: coord(:,:)
    integer,  intent(in) :: specie(:)
    type(TNeighborList), intent(in) :: neighList
    integer,  intent(in) :: img2CentCell(:)

    ASSERT(tInitialised_)

    call updateNNeigh_(specie, neighList)
    if (tPeriodic_) then
      call invR(invRMat_, nAtom_, coord, nNeighEwald_, neighList%iNeighbor, &
          &img2CentCell, gLatPoint_, alpha_, volume_)
    else
      call invR(invRMat_, nAtom_, coord)
    end if
    call initGamma_(coord, specie, neighList%iNeighbor)

    if (tExtChrg_) then
      if (tPeriodic_) then
        call updateCoords_ExtChrg(coord, gLatPoint_, alpha_, volume_)
      else
        call updateCoords_ExtChrg(coord)
      end if
    end if

    tCoordUp_ = .true.
    tChrgUp_ = .false.

  end subroutine updateCoords_SCC


  !!* Updates the SCC module, if the lattice vectors had been changed
  !!* @param latVec  New lattice vectors
  !!* @param recVec  New reciprocal lattice vectors
  !!* @param vol  New volume
  subroutine updateLatVecs_SCC(latVec, recVec, vol)
    real(dp), intent(in) :: latVec(:,:), recVec(:,:)
    real(dp), intent(in) :: vol

    real(dp) :: maxGEwald

    ASSERT(tInitialised_)
    ASSERT(tPeriodic_)
    
    volume_ = vol
    if (tAutoEwald_) then
      alpha_ = getOptimalAlphaEwald(latVec, recVec, volume_, tolEwald)
      maxREwald_ = getMaxREwald(alpha_, tolEwald)
    end if
    maxGEwald = getMaxGEwald(alpha_, volume_, tolEwald)
    DEALLOCATE_PARR(gLatPoint_)
    call getLatticePoints(gLatPoint_, recVec, latVec/(2.0_dp*pi), maxGEwald, &
        &onlyInside=.true., reduceByInversion=.true., withoutOrigin=.true.)
    gLatPoint_ = matmul(recVec, gLatPoint_)
    cutoff_ = max(cutoff_, maxREwald_)

    if (tExtChrg_) then
      call updateLatVecs_extChrg(latVec, recVec, maxREwald_)
    end if

    tCoordUp_ = .false.
    tChrgUp_ = .false.

  end subroutine updateLatVecs_SCC


  !!* Up-dates the SCC module, if the charges had been changed
  !!* @param qOrbital Orbital resolved charges
  !!* @param q0 Reference charge distribution (neutral atoms)
  !!* @param mAngSpecie Maximal angular momentum of the species (should not
  !!* change during run)
  !!* @param orb Contains information about the atomic orbitals in the system
  !!* @param specie Specie of the atoms (should not change during run)
  !!* @param iNeighbor Neighbor indexes
  !!* @param img2CentCell Mapping on atoms in the centrall cell
  subroutine updateCharges_SCC(qOrbital, q0, orb, specie, iNeighbor, &
      &img2CentCell)
    real(dp), intent(in) :: qOrbital(:,:,:)
    real(dp), intent(in) :: q0(:,:,:)
    type(TOrbitals), intent(in) :: orb
    integer, intent(in) :: specie(:)
    integer, intent(in) :: iNeighbor(0:,:)
    integer, intent(in) :: img2CentCell(:)

    integer :: iAt1, iSh1, iSp1

    ASSERT(tInitialised_)
    ASSERT(tCoordUp_)

    !! Build charge differences
    deltaQ_(:,:) = qOrbital(:,:,1) - q0(:,:,1)
    deltaQAtom_(:) = sum(deltaQ_(:,:), dim=1)
    deltaQUniqU_(:,:) = 0.0_dp
    deltaQPerLShell_(:,:) = 0.0_dp
    do iAt1 = 1, nAtom_
      iSp1 = specie(iAt1)
      do iSh1 = 1, orb%nShell(iSp1)
        deltaQPerLShell_(iSh1, iAt1) = &
            &sum(deltaQ_(orb%posShell(iSh1,iSp1):orb%posShell(iSh1+1,iSp1)-1,&
            &iAt1))
      end do
      do iSh1 = 1, orb%nShell(iSp1)
        deltaQUniqU_(iHubbU_(iSh1, iSp1), iAt1) = &
            &deltaQUniqU_(iHubbU_(iSh1, iSp1), iAt1) &
            &+ deltaQPerLShell_(iSh1, iAt1)
      end do
    end do

    !! Update charge dependent shift vectors
    call buildShifts_(orb, specie, iNeighbor, img2CentCell)

    tChrgUp_ = .true.

  end subroutine updateCharges_SCC




  !!* Set up the storage and internal values for the short range part of Gamma.
  !!* @param coord        List of coordinates
  !!* @param specie       List of the specie for each atom.
  !!* @param iNeighbor    Index of neighboring atoms for each atom.
  subroutine initGamma_(coord, specie, iNeighbor)
    real(dp), intent(in) :: coord(:,:)
    integer,  intent(in) :: specie(:)
    integer,  intent(in) :: iNeighbor(0:,:)
    
    integer :: iAt1, iAt2, iU1, iU2, iNeigh, iSp1, iSp2
    real(dp) :: rab, u1, u2

    ! Reallocate shortGamma, if it does not contain enough neighbors
    if (size(shortGamma_, dim=3) < maxval(nNeighShort_)+1) then
      DEALLOCATE_(shortGamma_)
      ALLOCATE_(shortGamma_,(maxval(nHubbU_(:)), maxval(nHubbU_(:)), \
          0:maxval(nNeighShort_), nAtom_))
    end if
    shortGamma_(:,:,:,:) = 0.0_dp

    ! some additional symmetry not used, as the value of gamma for atoms
    ! interacting with themselves is the same for all atoms of the same specie
    do iAt1 = 1, nAtom_
      iSp1 = specie(iAt1)
      do iNeigh = 0, maxval(nNeighShort_(:,:,:, iAt1))
        iAt2 = iNeighbor(iNeigh, iAt1)
        iSp2 = specie(iAt2)
        rab = sqrt(sum((coord(:,iAt1) - coord(:,iAt2))**2))
        do iU1 = 1, nHubbU_(specie(iAt1))
          u1 = uniqHubbU_(iU1, iSp1)
          do iU2 = 1, nHubbU_(specie(iAt2))
            u2 = uniqHubbU_(iU2, iSp2)
            if (iNeigh <= nNeighShort_(iU2,iU1,iSp2,iAt1)) then
              if (tDampedShort_(iSp1) .or. tDampedShort_(iSp2)) then
                shortGamma_(iU2 ,iU1, iNeigh, iAt1) = &
                    &expGammaDamped(rab, u2, u1, dampExp_)
              else
                shortGamma_(iU2 ,iU1, iNeigh, iAt1) = expGamma(rab, u2, u1)
              end if
            end if
          end do
        end do
      end do
    end do
    
  end subroutine initGamma_



  !!* Calculate  the derivative of the short range part of Gamma.
  !!* @param force        force vector to add the short-range part of gamma
  !!* contribution  
  !!* @param coord        list of coordinates
  !!* @param specie       List of the specie for each atom.
  !!* @param iNeighbor    Index of neighboring atoms for each atom.
  !!* @param img2CentCell Image of each atom in the central cell.
  subroutine addGammaPrime_(force, coord, specie, iNeighbor, img2CentCell)
    real(dp), intent(inout) :: force(:,:)
    real(dp), intent(in)    :: coord(:,:)
    integer,  intent(in)    :: specie(:)
    integer,  intent(in)    :: iNeighbor(0:,:)
    integer,  intent(in)    :: img2CentCell(:)
    
    integer :: iAt1, iAt2, iAt2f, iU1, iU2, iNeigh, ii, iSp1, iSp2
    real(dp) :: rab, tmpGammaPrime, u1, u2
    
    ASSERT(size(force,dim=1) == 3)
    ASSERT(size(force,dim=2) == nAtom_)
    
    ! some additional symmetry not used
    do iAt1 = 1, nAtom_
      iSp1 = specie(iAt1)
      do iNeigh = 1, maxval(nNeighShort_(:,:,:, iAt1))
        iAt2 = iNeighbor(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)
        iSp2 = specie(iAt2f)
        rab = sqrt(sum((coord(:,iAt1) - coord(:,iAt2))**2))
        do iU1 = 1, nHubbU_(specie(iAt1))
          u1 = uniqHubbU_(iU1, iSp1)
          do iU2 = 1, nHubbU_(specie(iAt2f))
            u2 = uniqHubbU_(iU2, specie(iAt2f))
            if (iNeigh <= nNeighShort_(iU2,iU1,specie(iAt2f),iAt1)) then
              if (tDampedShort_(iSp1) .or. tDampedShort_(iSp2)) then
                tmpGammaPrime = expGammaDampedPrime(rab, u2, u1, dampExp_)
              else
                tmpGammaPrime = expGammaPrime(rab, u2, u1)
              end if
              do ii = 1,3
                force(ii,iAt1) = force(ii,iAt1) - & 
                    & deltaQUniqU_(iU1,iAt1) * &
                    & deltaQUniqU_(iU2,iAt2f) &
                    & *tmpGammaPrime*(coord(ii,iAt1) - coord(ii,iAt2))/rab
                force(ii,iAt2f) = force(ii,iAt2f)&
                    & +deltaQUniqU_(iU1,iAt1) * &
                    & deltaQUniqU_(iU2,iAt2f) &
                    & *tmpGammaPrime*(coord(ii,iAt1) - coord(ii,iAt2))/rab
              end do
            end if
          end do
        end do
      end do
    end do
    
  end subroutine addGammaPrime_


  !!* Calculate  the derivative of the short range part of Gamma.
  !!* @param st           Stress tensor component to add the short-range
  !!*  part of the gamma contribution
  !!* @param coord        list of coordinates
  !!* @param specie       List of the specie for each atom.
  !!* @param iNeighbor    Index of neighboring atoms for each atom.
  !!* @param img2CentCell Image of each atom in the central cell.
  subroutine addSTGammaPrime_(st, coord, specie, iNeighbor, img2CentCell)
    real(dp), intent(out)   :: st(:,:)
    real(dp), intent(in)    :: coord(:,:)
    integer,  intent(in)    :: specie(:)
    integer,  intent(in)    :: iNeighbor(0:,:)
    integer,  intent(in)    :: img2CentCell(:)
    
    integer  :: iAt1, iAt2, iAt2f, iU1, iU2, iNeigh, ii, jj, iSp1, iSp2
    real(dp) :: rab, tmpGammaPrime, u1, u2
    real(dp) :: intermed(3), vect(3)

    ASSERT(all(shape(st)==(/3,3/)))
    
    st(:,:) = 0.0_dp
    ! some additional symmetry not used
    do iAt1 = 1, nAtom_
      iSp1 = specie(iAt1)
      do iNeigh = 1, maxval(nNeighShort_(:,:,:, iAt1))
        iAt2 = iNeighbor(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)
        iSp2 = specie(iAt2f)
        vect(:) = coord(:,iAt1) - coord(:,iAt2)
        rab = sqrt(sum((vect)**2))
        intermed(:) = 0.0_dp
        do iU1 = 1, nHubbU_(specie(iAt1))
          u1 = uniqHubbU_(iU1, iSp1)
          do iU2 = 1, nHubbU_(specie(iAt2f))
            u2 = uniqHubbU_(iU2, specie(iAt2f))
            if (iNeigh <= nNeighShort_(iU2,iU1,specie(iAt2f),iAt1)) then
              if (tDampedShort_(iSp1) .or. tDampedShort_(iSp2)) then
                tmpGammaPrime = expGammaDampedPrime(rab, u2, u1, dampExp_)
              else
                tmpGammaPrime = expGammaPrime(rab, u2, u1)
              end if
              do ii = 1,3
                intermed(ii) = intermed(ii) - & 
                    & deltaQUniqU_(iU1,iAt1) * &
                    & deltaQUniqU_(iU2,iAt2f) &
                    & *tmpGammaPrime*vect(ii)/rab
              end do
            end if
          end do
        end do
        if (iAt2f /= iAt1) then
          do ii = 1, 3
            do jj = 1, 3
              st(jj,ii) = st(jj,ii) &
                  & + (vect(jj)*intermed(ii) + intermed(jj)*vect(ii))
            end do
          end do
        else
          do ii = 1, 3
            do jj = 1, 3
              st(jj,ii) = st(jj,ii) + 0.5_dp * &
                  & (vect(jj)*intermed(ii) + intermed(jj)*vect(ii))
            end do
          end do
        end if
      end do
    end do

    st(:,:) = st(:,:) / volume_
    
  end subroutine addSTGammaPrime_



  !!* Calculates the contribution of the charge consitent part to the energy.
  !!* @param eSCC         The SCC contribution to the energy
  !!* @param specie       Specie for each atom.
  !!* @param iNeighbor    List of neighbors for each atom.
  !!* @param img2CentCell Image of the atoms in the central cell.
  subroutine getEnergy_SCC(eSCC, specie, iNeighbor, img2CentCell)
    real(dp), intent(out) :: eSCC
    integer,  intent(in)  :: specie(:)
    integer,  intent(in)  :: iNeighbor(0:,:)
    integer,  intent(in)  :: img2CentCell(:)

    real(dp), allocatable :: eSCC_atom(:)

    
    ALLOCATE_(eSCC_atom, (nAtom_))

    call getEnergyPerAtom_SCC(eSCC_atom, specie, iNeighbor,img2CentCell)
    
    eSCC = sum(eSCC_atom)

    DEALLOCATE_(eSCC_atom)

  end subroutine getEnergy_SCC
  
  

  !!* Calculates the contribution of the charge consitent part to the energy
  !!* per atom.
  !!* @param eSCC         The SCC contribution to the energy
  !!* @param specie       Specie for each atom.
  !!* @param iNeighbor    List of neighbors for each atom.
  !!* @param img2CentCell Image of the atoms in the central cell.
  subroutine getEnergyPerAtom_SCC(eSCC, specie, iNeighbor, img2CentCell)
    real(dp), intent(out) :: eSCC(:)
    integer,  intent(in)  :: specie(:)
    integer,  intent(in)  :: iNeighbor(0:,:)
    integer,  intent(in)  :: img2CentCell(:)

    integer               :: iAt1, iAt2f, iNeigh, iU1, iU2
    real(dp)              :: rTmp

    ASSERT(tInitialised_)
    ASSERT(size(eSCC) == nAtom_)

    !! 1/R contribution
    eSCC = 0.0_dp
    eSCC = eSCC + 0.5_dp * shiftPerAtom_ * deltaQAtom_
    if (tExtChrg_) then
      call addEnergyPerAtom_ExtChrg(deltaQAtom_, eSCC)
    end if
    if (tChrgConstr_) then
      call addEnergyPerAtom(chrgConstr_, eSCC, deltaQAtom_)
    end if
    if (tThirdOrder_) then
      call addEnergyPerAtom(thirdOrder_, eSCC, deltaQAtom_)
    end if
    
    !! small gamma shifts
    do iAt1 = 1, nAtom_
      do iU1 = 1, nHubbU_(specie(iAt1))
        do iNeigh = 0, maxval(nNeighShort_(:,:,:,iAt1))
          iAt2f = img2CentCell(iNeighbor(iNeigh, iAt1))
          do iU2 = 1, nHubbU_(specie(iAt2f))
            ! should re-write using shift vectors :
            rTmp = deltaQUniqU_(iU2, iAt2f) * deltaQUniqU_(iU1, iAt1) &
                & * shortGamma_(iU2, iU1, iNeigh, iAt1)
            if (iAt1 == iAt2f) then
              eSCC(iAt1) = eSCC(iAt1) - 0.5_dp * rTmp
            else
              eSCC(iAt1) = eSCC(iAt1) - 0.5_dp * rTmp
              eSCC(iAt2f) = eSCC(iAt2f) - 0.5_dp * rTmp
            end if
          end do
        end do
      end do
    end do

  end subroutine getEnergyPerAtom_SCC

  
  !!* Calculates the contribution of the charge consitent part to the forces
  !!* for molecules/clusters, which is not covered in the term with the shift
  !!* vectors.
  !!* @param force        Force contribution
  !!* @param mAngSpecie   Maximal angular momentum per specie.
  !!* @param specie       Specie for each atom.
  !!* @param iNeighbor    List of neighbors for each atom.
  !!* @param img2CentCell Image of the atoms in the central cell.
  !!* @param coord        List of coordinates
  !!* @param chrgForce    Force contribution due to the external charges, which
  !!* is not contained in the term with the shift vectors.
  subroutine addForceDCSCC(force, specie, iNeighbor, img2CentCell, &
      & coord,chrgForce)
    real(dp), intent(inout) :: force(:,:)
    integer,  intent(in)    :: specie(:)
    integer,  intent(in)    :: iNeighbor(0:,:)
    integer,  intent(in)    :: img2CentCell(:)
    real(dp), intent(in)    :: coord(:,:)
    real(dp), intent(inout), optional :: chrgForce(:,:)
    

    ASSERT(size(force,dim=1) == 3)
    ASSERT(size(force,dim=2) == nAtom_)
    ASSERT(present(chrgForce) .eqv. tExtChrg_)
    
    ! Short-range part of gamma contribution
    call addGammaPrime_(force,coord,specie,iNeighbor,img2CentCell)
    
    ! 1/R contribution
    if (tPeriodic_) then
      call addInvRPrime(force, nAtom_, coord, nNeighEwald_, iNeighbor, &
          & img2CentCell, gLatPoint_, alpha_, volume_, deltaQAtom_)
      if (tExtChrg_) then
        call addForceDCSCC_ExtChrg(force, chrgForce, coord, deltaQAtom_, &
            &gLatPoint_, alpha_, volume_)
      end if
    else
      call addInvRPrime(force, nAtom_, coord, deltaQAtom_)
      if (tExtChrg_) then
        call addForceDCSCC_ExtChrg(force, chrgForce, coord, deltaQAtom_)
      end if
    end if

  end subroutine addForceDCSCC

  
  !!* Calculates the contribution of the charge consitent part to the forces
  !!* for molecules/clusters, which is not covered in the term with the shift
  !!* vectors.
  !!* @param force        Force contribution
  !!* @param mAngSpecie   Maximal angular momentum per specie.
  !!* @param specie       Specie for each atom.
  !!* @param iNeighbor    List of neighbors for each atom.
  !!* @param img2CentCell Image of the atoms in the central cell.
  !!* @param coord        List of coordinates
  subroutine addStressDCSCC(st, specie, iNeighbor, img2CentCell, &
      & coord)!,chrgForce)
    real(dp), intent(inout) :: st(:,:)
    integer,  intent(in)    :: specie(:)
    integer,  intent(in)    :: iNeighbor(0:,:)
    integer,  intent(in)    :: img2CentCell(:)
    real(dp), intent(in)    :: coord(:,:)
    
    real(dp) :: stTmp(3,3)
        
    ASSERT(tPeriodic_)
    ASSERT(all(shape(st)==(/3,3/)))
    
    stTmp = 0.0_dp
    
    ! Short-range part of gamma contribution
    call addSTGammaPrime_(stTmp,coord,specie,iNeighbor,img2CentCell)

    st(:,:) = st(:,:) - 0.5_dp * stTmp(:,:)
    
    ! 1/R contribution
    ! call invRstress

    stTmp = 0.0_dp
    call invR_stress(stTmp, nAtom_, coord, nNeighEwald_, &
        & iNeighbor,img2CentCell, gLatPoint_, &
        & alpha_, volume_, deltaQAtom_)

    st(:,:) = st(:,:) - 0.5_dp * stTmp(:,:)
    
    ! if (tExtChrg_) then
    ! ????
    ! end if
    
  end subroutine addStressDCSCC

  
  !!* Constructs the shift vectors for the SCC contributions
  !!* @param orb Contains information about the atomic orbitals in the system
  !!* @param specie       List of the specie for each atom.
  !!* @param iNeighbor    List of surrounding neighbours for each atom.
  !!* @parma img2CentCell Image of each atom in the central cell.
  subroutine buildShifts_(orb, specie, iNeighbor, img2CentCell)
    type(TOrbitals), intent(in) :: orb
    integer,  intent(in)    :: specie(:)
    integer,  intent(in)    :: iNeighbor(0:,:)
    integer,  intent(in)    :: img2CentCell(:)

    integer :: iAt1, iSp1, iSh1, iU1, iNeigh, iAt2f, iSp2, iSh2, iU2

    ASSERT(tInitialised_)
    
    !! 1/R contribution [shiftPerAtom(A) = \sum_B 1/R_AB * (Q_B - Q0_B)]
    shiftPerAtom_(:) = 0.0_dp
    call hemv(shiftPerAtom_, invRMat_, deltaQAtom_,'L')
    
    !! Contribution of the short range part of gamma to the shift
    !! sgamma'_{A,l} = sum_B sum_{u\in B} S(U(A,l),u)*q_u
    shiftPerL_(:,:) = 0.0_dp
    do iAt1 = 1, nAtom_
      iSp1 = specie(iAt1)
      do iSh1 = 1, orb%nShell(iSp1)
        iU1 = iHubbU_(iSh1, iSp1)
        do iNeigh = 0, maxval(nNeighShort_(:,:,:,iAt1))
          iAt2f = img2CentCell(iNeighbor(iNeigh, iAt1))
          iSp2 = specie(iAt2f)
          do iSh2 = 1, orb%nShell(iSp2)
            iU2 = iHubbU_(iSh2, iSp2)
            shiftPerL_(iSh1, iAt1) = shiftPerL_(iSh1, iAt1) &
                &- shortGamma_(iU2, iU1, iNeigh, iAt1) &
                &* deltaQPerLShell_(iSh2, iAt2f)
            if (iAt2f /= iAt1) then
              shiftPerL_(iSh2, iAt2f) = shiftPerL_(iSh2, iAt2f) &
                  &- shortGamma_(iU2, iU1, iNeigh, iAt1) &
                  &* deltaQPerLShell_(iSh1, iAt1)
            end if
          end do
        end do
      end do
    end do

    if (tChrgConstr_) then
      call buildShift(chrgConstr_, deltaQAtom_)
    end if
    if (tThirdOrder_) then
      call buildShift(thirdOrder_, deltaQAtom_)
    end if
    
  end subroutine buildShifts_



  !!* Returns the shift per atom coming from the SCC part (with a spin
  !!* index)
  !!* @param shift Contains the shift on exit.
  subroutine getShiftPerAtom(shift)
    real(dp), intent(out) :: shift(:,:)

    ASSERT(size(shift,dim=1) == size(shiftPerAtom_,dim=1))

    shift(:,1) = shiftPerAtom_
    if (tExtChrg_) then
      call addShiftPerAtom_ExtChrg(shift(:,1))
    end if
    if (tChrgConstr_) then
      call addShiftPerAtom(chrgConstr_, shift(:,1))
    end if
    if (tThirdOrder_) then
      call addShiftPerAtom(thirdOrder_, shift(:,1))
    end if
    
  end subroutine getShiftPerAtom

  
  
  !!* Returns the shift per L contribution of the SCC. (with a spin index)
  !!* @param shift Contains the shift on exit.
  subroutine getShiftPerL(shift)
    real(dp), intent(out) :: shift(:,:,:)

    ASSERT(size(shift,dim=1) == size(shiftPerL_,dim=1))
    ASSERT(size(shift,dim=2) == size(shiftPerL_,dim=2))
    ASSERT(size(shift,dim=3) > 0)

    shift(:,:,1) = shiftPerL_(:,:)

  end subroutine getShiftPerL

  

  !!* Returns the equivalency relations between the orbitals of the atoms.
  !!* @param orb Contains information about the atomic orbitals in the system
  !!* @param specie Type of each atom (nAtom).
  !!* @param equiv The vector describing the equivalence on return.
  subroutine SCC_getOrbitalEquiv(orb, specie, equiv)
    type(TOrbitals), intent(in) :: orb
    integer, intent(in) :: specie(:)
    integer, intent(out) :: equiv(:,:,:)

    integer :: iAt, iOrb, iS, iSp
    integer :: nSpin, shift

    nSpin = size(equiv, dim=3)

    ASSERT(tInitialised_)
    ASSERT(size(specie) == nAtom_)
    ASSERT(size(equiv, dim=1) == orb%mOrb)
    ASSERT(all(shape(equiv) == (/ orb%mOrb, nAtom_, nSpin /)))

    equiv(:,:,:) = 0
    shift = 0
    do iAt = 1, nAtom_
      iSp = specie(iAt)
      do iOrb = 1, orb%nOrbSpecie(iSp)
        equiv(iOrb, iAt, 1) = iHubbU_(orb%iShellOrb(iOrb, iSp), iSp) + shift
      end do
      shift = shift + maxval(iHubbU_(:, iSp))
    end do
    do iS = 2, nSpin
      equiv(:,:,iS) = equiv(:,:,1)
    end do

  end subroutine SCC_getOrbitalEquiv


end module SCC
