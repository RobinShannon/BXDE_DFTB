!!* Contains code to calculate the H0 Hamiltonian and overlap then pack it into
!!* the sparse storage format. Also code to differentiate H0 and S numerically.
!!* Finally, also the code to calculate the non-SCC total energy contribution.
!!* @author B. Hourahine, B. Aradi
!!* @desc Fills up the 1-dimensional sparse overlap and Hamiltonian (note the
!!* Hamiltonian is a 2-D array, with the second index reserved for spin and
!!* future extensions). Also suplies numerical derivatives of diatomic blocks
module nonscc
#include "assert.h"  
  use accuracy, only : dp, deltaXDiff
  use SK
  use SlakoCont
  use CommonTypes
  implicit none
  private

  public :: buildH0S, H0Sprime

  
contains

  !!* Driver for making the non-SCC Hamiltonian and the overlap matrices in the
  !!* primitive sparse format from the square form returned by rotateH0
  !!* @param ham  Returns the non-self-consistent Hamiltonian.
  !!* @param over  Returns the overlap matrix
  !!* @param skHamCont  Container for the SlaKo Hamiltonian integrals.
  !!* @param skOverCont  Container for the SlaKo overlap integrals.
  !!* @param selfegy  On-site energies for each species
  !!* @param coords list of all coordinates, including possible periodic images
  !!* of atoms
  !!* @param nNeighbors number of surrounding neighbors for each atom
  !!* @param iNeighbors list of surrounding neighbors for each atom
  !!* @param species chemical species of each atom
  !!* @param iPair Shift vector, where the interaction between two atoms
  !!* starts in the sparse format.
  !!* @param orb Information about the orbitals in the system.
  subroutine buildH0S(ham, over, skHamCont, skOverCont, selfegy, coords, &
      &nNeighbors, iNeighbors, species, iPair, orb)
    real(dp), intent(out) :: ham(:), over(:)
    type(OSlakoCont), intent(in) :: skHamCont, skOverCont
    real(dp), intent(in) :: selfegy(:,:)
    real(dp), intent(in) :: coords(:,:)
    integer, intent(in) :: nNeighbors(:)
    integer, intent(in) :: iNeighbors(0:,:)
    integer, intent(in) :: species(:)
    integer, intent(in) :: iPair(0:,:)
    type(TOrbitals), intent(in) :: orb

    real(dp) :: Htmp(orb%mOrb,orb%mOrb), Stmp(orb%mOrb,orb%mOrb)
    real(dp) :: vect(3), dist
    real(dp) :: interSKHam(getMIntegrals(skHamCont))
    real(dp) :: interSKOver(getMIntegrals(skOverCont))

    integer :: nAtom, nOrb1, nOrb2
    integer :: iAt1, iAt2, iSp1, iSp2, iNeigh1, iOrb1, ind

    nAtom = size(nNeighbors)

    ham(:) = 0.0_dp
    over(:) = 0.0_dp

    !! Put the on-site energies into the Hamiltonian,
    !! and <lm|l'm'> = delta_l,l' * delta_m',m' for the overlap	!'
    do iAt1 = 1, nAtom
      iSp1 = species(iAt1)
      ind = iPair(0,iAt1) + 1
      ! loop over the l range for the ith atom
      do iOrb1 = 1, orb%nOrbAtom(iAt1)
        over(ind) = 1.0_dp
        ham(ind) = selfegy(orb%iShellOrb(iOrb1, iSp1), iSp1)
        ind = ind + orb%nOrbAtom(iAt1) + 1
      end do
    end do

    !! Do the diatomic blocks for each of the atoms with its nNeighbors
    do iAt1 = 1, nAtom
      iSp1 = species(iAt1)
      nOrb1 = orb%nOrbSpecie(iSp1)
      do iNeigh1 = 1, nNeighbors(iAt1)
        iAt2 = iNeighbors(iNeigh1, iAt1)
        iSp2 = species(iAt2)
        nOrb2 = orb%nOrbSpecie(iSp2)
        ind = iPair(iNeigh1,iAt1)
        vect(:) = coords(:,iAt2) - coords(:,iAt1)
        dist = sqrt(sum(vect**2))
        vect(:) = vect(:)/dist
        call getSKIntegrals(skHamCont, interSKHam, dist, iSp1, iSp2)
        call rotateH0(Htmp, interSKHam, vect(1), vect(2), vect(3), &
            &iSp1, iSp2, orb)
        call getSKIntegrals(skOverCont, interSKOver, dist, iSp1, iSp2)
        call rotateH0(Stmp, interSKOver, vect(1), vect(2), vect(3), &
            &iSp1, iSp2, orb)
        ham(ind+1:ind+nOrb2*nOrb1) = &
            &reshape(Htmp(1:nOrb2,1:nOrb1), (/ nOrb2 * nOrb1 /))
        over(ind+1:ind+nOrb2*nOrb1) = &
            &reshape(Stmp(1:nOrb2,1:nOrb1), (/ nOrb2 * nOrb1 /))
      end do
    end do
    
  end subroutine buildH0S


  
  !!* Contains code to calculate the numerical derivative of a diatomic block
  !!* of the H0 Hamiltonian and overlap
  !!* @param Hprime returns the derivative of the non-self-consistent
  !!*   Hamiltonian diatomic block, with respect to x,y,z (last index)
  !!* @param Sprime returns the derivative of the overlap matrix diatomic
  !!*   block, with respect to x,y,z (last index)
  !!* @param skHamCont  Container for SK Hamiltonian integrals.
  !!* @param skOverCont  Container for SK overlap integrals.
  !!* @param coords list of all coordinates, including possible periodic images
  !!*   of atoms
  !!* @param species chemical species of each atom
  !!* @param atomI the first atom in the diatomic block
  !!* @param atomJ the second atom in the diatomic block
  !!* @param orb  Orbital informations.
  subroutine H0Sprime(Hprime, Sprime, skHamCont, skOverCont, coords, species, &
      &atomI, atomJ, orb)
    real(dp), intent(out) :: Hprime(:,:,:), Sprime(:,:,:)
    type(OSlakoCont), intent(in) :: skHamCont, skOverCont
    real(dp), intent(in) :: coords(:,:)
    integer, intent(in) :: species(:)
    integer, intent(in) :: atomI, atomJ
    type(TOrbitals), intent(in) :: orb

    real(dp) :: interSKHam(getMIntegrals(skHamCont))    ! interpolated H integs.
    real(dp) :: interSKOver(getMIntegrals(skOverCont))  ! interpolated S integs.
    real(dp) :: vect(3), dist
    integer  :: ii, jj
    integer :: sp1, sp2
    real(dp) :: Htmp(size(Hprime,dim=1), size(Hprime,dim=2),2,3)
    real(dp) :: Stmp(size(Sprime,dim=1), size(Sprime,dim=2),2,3)


    ASSERT(size(Hprime, dim=3) == 3)
    ASSERT(size(Sprime, dim=3) == 3)

    Hprime(:,:,:) = 0.0_dp
    Sprime(:,:,:) = 0.0_dp

    sp1 = species(atomI)
    sp2 = species(atomJ)
    
    do jj = 1, 3
      do ii = 1, 2
        vect(:) = coords(:,atomJ) - coords(:,atomI)
        vect(jj) = vect(jj) - real(2*ii-3,dp) * deltaXDiff
        dist = sqrt(sum(vect(:)**2))
        vect(:) = vect(:) / dist
        call getSKIntegrals(skHamCont, interSKHam, dist, sp1, sp2)
        call rotateH0(Htmp(:,:,ii,jj), interSKHam, vect(1), vect(2), &
            &vect(3), sp1, sp2, orb)
        call getSKIntegrals(skOverCont, interSKOver, dist, sp1, sp2)
        call rotateH0(Stmp(:,:,ii,jj), interSKOver, vect(1), vect(2), &
            &vect(3), sp1, sp2, orb)
      end do
    end do

    do ii = 1, 3
      Hprime(:,:,ii) = (Htmp(:,:,2,ii) - Htmp(:,:,1,ii)) / (deltaXDiff)
      Sprime(:,:,ii) = (Stmp(:,:,2,ii) - Stmp(:,:,1,ii)) / (deltaXDiff)
    end do
    
  end subroutine H0Sprime


  
!!  !!* Calculates the contribution of the non-scc part to the total energy.
!!  !!* @param h0           Zero order (non-scc) Hamiltonian.
!!  !!* @param dens         Density matrix in packed form.
!!  !!* @param nNeighbors    Nr. of neighbors for SK-interaction for each atom.
!!  !!* @param iNeighbor    List of neighbors for each atom.
!!  !!* @param img2CentCell Image of the atoms in the central cell.
!!  !!* @param iPair        Indexing array for the accessing h0 and density
!!  !!* matrix.  
!!  !!* @param specie       Specie of each atom.
!!  !!* @param orb  Orbitals
!!  !!* @note This function is not used any more!!! The main code (mis)uses the
!!  !!*   mullikan routine instead.
!!  function getEnergyH0(h0, dens, nNeighbors, iNeighbor, img2CentCell, iPair, &
!!      &specie, orb) result(eH0)
!!    real(dp), intent(in)  :: h0(:)
!!    real(dp), intent(in)  :: dens(:)
!!    integer,  intent(in)  :: nNeighbors(:)
!!    integer,  intent(in)  :: iNeighbor(0:,:)
!!    integer,  intent(in)  :: img2CentCell(:)
!!    integer,  intent(in)  :: iPair(0:,:)
!!    integer,  intent(in)  :: specie(:)
!!    type(TOrbitals), intent(in) :: orb
!!    real(dp) :: eH0
!!
!!    integer  :: nAtom
!!    integer  :: iAt1, iAt2f, iNeigh, iOrig, nOrb1, nOrb2
!!    real(dp) :: rTmp
!!    
!!    nAtom = size(nNeighbors)
!!    eH0 = 0.0_dp
!!    do iAt1 = 1, nAtom
!!      nOrb1 = orb%nOrbAtom(iAt1)
!!      do iNeigh = 0, nNeighbors(iAt1)
!!        iAt2f = img2CentCell(iNeighbor(iNeigh, iAt1))
!!        iOrig = iPair(iNeigh, iAt1)
!!        nOrb2 = orb%nOrbAtom(iAt2f)
!!        rTmp = sum(h0(iOrig+1:iOrig+nOrb1*nOrb2) &
!!            &* dens(iOrig+1:iOrig+nOrb1*nOrb2))
!!        if (iAt2f == iAt1) then
!!          eH0 = eH0 + rTmp  ! on diagonal block
!!        else
!!          ! off diagonal, so a contribution for the other triangle
!!          eH0 = eH0 + 2.0_dp * rTmp 
!!        end if
!!      end do
!!    end do
!!
!!  end function getEnergyH0

end module nonscc
