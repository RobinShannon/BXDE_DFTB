!!* Contains routines to calculate contributions to typical DFTB Hamiltonian
!!* parts using various generalisations of H_mu,nu = 0.5*S_mu,nu*(V_mu + V_nu)
module shift
#include "assert.h"
  use accuracy
  use CommonTypes

  implicit none

  private
  public :: add_shift, total_shift

  ! add shifts to a given Hamiltonian
  interface add_shift
    module procedure add_shift_atom
    module procedure add_shift_lshell
    module procedure add_shift_block
  end interface

  ! add the old style H_mu,nu = 0.5*S_mu,nu*(V_mu + V_nu) $\mu\in a,\nu\in a$
  ! not sure if these are correct for periodic boundaries !
  interface add_old_shift
    module procedure add_oldshift_atom
    module procedure add_oldshift_lshell
    module procedure add_oldshift_block
  end interface

  !!* Totals together shifts to get composites
  interface total_shift
    module procedure addatom_shell
    module procedure addshell_block
  end interface
  
contains

  !!* Regular atomic shift (potential is only dependent on number of atom)
  !!* @param ham          The resulting Hamiltonian contribution.
  !!* @param over         The overlap matrix.
  !!* @param nNeighbour   Number of neighbours surrounding each atom.
  !!* @param iNeighbour   List of neighbours for each atom.
  !!* @param species      List of the species of each atom.
  !!* @param orb Contains Information about the atomic orbitals in the system
  !!* @param iPair        Indexing array for the Hamiltonian.
  !!* @parma img2CentCell Index mapping atoms onto the central cell atoms.
  !!* @param shift        Shift to add at atom sites
  subroutine add_shift_atom(ham,over,nNeighbour,iNeighbour,species,orb,iPair, &
      & nAtom,img2CentCell,shift)
    real(dp), intent(inout)     :: ham(:,:)
    real(dp), intent(in)        :: over(:)
    integer, intent(in)         :: nNeighbour(:)
    integer, intent(in)         :: iNeighbour(0:,:)
    integer, intent(in)         :: species(:)
    type(TOrbitals), intent(in) :: orb
    integer, intent(in)         :: iPair(0:,:)
    integer, intent(in)         :: nAtom
    integer, intent(in)         :: img2CentCell(:)
    real(dp), intent(in)        :: shift(:,:)

    integer :: iAt1, iAt2, iAt2f, iOrig, iSp1, iSp2, nOrb1, nOrb2
    integer :: iNeigh, iSpin, nSpin

    ASSERT(size(ham,dim=1)==size(over))
    ASSERT(size(ham,dim=2)==size(shift,dim=2))
    ASSERT(size(nNeighbour)==nAtom)
    ASSERT(size(iNeighbour,dim=2)==nAtom)
    ASSERT(size(species)>=maxval(iNeighbour))
    ASSERT(size(species)<=size(img2CentCell))
    ASSERT(size(iPair,dim=1)>=(maxval(nNeighbour)+1))
    ASSERT(size(iPair,dim=2)==nAtom)
    ASSERT(size(shift,dim=1)==nAtom)

    nSpin = size(shift,dim=2)
    ASSERT(nSpin == 1 .or. nSpin == 2 .or. nSpin == 4)

    do iSpin = 1, nSpin
      do iAt1 = 1, nAtom
        iSp1 = species(iAt1)
        nOrb1 = orb%nOrbSpecie(iSp1)
        do iNeigh = 0, nNeighbour(iAt1)
          iAt2 = iNeighbour(iNeigh, iAt1)
          iAt2f = img2CentCell(iAt2)
          iSp2 = species(iAt2f)
          nOrb2 = orb%nOrbSpecie(iSp2)
          iOrig = iPair(iNeigh, iAt1)
          ham(iOrig+1:iOrig+nOrb2*nOrb1,iSpin) = &
              & ham(iOrig+1:iOrig+nOrb2*nOrb1,iSpin) + &
              & over(iOrig+1:iOrig+nOrb2*nOrb1) * 0.5_dp * &
              & ( shift(iAt1,iSpin) + shift(iAt2f,iSpin) )
        end do
      end do
    end do
    
  end subroutine add_shift_atom

  !!* l-dependent shift (potential is dependent on number of atom and l-shell)
  !!* @param ham          The resulting Hamiltonian contribution.
  !!* @param over         The overlap matrix.
  !!* @param nNeighbour   Number of neighbours surrounding each atom.
  !!* @param iNeighbour   List of neighbours for each atom.
  !!* @param species      List of the species of each atom.
  !!* @param orb Contains Information about the atomic orbitals in the system
  !!* @param iPair        Indexing array for the Hamiltonian.
  !!* @parma img2CentCell Index mapping atoms onto the central cell atoms.
  !!* @param shift        Shift to add for each l-shell on all atom sites,
  !!* (0:lmax,1:nAtom)
  subroutine add_shift_lshell( ham,over,nNeighbour,iNeighbour,species,orb, &
      & iPair,nAtom,img2CentCell,shift )
    real(dp), intent(inout)     :: ham(:,:)
    real(dp), intent(in)        :: over(:)
    integer, intent(in)         :: nNeighbour(:)
    integer, intent(in)         :: iNeighbour(0:,:)
    integer, intent(in)         :: species(:)
    type(TOrbitals), intent(in) :: orb
    integer, intent(in)         :: iPair(0:,:)
    integer, intent(in)         :: nAtom
    integer, intent(in)         :: img2CentCell(:)
    real(dp), intent(in)        :: shift(:,:,:)

    integer  :: iAt1, iAt2f, iOrig, iSp1, iSp2, nOrb1, nOrb2
    integer  :: iSh1, iSh2, iNeigh, iSpin, nSpin
    real(dp) :: tmpH(orb%mOrb,orb%mOrb), rTmp

    ASSERT(size(ham,dim=1)==size(over))
    ASSERT(size(nNeighbour)==nAtom)
    ASSERT(size(iNeighbour,dim=2)==nAtom)
    ASSERT(size(species)>=maxval(iNeighbour))
    ASSERT(size(species)<=size(img2CentCell))
    ASSERT(size(iPair,dim=1)>=(maxval(nNeighbour)+1))
    ASSERT(size(iPair,dim=2)==nAtom)
    ASSERT(size(shift,dim=1)==orb%mShell)
    ASSERT(size(shift,dim=2)==nAtom)
    ASSERT(size(ham,dim=2)==size(shift,dim=3))

    nSpin = size(shift,dim=3)

    do iSpin = 1, nSpin
      do iAt1= 1, nAtom
        iSp1 = species(iAt1)
        nOrb1 = orb%nOrbSpecie(iSp1)
        do iNeigh = 0, nNeighbour(iAt1)
          iAt2f = img2CentCell(iNeighbour(iNeigh, iAt1))
          iSp2 = species(iAt2f)
          nOrb2 = orb%nOrbSpecie(iSp2)
          iOrig = iPair(iNeigh, iAt1)
          do iSh1 = 1, orb%nShell(iSp1)
            rTmp = shift(iSh1, iAt1, iSpin)
            do iSh2 = 1, orb%nShell(iSp2)
              tmpH(orb%posShell(iSh2,iSp2):orb%posShell(iSh2+1,iSp2)-1, &
                  & orb%posShell(iSh1,iSp1):orb%posShell(iSh1+1,iSp1)-1) = &
                  & rTmp + shift(iSh2, iAt2f,iSpin)
            end do
          end do
          ham(iOrig+1:iOrig+nOrb2*nOrb1,iSpin) = &
              & ham(iOrig+1:iOrig+nOrb2*nOrb1,ispin) + &
              & 0.5_dp * over(iOrig+1:iOrig+nOrb2*nOrb1) &
              & * reshape(tmpH(1:nOrb2, 1:nOrb1), (/nOrb2*nOrb1/))
        end do
      end do
    end do
    
  end subroutine add_shift_lshell

  !!* shift depending on occupation-matrix like potentials. To use this for
  !!* lm-dependent potentials, use a diagonal shift matrix
  !!* @param ham          The resulting Hamiltonian contribution.
  !!* @param over         The overlap matrix.
  !!* @param nNeighbour   Number of neighbours surrounding each atom.
  !!* @param iNeighbour   List of neighbours for each atom.
  !!* @param species      List of the species of each atom.
  !!* @param orb Contains Information about the atomic orbitals in the system
  !!* @param iPair        Indexing array for the Hamiltonian.
  !!* @parma img2CentCell Index mapping atoms onto the central cell atoms.
  !!* @param shift        Shift to add at atom sites, listed as
  !!* (0:nOrb,0:nOrb,1:nAtom)
  subroutine add_shift_block( ham,over,nNeighbour,iNeighbour,species,orb, &
      & iPair,nAtom,img2CentCell,shift )
    real(dp), intent(inout)     :: ham(:,:)
    real(dp), intent(in)        :: over(:)
    integer, intent(in)         :: nNeighbour(:)
    integer, intent(in)         :: iNeighbour(0:,:)
    integer, intent(in)         :: species(:)
    type(TOrbitals), intent(in) :: orb
    integer, intent(in)         :: iPair(0:,:)
    integer, intent(in)         :: nAtom
    integer, intent(in)         :: img2CentCell(:)
    real(dp), intent(in)        :: shift(:,:,:,:)

    integer  :: iAt1, iAt2, iAt2f, iOrig, iSp1, iSp2, nOrb1, nOrb2
    integer :: iNeigh, iSpin, nSpin
    real(dp) :: tmpH(orb%mOrb,orb%mOrb), tmpS(orb%mOrb,orb%mOrb)
    
    ASSERT(size(ham,dim=1)==size(over))
    ASSERT(size(nNeighbour)==nAtom)
    ASSERT(size(iNeighbour,dim=2)==nAtom)
    ASSERT(size(species)>=maxval(iNeighbour))
    ASSERT(size(species)<=size(img2CentCell))
    ASSERT(size(iPair,dim=1)>=(maxval(nNeighbour)+1))
    ASSERT(size(iPair,dim=2)==nAtom)
    ASSERT(size(shift,dim=1)==orb%mOrb)
    ASSERT(size(shift,dim=2)==orb%mOrb)
    ASSERT(size(shift,dim=3)==nAtom)    
    ASSERT(size(ham,dim=2)==size(shift,dim=4))

    nSpin = size(shift,dim=4)

    do iSpin = 1, nSpin
      do iAt1 = 1, nAtom
        iSp1 = species(iAt1)
        nOrb1 = orb%nOrbSpecie(iSp1)
        do iNeigh = 0, nNeighbour(iAt1)
          iAt2 = iNeighbour(iNeigh, iAt1)
          iAt2f = img2CentCell(iAt2)
          iSp2 = species(iAt2f)
          nOrb2 = orb%nOrbSpecie(iSp2)
          iOrig = iPair(iNeigh, iAt1)
          tmpS(1:nOrb2,1:nOrb1) = reshape( &
              & over(iOrig+1:iOrig+nOrb2*nOrb1),(/nOrb2,nOrb1/) )
          tmpH(1:nOrb2,1:nOrb1) = 0.5_dp * ( &
              & matmul(tmpS(1:nOrb2,1:nOrb1), &
              & shift(1:nOrb1,1:nOrb1,iAt1,iSpin)) + &
              & matmul(shift(1:nOrb2,1:nOrb2,iAt2f,iSpin), &
              & tmpS(1:nOrb2,1:nOrb1)) )
          ham(iOrig+1:iOrig+nOrb2*nOrb1,iSpin) = &
              & ham(iOrig+1:iOrig+nOrb2*nOrb1,iSpin) + &
              & reshape(tmpH(1:nOrb2, 1:nOrb1), (/nOrb2*nOrb1/))
        end do
      end do
    end do
    
  end subroutine add_shift_block

  !!* Atomic shift (potential is only dependent on atomic label) using the
  !!* old non-interatomic form that was applied in original spin formalism
  !!* and original LDA+U
  !!* @param ham          The resulting Hamiltonian contribution.
  !!* @param over         The overlap matrix.
  !!* @param nNeighbour   Number of neighbours surrounding each atom.
  !!* @param iNeighbour   List of neighbours for each atom.
  !!* @param species      List of the species of each atom.
  !!* @param orb Contains Information about the atomic orbitals in the system
  !!* @param iPair        Indexing array for the Hamiltonian.
  !!* @parma img2CentCell Index mapping atoms onto the central cell atoms.
  !!* @param shift        Shift to add at atom sites
  subroutine add_oldshift_atom(ham,over,nNeighbour,iNeighbour,species,orb, &
      & iPair,nAtom,img2CentCell,shift)
    real(dp), intent(inout)     :: ham(:)
    real(dp), intent(in)        :: over(:)
    integer, intent(in)         :: nNeighbour(:)
    integer, intent(in)         :: iNeighbour(0:,:)
    integer, intent(in)         :: species(:)
    type(TOrbitals), intent(in) :: orb
    integer, intent(in)         :: iPair(0:,:)
    integer, intent(in)         :: nAtom
    integer, intent(in)         :: img2CentCell(:)
    real(dp), intent(in)        :: shift(:)

    integer :: iAt1, iAt2, iAt2f, iOrig, iSp1, nOrb1
    integer :: iNeigh

    do iAt1 = 1, nAtom
      iSp1 = species(iAt1)
      nOrb1 = orb%nOrbSpecie(iSp1)
      do iNeigh = 0, nNeighbour(iAt1)
        iAt2 = iNeighbour(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)
        if (iAt1==iAt2f) then
          iOrig = iPair(iNeigh, iAt1)
          ham(iOrig+1:iOrig+nOrb1*nOrb1) = ham(iOrig+1:iOrig+nOrb1*nOrb1) + &
              & over(iOrig+1:iOrig+nOrb1*nOrb1) * 0.5_dp * &
              & ( shift(iAt1) + shift(iAt2f) )
        end if
      end do
    end do

  end subroutine add_oldshift_atom

  !!* l-dependent shift (potential is dependent on atomic label and l-shell)
  !!* using the old non-interatomic form that was applied in original spin
  !!* formalism and original LDA+U
  !!* @param ham          The resulting Hamiltonian contribution.
  !!* @param over         The overlap matrix.
  !!* @param nNeighbour   Number of neighbours surrounding each atom.
  !!* @param iNeighbour   List of neighbours for each atom.
  !!* @param species      List of the species of each atom.
  !!* @param orb Contains Information about the atomic orbitals in the system
  !!* @param iPair        Indexing array for the Hamiltonian.
  !!* @parma img2CentCell Index mapping atoms onto the central cell atoms.
  !!* @param shift        Shift to add at l-shells on atom sites
  subroutine add_oldshift_lshell( ham,over,nNeighbour,iNeighbour,species,orb, &
      & iPair,nAtom,img2CentCell,shift )
    real(dp), intent(inout)     :: ham(:)
    real(dp), intent(in)        :: over(:)
    integer, intent(in)         :: nNeighbour(:)
    integer, intent(in)         :: iNeighbour(0:,:)
    integer, intent(in)         :: species(:)
    type(TOrbitals), intent(in) :: orb
    integer, intent(in)         :: iPair(0:,:)
    integer, intent(in)         :: nAtom
    integer, intent(in)         :: img2CentCell(:)
    real(dp), intent(in)        :: shift(:,:)


    integer  :: iAt1, iAt2f, iOrig, iSp1, nOrb1
    integer  :: iSh1, iSh2, iNeigh
    real(dp) :: tmpH(orb%mOrb,orb%mOrb), rTmp

    do iAt1= 1, nAtom
      iSp1 = species(iAt1)
      nOrb1 = orb%nOrbSpecie(iSp1)
      do iNeigh = 0, nNeighbour(iAt1)
        iAt2f = img2CentCell(iNeighbour(iNeigh, iAt1))
        if (iAt1==iAt2f) then
          iOrig = iPair(iNeigh, iAt1)
          do iSh1 = 1, orb%nShell(iSp1)
            rTmp = shift(iSh1, iAt1)
            do iSh2 = 1, orb%nShell(iSp1)
              tmpH(orb%posShell(iSh2,iSp1):orb%posShell(iSh2+1,iSp1)-1, &
                  & orb%posShell(iSh1,iSp1):orb%posShell(iSh1+1,iSp1)-1) = &
                  & rTmp + shift(iSh2, iAt2f)
            end do
          end do
          ham(iOrig+1:iOrig+nOrb1*nOrb1) = ham(iOrig+1:iOrig+nOrb1*nOrb1) + &
              & 0.5_dp * over(iOrig+1:iOrig+nOrb1*nOrb1) &
              & * reshape(tmpH(1:nOrb1, 1:nOrb1), (/nOrb1*nOrb1/))
        end if
      end do
    end do

  end subroutine add_oldshift_lshell

  !!* shift depending on occupation-matrix like potentials. To use this for
  !!* lm-dependent potentials, use a diagonal shift matrix
  !!* using the old non-interatomic form that was applied in original spin
  !!* formalism and original LDA+U
  !!* @param ham          The resulting Hamiltonian contribution.
  !!* @param over         The overlap matrix.
  !!* @param nNeighbour   Number of neighbours surrounding each atom.
  !!* @param iNeighbour   List of neighbours for each atom.
  !!* @param species      List of the species of each atom.
  !!* @param orb Contains Information about the atomic orbitals in the system
  !!* @param iPair        Indexing array for the Hamiltonian.
  !!* @parma img2CentCell Index mapping atoms onto the central cell atoms.
  !!* @param shift        Shift to add at atom sites
  subroutine add_oldshift_block( ham,over,nNeighbour,iNeighbour,species,orb, &
      & iPair,nAtom,img2CentCell,shift )
    real(dp), intent(inout)     :: ham(:)
    real(dp), intent(in)        :: over(:)
    integer, intent(in)         :: nNeighbour(:)
    integer, intent(in)         :: iNeighbour(0:,:)
    integer, intent(in)         :: species(:)
    type(TOrbitals), intent(in) :: orb
    integer, intent(in)         :: iPair(0:,:)
    integer, intent(in)         :: nAtom
    integer, intent(in)         :: img2CentCell(:)
    real(dp), intent(in)        :: shift(:,:,:)

    integer  :: iAt1, iAt2, iAt2f, iOrig, iSp1, nOrb1, iNeigh
    real(dp) :: tmpH(orb%mOrb,orb%mOrb), tmpS(orb%mOrb,orb%mOrb)

    do iAt1 = 1, nAtom
      iSp1 = species(iAt1)
      nOrb1 = orb%nOrbSpecie(iSp1)
      do iNeigh = 0, nNeighbour(iAt1)
        iAt2 = iNeighbour(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)
        if (iAt1==iAt2f) then
          iOrig = iPair(iNeigh, iAt1)
          tmpS(1:nOrb1,1:nOrb1) = reshape( &
              & over(iOrig+1:iOrig+nOrb1*nOrb1),(/nOrb1,nOrb1/) )
          tmpH(1:nOrb1,1:nOrb1) = ( &
              & matmul(tmpS(1:nOrb1,1:nOrb1),shift(1:nOrb1,1:nOrb1,iAt1)) + &
              & matmul(shift(1:nOrb1,1:nOrb1,iAt1),tmpS(1:nOrb1,1:nOrb1)) )
          ham(iOrig+1:iOrig+nOrb1*nOrb1) = ham(iOrig+1:iOrig+nOrb1*nOrb1) + &
              &  0.5_dp * reshape(tmpH(1:nOrb1, 1:nOrb1), (/nOrb1*nOrb1/))
        end if
      end do
    end do

  end subroutine add_oldshift_block

  subroutine addatom_shell(out,atom,shell,orb,species)
    real(dp), intent(out) :: out(:,:,:)
    real(dp), intent(in)  :: atom(:,:)
    real(dp), intent(in)  :: shell(:,:,:)
    type(TOrbitals), intent(in) :: orb
    integer, intent(in)   :: species(:)

    integer iAtom, iSpin, nAtom, nSpin

    nAtom = size(atom,dim=1)
    nSpin = size(atom,dim=2)
    
    ASSERT(all(shape(out)==shape(shell)))
    ASSERT(size(shell,dim=1)==orb%mShell)
    ASSERT(size(shell,dim=2)==nAtom)
    ASSERT(size(shell,dim=3)==nSpin)
    ASSERT(size(species)>=nAtom)

    do iSpin = 1, nSpin
      do iAtom = 1, nAtom
        out(1:orb%nShell(species(iAtom)),iAtom,iSpin) = &
            & shell(1:orb%nShell(species(iAtom)),iAtom,iSpin) &
            & + atom(iAtom,iSpin)
      end do
    end do
    
  end subroutine addatom_shell

  subroutine addshell_block(out,shell,atblock,orb,species)
    real(dp), intent(out) :: out(:,:,:,:)
    real(dp), intent(in)  :: shell(:,:,:)
    real(dp), intent(in)  :: atblock(:,:,:,:)
    type(TOrbitals), intent(in) :: orb
    integer, intent(in)   :: species(:)

    integer iAt, iSpin, nAtom, nSpin, iSh, iSp, iOrb

    nAtom = size(atblock,dim=3)
    nSpin = size(atblock,dim=4)
    
    ASSERT(all(shape(out)==shape(atblock)))
    ASSERT(size(atblock,dim=1)==orb%mOrb)
    ASSERT(size(atblock,dim=2)==orb%mOrb)
    ASSERT(size(shell,dim=1)==orb%mShell)
    ASSERT(size(shell,dim=2)==nAtom)
    ASSERT(size(shell,dim=3)==nSpin)
    ASSERT(size(species)>=nAtom)

    out(:,:,:,:) = atblock(:,:,:,:)

    do iSpin = 1, nSpin
      do iAt = 1, nAtom
        iSp = species(iAt)
        do iSh = 1, orb%nShell(iSp)
          do iOrb = orb%posShell(iSh,iSp),orb%posShell(iSh+1,iSp)-1
            out(iOrb,iOrb,iAt,iSpin) = out(iOrb,iOrb,iAt,iSpin) + &
                & shell(iSh,iAt,iSpin)
          end do
        end do
      end do
    end do
    
  end subroutine addshell_block

  
end module shift
