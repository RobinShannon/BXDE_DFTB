!!* Various I/O routines for the main program.
module MainIO
#include "allocate.h"
#include "assert.h"  
  use Accuracy
  use Constants
  use Periodic
  use CommonTypes
  use Fifo
  use Sparse2Dense
  use BlasRoutines
  use CharManip, only : i2c
  use LinkedList
  use FileId
  implicit none
  private

  public :: writeRealEigvecs, writeCplxEigvecs
  public :: writeProjRealEigvecs, writeProjCplxEigvecs
  public :: getH
  
  character(*), parameter :: eigvecOut = "eigenvec.out"
  character(*), parameter :: eigvecBin = "eigenvec.bin"
  character(*), parameter :: regionOut = "region_"

  !!* Routines to get eigenvectors out of storage/memory
  interface getH
    module procedure getHreal
    module procedure getHcmplx
  end interface
  
contains

  !!* Write the real eigenvectors into text and binary output files.
  !!* @param fdEigvec  Fileid (file not yet opened) to use.
  !!* @param runId  Id of the current program run.
  !!* @param nAtom  Nr. of atoms in the system.
  !!* @param nSpin  Nr. of spin channels.
  !!* @param neighlist  Neighbor list.
  !!* @param nNeighbor  Nr. of neighbors for SK-interaction.
  !!* @param iAtomStart  Positions of atoms int the dense matrices.
  !!* @param iPair  Positions of interactions in the sparse matrices.
  !!* @param img2CentCell  Mapping of atoms into the central cell.
  !!* @param orb  Orbital information.
  !!* @param specie  Species.
  !!* @param specieName  Name of the species.
  !!* @param over  Sparse overlap matrix.
  !!* @param HSqrReal  Square Hamiltonian (or work array)
  !!* @param SSqrReal  Work array.
  !!* @param storeEigvecs  If present, Hamiltonian(s) are fetched from this
  !!* storage into HSqrReal, instead of using whatever is already there.
  subroutine writeRealEigvecs(fdEigvec, runId, nAtom, nSpin, neighlist, &
      &nNeighbor, iAtomStart, iPair, img2CentCell, orb, specie, specieName, &
      &over, HSqrReal, SSqrReal, storeEigvecs)
    integer, intent(in) :: fdEigvec, runId, nAtom, nSpin
    type(TNeighborList), intent(in) :: neighlist
    integer, intent(in) :: nNeighbor(:)
    integer, intent(in) :: iAtomStart(:), iPair(:,:), img2CentCell(:)
    type(TOrbitals), intent(in) :: orb
    integer, intent(in) :: specie(:)
    character(mc), intent(in) :: specieName(:)
    real(dp), intent(in) :: over(:)
    real(dp), intent(inout) :: HSqrReal(:,:,:), SSqrReal(:,:)
    type(OFifoRealR2), intent(inout), optional :: storeEigvecs(:)
    
    character(lc) :: tmpStr
    integer :: iSpin, iSpin2, iAtom, iSp1, iSh1, iOrb, ang
    integer :: ii, jj
    real(dp), allocatable :: rVecTemp(:)

    close(fdEigvec) ! just to be on the safe side
    ! Write eigenvalues in binary form
    open(fdEigvec, file=eigvecBin, action="write", status="replace", &
        &position="rewind", form="unformatted")
    write (fdEigVec) runId
    do iSpin = 1, nSpin
      call getH(iSpin, HSqrReal, iSpin2, storeEigvecs)
      do ii = 1, size(HSqrReal, dim=2)                
        write (fdEigvec) HSqrReal(:,ii,iSpin2)
      end do
    end do
    close(fdEigvec)

    ! Write eigenvalues (together with Mulliken populations) in text form
    open(fdEigvec, file=eigvecOut, action="write", status="replace", &
        &position="rewind")
    write (fdEigvec,"(A/)") "Coefficients and Mulliken populations&
        & of the atomic orbitals"

    ALLOCATE_(rVecTemp, (size(HSqrReal, dim=1)))
    call unpackHS(SSqrReal, over, neighlist%iNeighbor, nNeighbor, &
        &iAtomStart, iPair, img2CentCell)
    do iSpin = 1, nSpin
      call getH(iSpin, HSqrReal, iSpin2, storeEigvecs)
      do ii = 1, orb%nOrb
        call hemv(rVecTemp, SSqrReal, HSqrReal(:,ii,iSpin2))
        write(fdEigvec, "('Eigenvector:',I4,4X,'(',A,')'/)") ii, &
            &trim(spinName(iSpin2))
        jj = 0
        do iAtom = 1, nAtom
          iSp1 = specie(iAtom)
          do iSh1 = 1, orb%nShell(iSp1)
            ang = orb%angShell(iSh1, iSp1)
            if (iSh1 == 0) then
              write(tmpStr, "(I5,1X,A2,2X,A1)") iAtom, specieName(iSp1), &
                  &orbitalNames(ang+1)
            else
              write(tmpStr, "(10X,A1)") orbitalNames(ang+1)
            end if
            do iOrb = 1, 2 * ang + 1
              jj = jj + 1
              write(fdEigvec,"(A,I1,T15,F12.6,3X,F12.6)") trim(tmpStr),&
                  &iOrb, HSqrReal(jj, ii, iSpin2), &
                  & HSqrReal(jj, ii, iSpin2) * rVecTemp(jj)
            end do
          end do
          write (fdEigvec,*)
        end do
      end do
    end do
    close(fdEigvec)
    DEALLOCATE_(rVecTemp)
    
  end subroutine writeRealEigvecs


  
  !!* Write the complex eigenvectors into text and binary output files.
  !!* @param fdEigvec  Fileid (file not yet opened) to use.
  !!* @param runId  Id of the current program run.
  !!* @param nAtom  Nr. of atoms in the system.
  !!* @param nSpin  Nr. of spin channels.
  !!* @param neighlist  Neighbor list.
  !!* @param nNeighbor  Nr. of neighbors for SK-interaction.
  !!* @param cellVec  Cell vectors of shifted cells.
  !!* @param iCellVec  Cell vector index of every atom.
  !!* @param iAtomStart  Positions of atoms int the dense matrices.
  !!* @param iPair  Positions of interactions in the sparse matrices.
  !!* @param img2CentCell  Mapping of atoms into the central cell.
  !!* @param orb  Orbital information.
  !!* @param specie  Species.
  !!* @param specieName  Name of the species.
  !!* @param over  Sparse overlap matrix.
  !!* @param kpoint  KPoints.
  !!* @param HSqrCplx  Square Hamiltonian (or work array)
  !!* @param SSqrCplx  Work array.
  !!* @param storeEigvecs  If present, Hamiltonian(s) are fetched from this
  !!* storage into HSqrCplx, instead of using whatever is already there.
  subroutine writeCplxEigvecs(fdEigvec, runId, nAtom, nSpin, neighlist, &
      &nNeighbor, cellVec, iCellVec, iAtomStart, iPair, img2CentCell, orb, &
      &specie, specieName, over, kpoint, HSqrCplx, SSqrCplx, storeEigvecs)
    integer, intent(in) :: fdEigvec, runId, nAtom, nSpin
    type(TNeighborList), intent(in) :: neighlist
    integer, intent(in) :: nNeighbor(:)
    real(dp), intent(in) :: cellVec(:,:)
    integer, intent(in) :: iCellVec(:)
    integer, intent(in) :: iAtomStart(:), iPair(:,:), img2CentCell(:)
    type(TOrbitals), intent(in) :: orb
    integer, intent(in) :: specie(:)
    character(mc), intent(in) :: specieName(:)
    real(dp), intent(in) :: over(:), kpoint(:,:)    
    complex(dp), intent(inout) :: HSqrCplx(:,:,:,:), SSqrCplx(:,:)
    type(OFifoCplxR2), intent(inout), optional :: storeEigvecs(:)
    
    character(lc) :: tmpStr
    integer :: iSpin, iSpin2, iAtom, iSp1, iSh1, iOrb, ang, iK, iK2, nK
    integer :: ii, jj
    complex(dp), allocatable :: cVecTemp(:)

    nK = size(kPoint, dim=2)
    close(fdEigvec) ! just to be on the safe side
    ! Write eigenvalues in binary form
    open(fdEigvec, file=eigvecBin, action="write", status="replace",  &
        &position="rewind", form="unformatted")
    write (fdEigVec) runId
    do iSpin = 1, nSpin
      do iK = 1, nK
        call getH(iSpin, iK, HSqrCplx, iSpin2, iK2, storeEigvecs)
        do ii = 1, size(HSqrCplx, dim=2)
          write (fdEigvec) HSqrCplx(:,ii,iK2, iSpin2)
        end do
      end do
    end do
    close(fdEigvec)
    
    ! Write eigenvalues (together with Mullikan populations) in text form
    open(fdEigvec, file=eigvecOut, action="write", status="replace")
    write (fdEigvec,"(A/)") "Coefficients and Mulliken populations of the &
        &atomic orbitals"
    ALLOCATE_(cVecTemp,(size(HSqrCplx, dim=1)))
    do iSpin = 1, nSpin
      do iK = 1, nK
        call getH(iSpin, iK, HSqrCplx, iSpin2, iK2, storeEigvecs)
        call unpackHS(SSqrCplx, over, kPoint(:,iK), neighlist%iNeighbor, &
            &nNeighbor, iCellVec, cellVec, iAtomStart, iPair, img2CentCell)
        do ii = 1, orb%nOrb
          call hemv(cVecTemp, SSqrCplx, HSqrCplx(:,ii,iK2,iSpin2))
          write(fdEigvec, "(A,I4,4X,A,I4,4X,'(',A,')'/)") "K-point: ", ik, &
              &"Eigenvector: ", ii, trim(spinName(iSpin))
          jj = 0
          do iAtom = 1, nAtom
            iSp1 = specie(iAtom)
            do iSh1 = 1, orb%nShell(iSp1)
              ang = orb%angShell(iSh1,iSp1)
              if (iSh1 == 1) then
                write(tmpStr, "(I5,1X,A2,2X,A1)") iAtom, specieName(iSp1), &
                    &orbitalNames(ang+1)
              else
                write(tmpStr, "(10X,A1)") orbitalNames(ang+1)
              end if
              do iOrb = 1, 2*ang+1
                jj = jj + 1
                write(fdEigvec,&
                    &"(A,I1,T15,'(',F12.6,',',F12.6,')',3X,F12.6)") &
                    &trim(tmpStr), iOrb, &
                    &real(HSqrCplx(jj, ii, iK2, iSpin2)), &
                    &aimag(HSqrCplx(jj, ii, iK2, iSpin2)), &
                    & real( conjg(HSqrCplx(jj, ii, iK2, iSpin2)) &
                    & * cVecTemp(jj))
              end do
            end do
            write (fdEigvec,*)
          end do
        end do
      end do
    end do
    close(fdEigvec)
    DEALLOCATE_(cVecTemp)

  end subroutine writeCplxEigvecs

  !!* Write the projected eigenstates into text files.
  !!* @param filenames List with filenames for each region.
  !!* @param ei eigenvalues
  !!* @param nSpin  Nr. of spin channels.
  !!* @param neighlist  Neighbor list.
  !!* @param nNeighbor  Nr. of neighbors for SK-interaction.
  !!* @param iAtomStart  Positions of atoms int the dense matrices.
  !!* @param iPair  Positions of interactions in the sparse matrices.
  !!* @param img2CentCell  Mapping of atoms into the central cell.
  !!* @param orb  Orbital information.
  !!* @param over  Sparse overlap matrix.
  !!* @param HSqrReal  Square Hamiltonian (or work array)
  !!* @param SSqrReal  Work array.
  !!* @param iOrbRegion orbital number in each region
  !!* @param storeEigvecs  If present, Hamiltonian(s) are fetched from this
  !!* storage into HSqrReal, instead of using whatever is already there.
  subroutine writeProjRealEigvecs(filenames, ei, nSpin, neighlist, &
      &nNeighbor, iAtomStart, iPair, img2CentCell, orb, &
      &over, HSqrReal, SSqrReal, iOrbRegion, storeEigvecs)
    type(listCharLc), intent(inout) :: filenames
    real(dp), intent(in) :: ei(:,:,:)
    integer, intent(in) :: nSpin
    type(TNeighborList), intent(in) :: neighlist
    integer, intent(in) :: nNeighbor(:)
    integer, intent(in) :: iAtomStart(:), iPair(:,:), img2CentCell(:)
    type(TOrbitals), intent(in) :: orb
    real(dp), intent(in) :: over(:)
    real(dp), intent(inout) :: HSqrReal(:,:,:), SSqrReal(:,:)
    type(listIntR1), intent(inout) :: iOrbRegion
    type(OFifoRealR2), intent(inout), optional :: storeEigvecs(:)
    
    integer, allocatable :: fdProjEig(:), iOrbs(:)
    integer :: iSpin, iSpin2, iLev, ii, nReg, dummy
    integer :: valshape(1)
    real(dp) :: qState
    real(dp), allocatable :: rVecTemp(:)
    character(lc) :: tmpStr

    nReg = len(iOrbRegion)
    ASSERT(len(filenames) == nReg)

    ALLOCATE_(fdProjEig, (nReg))
    do ii = 1, nReg
      fdProjEig(ii) = getFileId()
      call get(filenames, tmpStr, ii)
      open(fdProjEig(ii), file=tmpStr, action="write", status="replace")
    end do
    
    ALLOCATE_(rVecTemp, (size(HSqrReal, dim=1)))
    call unpackHS(SSqrReal, over, neighlist%iNeighbor, nNeighbor, &
        &iAtomStart, iPair, img2CentCell)
    do iSpin = 1, nSpin
      do ii = 1, nReg
        if (nSpin <= 2) then
          write(fdProjEig(ii),*)' KPT',1,' SPIN ', iSpin
        else
          write(fdProjEig(ii),*)' KPT',1
        endif
      end do
      call getH(iSpin, HSqrReal, iSpin2, storeEigvecs)
      do iLev = 1, orb%nOrb
        call hemv(rVecTemp, SSqrReal, HSqrReal(:,iLev,iSpin2))
        rVecTemp = rVecTemp * HSqrReal(:,iLev,iSpin2)        
        do ii = 1, nReg
          call elemShape(iOrbRegion, valshape, ii)
          ALLOCATE_(iOrbs, (valshape(1)))
          call intoArray(iOrbRegion, iOrbs, dummy, ii)
          qState = sum(rVecTemp(iOrbs))
          write(fdProjEig(ii), "(f13.6,f10.6)") Hartree__eV * ei(iLev,1,iSpin),&
              & qState
          DEALLOCATE(iOrbs)
        end do
      end do
      if (iSpin < nSpin) then
        do ii = 1, nReg
          write(fdProjEig(ii),*)
        end do
      end if
    end do
        
    do ii = 1, nReg
      close(fdProjEig(ii))
    end do
    
    DEALLOCATE_(rVecTemp)
    DEALLOCATE_(fdProjEig)
    
  end subroutine writeProjRealEigvecs

  
  !!* Write the projected complex eigenstates into text files.
  !!* @param fdProjEig  Fileid (file not yet opened) to use.
  !!* @param filenames array of files to print out the regions
  !!* @param ei eigenvalues
  !!* @param nSpin  Nr. of spin channels.
  !!* @param neighlist  Neighbor list.
  !!* @param nNeighbor  Nr. of neighbors for SK-interaction.
  !!* @param cellVec  Cell vectors of shifted cells.
  !!* @param iCellVec  Cell vector index of every atom.
  !!* @param iAtomStart  Positions of atoms int the dense matrices.
  !!* @param iPair  Positions of interactions in the sparse matrices.
  !!* @param img2CentCell  Mapping of atoms into the central cell.
  !!* @param orb  Orbital information.
  !!* @param over  Sparse overlap matrix.
  !!* @param kpoint  KPoints.
  !!* @param kweight KPoints weights
  !!* @param HSqrCplx  Square Hamiltonian (or work array)
  !!* @param SSqrCplx  Work array.
  !!* @param nRegion number of orbitals for each region
  !!* @param iOrbRegion orbital number in each region
  !!* @param storeEigvecs  If present, Hamiltonian(s) are fetched from this
  !!* storage into HSqrReal, instead of using whatever is already there.
  subroutine writeProjCplxEigvecs(filenames, ei, nSpin, neighlist, &
      & nNeighbor, cellVec, iCellVec, iAtomStart, iPair, img2CentCell, orb, &
      & over, kpoint, kweight, HSqrCplx, SSqrCplx, iOrbRegion, storeEigvecs)
    type(ListCharLc), intent(inout) :: filenames
    real(dp), intent(in) :: ei(:,:,:)
    integer, intent(in) :: nSpin
    type(TNeighborList), intent(in) :: neighlist
    integer, intent(in) :: nNeighbor(:)
    real(dp), intent(in) :: cellVec(:,:)
    integer, intent(in) :: iCellVec(:)
    integer, intent(in) :: iAtomStart(:), iPair(:,:), img2CentCell(:)
    type(TOrbitals), intent(in) :: orb
    real(dp), intent(in) :: over(:), kpoint(:,:), kweight(:)
    complex(dp), intent(inout) :: HSqrCplx(:,:,:,:), SSqrCplx(:,:)
    type(listIntR1), intent(inout) :: iOrbRegion
    type(OFifoCplxR2), intent(inout), optional :: storeEigvecs(:)

    integer, allocatable :: fdProjEig(:), iOrbs(:)
    integer :: iSpin, iSpin2, iK, iK2, nK, iLev, ii, nReg, dummy
    integer :: valshape(1)
    real(dp) :: qState
    complex(dp), allocatable :: cVecTemp(:)
    character(lc) :: tmpStr

    nK = size(kPoint, dim=2)
    nReg = len(iOrbRegion)
    ASSERT(len(filenames) == nReg)
    ASSERT(size(kweight) == nK)

    ALLOCATE_(fdProjEig, (nReg))
    do ii = 1, nReg
      fdProjEig(ii) = getFileId()
      call get(filenames, tmpStr, ii)
      open(fdProjEig(ii), file=tmpStr, action="write", status="replace")
    end do
    
    ALLOCATE_(cVecTemp,(size(HSqrCplx, dim=1)))
    do iSpin = 1, nSpin
      do iK = 1, nK
        if (nSpin <= 2) then
          do ii = 1, nReg            
            write(fdProjEig(ii),*)'KPT ',iK,' SPIN ', iSpin, &
                &' KWEIGHT ', kweight(iK)
          end do
        else
          do ii = 1, nReg
            write(fdProjEig(ii),*)'KPT ',iK, ' KWEIGHT ', kweight(iK)
          end do
        end if
        call getH(iSpin, iK, HSqrCplx, iSpin2, iK2, storeEigvecs)
        call unpackHS(SSqrCplx, over, kPoint(:,iK), neighlist%iNeighbor, &
            & nNeighbor, iCellVec, cellVec, iAtomStart, iPair, img2CentCell)
        do iLev = 1, orb%nOrb
          call hemv(cVecTemp, SSqrCplx, HSqrCplx(:,iLev,iK2,iSpin2))
          cVecTemp = conjg(HSqrCplx(:,iLev,iK2,iSpin2)) * cVecTemp
          do ii = 1, nReg
            call elemShape(iOrbRegion, valshape, ii)
            ALLOCATE_(iOrbs, (valshape(1)))
            call intoArray(iOrbRegion, iOrbs, dummy, ii)
            qState = real(sum(cVecTemp(iOrbs)), dp)
            write(fdProjEig(ii), "(f13.6,f10.6)") &
                & Hartree__eV * ei(iLev,iK,iSpin), qState
            DEALLOCATE_(iOrbs)
          end do
        end do
        if (iK < nK .or. iSpin < nSpin) then
          do ii = 1, nReg
            write(fdProjEig(ii),*)
          end do
        end if
      end do
    end do
        
    do ii = 1, nReg
      close(fdProjEig(ii))
    end do
    
    DEALLOCATE_(cVecTemp)
    DEALLOCATE_(fdProjEig)

    
  end subroutine writeProjCplxEigvecs

  
  subroutine getHreal(iSpin, HSqrReal, iSpin2, storeEigvecs)
    integer, intent(in) :: iSpin
    real(dp), intent(inout) :: HSqrReal(:,:,:)
    integer, intent(out) :: iSpin2
    type(OFifoRealR2), intent(inout), optional :: storeEigvecs(:)
    
    if (present(storeEigvecs)) then
      iSpin2 = 1
      call get(storeEigvecs(iSpin), HSqrReal(:,:,iSpin2))
    else
      iSpin2 = iSpin
    end if
    
  end subroutine getHreal

  
  subroutine getHcmplx(iSpin, iK, HSqrCplx, iSpin2, iK2, storeEigvecs)
    integer, intent(in) :: iSpin, iK
    complex(dp), intent(inout) :: HSqrCplx(:,:,:,:)
    integer, intent(out) :: iSpin2, iK2
    type(OFifoCplxR2), intent(inout), optional :: storeEigvecs(:)
        
    if (present(storeEigvecs)) then
      iSpin2 = 1
      iK2 = 1
      call get(storeEigvecs(iSpin), HSqrCplx(:,:,iK2,iSpin2))
    else
      iSpin2 = iSpin
      iK2 = iK
    end if
    
  end subroutine getHcmplx
  
end module MainIO
