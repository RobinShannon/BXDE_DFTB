!!* Common routines for the dispersion modules.
!!* Periodic summation from the following references:
!!* @ref N. Karasawa et al., J. Phys. Chem. 93, 7320-7327 (1989)
!!* @ref Zhou-Min Chen et al., J. Comp. Chem. 18, 1365 (1997)
module DispCommon
#include "assert.h"  
  use Accuracy
  use Constants, only : pi
  use Message
  use Sorting
  use SimpleAlgebra, only : cross3
#ifdef EXTERNALERFC
  use external_erfc
#endif
  implicit none
  private

  public :: addDispEGr_per_species, addDispEGr_per_atom
  public :: getOptimalEta, getMaxRDispersion, getMaxGDispersion


contains

  !!* Adds the energy per atom and the gradients for periodic 1/r^6 summation.
  !!* @param nAtom Nr. of atoms (without periodic images)
  !!* @param coords Coordinates of the atoms (including images)
  !!* @param nNeighbors Nr. of neighbors for each atom
  !!* @param iNeighbor Neighborlist.
  !!* @param neighDist2 Square distances of the neighbours.
  !!* @param img2CentCell Mapping of periodic images onto the central cell
  !!* @param c6 Van der Waals coefficients (nAtom, nAtom)
  !!* @param eta Controling partitioning between real and reciprocal space sum
  !!* @param vol Volume of the unit cell
  !!* @param gLatVecs Set of reciprocal space vectors to include in the sum
  !!* @param energies Updated energy vector at return
  !!* @param gradients Updated gradient vector at return
  !!* @param st Updated stress tensor
  !!* @note Interaction parameter C6 is specified atomwise.
  !!* @desc Fast converging Ewald like summation on 1/r^6 type interactions.
  !!*   For details see the references in the module description.
  subroutine addDispEGr_per_atom(nAtom, coords, nNeighbors, &
      &iNeighbor, neighDist2, img2CentCell, c6, eta, vol, gLatVecs, &
      &energies, gradients,st)
    integer, intent(in) :: nAtom
    real(dp), intent(in) :: coords(:,:)
    integer, intent(in) :: nNeighbors(:)
    integer, intent(in) :: iNeighbor(0:,:)
    real(dp), intent(in) :: neighDist2(0:,:)
    integer, intent(in) :: img2CentCell(:)
    real(dp), intent(in) :: c6(:,:)
    real(dp), intent(in) :: eta
    real(dp), intent(in) :: vol
    real(dp), intent(in) :: gLatVecs(:,:)
    real(dp), intent(inout) :: energies(:)
    real(dp), intent(inout) :: gradients(:,:)
    real(dp), intent(inout) :: st(:,:)

    integer :: iAt1, iNeigh, iAt2, iAt2f, iG, ii
    real(dp) :: rSum, rSum3(3), gSum, gSum3(3), gg(3), ggAbs, vec(3)
    real(dp) :: aam2, bb, bbm2, rTmp, rTmp2, rTmp3, rc, r3c, gc, g3c, ddp
    real(dp) :: etam3, rTmp33, gsum33(3,3)

    ASSERT(size(energies) == nAtom)
    ASSERT(all(shape(gradients) == (/ 3, nAtom /)))
    ASSERT(all(shape(st) == (/ 3, 3 /)))

    etam3 = eta**(-3)
    rc = 0.5_dp * etam3 * etam3
    r3c = 2.0_dp * rc / (eta * eta)
    gc = pi**1.5_dp / (12.0_dp * vol)
    g3c = pi**1.5_dp / (6.0_dp * vol)
    do iAt1 = 1, nAtom
      do iNeigh = 1, nNeighbors(iAt1)
        iAt2 = iNeighbor(iNeigh, iAt1)
        vec = coords(:,iAt1) - coords(:,iAt2)
        iAt2f = img2CentCell(iAt2)
        aam2 = (sqrt(neighDist2(iNeigh, iAt1))/eta)**(-2)
        rSum =  rc * c6(iAt2f, iAt1) &
            &* ((aam2 + 1.0_dp)*aam2 + 0.5_dp)*aam2 * exp(-1.0_dp / aam2)
        rSum3(:) = r3c * c6(iAt2f, iAt1) * exp(-1.0_dp / aam2) &
            &* (((6.0_dp*aam2 + 6.0_dp)*aam2 + 3.0_dp)*aam2 + 1.0_dp)*aam2 &
            & * vec(:)
        energies(iAt1) = energies(iAt1) - rSum
        gradients(:,iAt1) = gradients(:,iAt1) + rSum3(:)
        if (iAt1 /= iAt2f) then
          energies(iAt2f) = energies(iAt2f) - rSum
          gradients(:,iAt2f) = gradients(:,iAt2f) - rSum3(:)
          do ii = 1, 3
            st(:,ii) = st(:,ii) - rSum3(:)*vec(ii) / vol
          end do
        else
          do ii = 1, 3
            st(:,ii) = st(:,ii) - 0.5_dp*rSum3(:)*vec(ii) / vol
          end do
        end if
      end do
      
      gSum = 0.0_dp
      gSum3(:) = 0.0_dp
      gSum33(:,:) = 0.0_dp
      do iG = 1, size(gLatVecs, dim=2)
        gg = gLatVecs(:, iG)
        ggAbs = sqrt(sum(gg**2))
        bb = 0.5_dp * ggAbs * eta
        bbm2 = bb**(-2)
        rTmp = 0.0_dp
        rTmp2 = 0.0_dp
        do iAt2 = 1, nAtom
          ddp = dot_product(gg, coords(:,iAt1) - coords(:,iAt2))
          rTmp = rTmp + cos(ddp) * c6(iAt1, iAt2)
          rTmp2 = rTmp2 + sin(ddp) * c6(iAt1, iAt2)
        end do
#ifdef EXTERNALERFC
        rTmp3 = sqrt(pi) * extErfc(bb) + (0.5_dp * bbm2 - 1.0_dp) / bb &
            &* exp(-1.0_dp / bbm2)
        rTmp33 = sqrt(pi) * extErfc(bb) - exp(-bb*bb)/ bb 
#else        
        rTmp3 = sqrt(pi) * erfc(bb) + (0.5_dp * bbm2 - 1.0_dp) / bb &
            &* exp(-1.0_dp / bbm2)
        rTmp33 = sqrt(pi) * erfc(bb) - exp(-bb*bb)/ bb 
#endif        
        gSum = gSum + rTmp * ggAbs**3 * rTmp3
        gSum3(:) = gSum3(:) + gg(:) * rTmp2 * ggAbs**3 * rTmp3
        do ii = 1, 3
          gSum33(:,ii) = gSum33(:,ii) + rTmp * 3.0_dp*ggAbs*rTmp33*gg(:)*gg(ii)
        end do
      end do
      gSum = gSum * gc
      gSum3(:) = gSum3(:) * g3c
      gSum33 = gSum33 * gc
      energies(iAt1) = energies(iAt1) - gSum &
          &+ (c6(iAt1,iAt1)/12.0_dp * etam3 &
          &- pi**1.5 * sum(c6(:,iAt1))/(6.0_dp * vol)) * etam3 
      do ii = 1, 3
        st(ii,ii) = st(ii,ii) - gSum/vol &
            &-(pi**1.5 * sum(c6(:,iAt1))/(6.0_dp * vol*vol)) * etam3
      end do
      st = st - gSum33 / vol
      gradients(:,iAt1) =  gradients(:,iAt1) +  gSum3
    end do

  end subroutine addDispEGr_per_atom


  
  !!* Adds the energy per atom and the gradients for periodic 1/r^6 summation.
  !!* @param nAtom Nr. of atoms (without periodic images)
  !!* @param coords Coordinates of the atoms (including images)
  !!* @param nNeighbors Nr. of neighbors for each atom
  !!* @param iNeighbor Neighborlist.
  !!* @param neighDist2 Square distances of the neighbours.
  !!* @param img2CentCell Mapping of periodic images onto the central cell
  !!* @param c6 Van der Waals coefficients (nSpecie, nSpecie)
  !!* @param eta Controling partitioning between real and reciprocal space sum
  !!* @param vol Volume of the unit cell
  !!* @param gLatVecs Set of reciprocal space vectors to include in the sum
  !!* @param energies Updated energy vector at return
  !!* @param gradients Updated gradient vector at return
  !!* @param st Updated stress tensor
  !!* @note Interaction coefficients (c6) is specified speciewise.
  !!* @desc Fast converging Ewald like summation on 1/r^6 type interactions.
  !!*   The 1/r^12 term is summed in direct space.
  subroutine addDispEGr_per_species(nAtom, coords, species0, &
      &nNeighbors, iNeighbor, neighDist2, img2CentCell, c6, eta, vol, &
      &gLatVecs, energies, gradients, st)
    integer, intent(in) :: nAtom
    real(dp), intent(in) :: coords(:,:)
    integer, intent(in) :: species0(:)
    integer, intent(in) :: nNeighbors(:)
    integer, intent(in) :: iNeighbor(0:,:)
    real(dp), intent(in) :: neighDist2(0:,:)
    integer, intent(in) :: img2CentCell(:)
    real(dp), intent(in) :: c6(:,:)
    real(dp), intent(in) :: eta
    real(dp), intent(in) :: vol
    real(dp), intent(in) :: gLatVecs(:,:)
    real(dp), intent(inout) :: energies(:)
    real(dp), intent(inout) :: gradients(:,:)
    real(dp), intent(inout) :: st(:,:)
    
    integer :: iAt1, iNeigh, iAt2, iAt2f, iG, iSp1, iSp2, ii
    real(dp) :: rSum, rSum3(3), gSum, gSum3(3),gsum33(3,3),gg(3), ggAbs, vec(3)
    real(dp) :: aam2, bb, bbm2, rTmp, rTmp2, rTmp3, rc, r3c, gc, g3c, ddp, etam3
    real(dp) :: rTmp33
    
    ASSERT(size(energies) == nAtom)
    ASSERT(all(shape(gradients) == (/ 3, nAtom /)))
    ASSERT(all(shape(st) == (/ 3, 3 /)))
    
    etam3 = eta**(-3)
    rc = 0.5_dp * etam3 * etam3
    r3c = 2.0_dp * rc / (eta * eta)
    gc = pi**1.5_dp / (12.0_dp * vol)
    g3c = pi**1.5_dp / (6.0_dp * vol)
    do iAt1 = 1, nAtom
      iSp1 = species0(iAt1)
      do iNeigh = 1, nNeighbors(iAt1)
        iAt2 = iNeighbor(iNeigh, iAt1)
        vec = coords(:,iAt1) - coords(:,iAt2)
        iAt2f = img2CentCell(iAt2)
        iSp2 = species0(iAt2f)
        aam2 = (sqrt(neighDist2(iNeigh, iAt1))/eta)**(-2)
        rSum =  rc * c6(iSp2, iSp1) &
            &* ((aam2 + 1.0_dp)*aam2 + 0.5_dp)*aam2 * exp(-1.0_dp / aam2)
        rSum3(:) = r3c * c6(iSp2, iSp1) * exp(-1.0_dp / aam2) &
            &* (((6.0_dp*aam2 + 6.0_dp)*aam2 + 3.0_dp)*aam2 + 1.0_dp)*aam2 &
            & * vec(:)
        energies(iAt1) = energies(iAt1) - rSum
        gradients(:,iAt1) = gradients(:,iAt1) + rSum3(:)
        if (iAt1 /= iAt2f) then
          energies(iAt2f) = energies(iAt2f) - rSum
          gradients(:,iAt2f) = gradients(:,iAt2f) - rSum3(:)
          do ii = 1, 3
            st(:,ii) = st(:,ii) - rSum3 * vec(ii) / vol
          end do
        else
          do ii = 1, 3
            st(:,ii) = st(:,ii) - 0.5_dp * rSum3 * vec(ii) / vol
          end do
        end if
      end do
      
      gSum = 0.0_dp
      gSum3(:) = 0.0_dp
      gSum33(:,:) = 0.0_dp
      do iG = 1, size(gLatVecs, dim=2)
        gg = gLatVecs(:, iG)
        ggAbs = sqrt(sum(gg**2))
        bb = 0.5_dp * ggAbs * eta
        bbm2 = bb**(-2)
        rTmp = 0.0_dp
        rTmp2 = 0.0_dp
        do iAt2 = 1, nAtom
          iSp2 = species0(iAt2)
          ddp = dot_product(gg, coords(:,iAt1) - coords(:,iAt2))
          rTmp = rTmp + cos(ddp) * c6(iSp1, iSp2)
          rTmp2 = rTmp2 + sin(ddp) * c6(iSp1, iSp2)
        end do
#ifdef EXTERNALERFC
        rTmp3 = sqrt(pi) * extErfc(bb) + (0.5_dp * bbm2 - 1.0_dp) / bb &
            &* exp(-1.0_dp / bbm2)
        rTmp33 = sqrt(pi) * extErfc(bb) - exp(-bb*bb)/ bb 
#else        
        rTmp3 = sqrt(pi) * erfc(bb) + (0.5_dp * bbm2 - 1.0_dp) / bb &
            &* exp(-1.0_dp / bbm2)
        rTmp33 = sqrt(pi) * erfc(bb) - exp(-bb*bb)/ bb 
#endif        
        gSum = gSum + rTmp * ggAbs**3 * rTmp3
        gSum3(:) = gSum3(:) + gg(:) * rTmp2 * ggAbs**3 * rTmp3
        do ii = 1, 3
          gSum33(:,ii) = gSum33(:,ii) + rTmp * 3.0_dp*ggAbs*rTmp33*gg(:)*gg(ii)
        end do
      end do
      gSum = gSum * gc
      gSum3(:) = gSum3(:) * g3c
      gSum33 = gSum33 * gc
      energies(iAt1) = energies(iAt1) - gSum &
          &+ (c6(iSp1,iSp1)/12.0_dp * etam3 &
          &- pi**1.5 * sum(c6(species0(1:nAtom),iSp1))/(6.0_dp * vol)) * etam3
      do ii = 1, 3
        st(ii,ii) = st(ii,ii) - gSum/vol &
            &-(pi**1.5 * sum(c6(species0(1:nAtom),iSp1))/(6.0_dp * vol*vol)) * etam3
      end do
      st = st - gSum33 / vol
      gradients(:,iAt1) =  gradients(:,iAt1) +  gSum3
    end do

  end subroutine addDispEGr_per_species



  !!* Delivers the optimal paramater eta for the Ewald-like summation
  !!* @param latVecs Lattice vectors
  !!* @param vol Volume
  !!* @return Optimal parameter.
  function getOptimalEta(latVecs, vol) result(eta)
    real(dp), intent(in) :: latVecs(:,:)
    real(dp), intent(in) :: vol
    real(dp) :: eta

    real(dp) :: vecLens(3), tmp(3)
    integer :: indx(3), i1, i2

    vecLens(:) = sqrt(sum(latVecs**2, dim=1))
    call index_heap_sort(indx, vecLens)
    !! *** Workaround for old compaq compilers ***
    i1 = indx(1)
    i2 = indx(2)
    call cross3(tmp, latVecs(:,i1), latVecs(:,i2))
    eta = sqrt(vecLens(indx(1)) * vol / (pi * sqrt(sum(tmp**2))))

  end function getOptimalEta

    

  !!* Returns the longest real space vector needed to achive a given accuracy
  !!* in the Ewald summation for the dispersion.
  !!* @param eta Parameter of the ewald summation.
  !!* @param c6sum Sum of the Van der Waals coefficients (for every possible
  !!*   interaction in the system).
  !!* @param vol Volume of the unit cell
  !!* @param minValue  Tolerance value.
  !!* @return Cutoff radius.
  function getMaxRDispersion(eta, c6sum, vol, minValue) result(xx)
    real(dp), intent(in) :: eta
    real(dp), intent(in) :: c6sum
    real(dp), intent(in) :: vol, minValue
    real(dp) :: xx

    real(dp), parameter :: rInit = 1.0e-8_dp
    real(dp) :: xLeft, xRight, yLeft, yRight, yy
    integer  :: iError, iIter
    character(lc) :: strError

    iError = 0
    xx = rInit
    yy = getDispRealError_(xx, c6sum, vol, eta)
    do while (yy > minValue .and. xx <= huge(1.0_dp))
      xx = 2.0_dp * xx
      yy = getDispRealError_(xx, c6sum, vol, eta)
    end do
    if (xx > huge(1.0_dp)) then
      iError = 1
    elseif (xx == rInit) then
      iError = 2
    end if

    if (iError == 0) then
      xLeft = 0.5_dp * xx
      yLeft = getDispRealError_(xLeft, c6sum, vol, eta)
      xRight = xx
      yRight = yy

      iIter = 1
      do while (yLeft - yRight > minValue .and. iIter <= nSearchIter)
        xx = 0.5_dp * (xLeft + xRight)
        yy = getDispRealError_(xx, c6sum, vol, eta)
        if (yy >= minValue) then
          xLeft = xx
          yLeft = yy
        else
          xRight = xx
          yRight = yy
        end if
        iIter = iIter + 1
      end do
      if (iIter > nSearchIter) then
        iError = 3
      end if
    end if

    if (iError /= 0) then
99020 format ('Failure in getMaxRDispersion_.', ' Error nr: ',I3)
      write(strError, 99020) iError
      call error(strError)
    end if

  end function getMaxRDispersion



  !!* Returns the longest reciprocal space vector needed to achive a given 
  !!* accuracy in the Ewald summation for the dispersion.
  !!* @param eta Parameter of the ewald summation.
  !!* @param c6sum Sum of the absolute values of the c6 coeffs for every atom
  !!*   pair.
  !!* @param minValue  Tolerance value.
  !!* @return Cutoff radius.
  function getMaxGDispersion(eta, c6sum, minValue) result(xx)
    real(dp), intent(in) :: eta
    real(dp), intent(in) :: c6sum
    real(dp), intent(in) :: minValue
    real(dp) :: xx

    real(dp), parameter :: gInit = 1.0e-8_dp
    real(dp) :: xLeft, xRight, yLeft, yRight, yy
    integer :: iError, iIter
    character(lc) :: strError

    iError = 0
    xx = gInit
    yy = getDispReciprocalError_(xx, c6sum, eta)
    do while (yy > minValue .and. xx <= huge(1.0_dp))
      xx = 2.0_dp * xx
      yy = getDispReciprocalError_(xx, c6sum, eta)
    end do
    if (xx > huge(1.0_dp)) then
      iError = 1
    elseif (xx == gInit) then
      iError = 2
    end if

    if (iError == 0) then
      xLeft = 0.5_dp * xx
      yLeft = getDispReciprocalError_(xLeft, c6sum, eta)
      xRight = xx
      yRight = yy

      iIter = 1
      do while (yLeft - yRight > minValue .and. iIter <= nSearchIter)
        xx = 0.5_dp * (xLeft + xRight)
        yy = getDispReciprocalError_(xx, c6sum, eta)
        if (yy >= minValue) then
          xLeft = xx
          yLeft = yy
        else
          xRight = xx
          yRight = yy
        end if
        iIter = iIter + 1
      end do
      if (iIter > nSearchIter) then
        iError = 3
      end if
    end if

    if (iError /= 0) then
99010 format ('Failure in getMaxGDispersion_.', ' Error nr: ',I3)
      write(strError, 99010) iError
      call error(strError)
    end if

  end function getMaxGDispersion

  
  !!* Returns the error of the real space summation for a certain cutoff
  !!* @param rr Cutoff radius
  !!* @param c6sum VdW coefficient sum.
  !!* @param vol Volume of the supercell
  !!* @param eta Summation parameter
  !!* @return Estimated error of the summation.
  function getDispRealError_(rr, c6sum, vol, eta) result(err)
    real(dp), intent(in) :: rr, c6sum, vol, eta
    real(dp) :: err

#ifdef EXTERNALERFC    
    err = (pi**1.5_dp * c6sum / vol) * eta * (1.0_dp/rr**4 &
        &+ 1.0_dp/(rr**2 * eta**2) + 1.0_dp / (2.0_dp * eta**4)) &
        &* extErfc(rr/eta)
#else
    err = (pi**1.5_dp * c6sum / vol) * eta * (1.0_dp/rr**4 &
        &+ 1.0_dp/(rr**2 * eta**2) + 1.0_dp / (2.0_dp * eta**4)) &
        &* erfc(rr/eta)
#endif    
    
  end function getDispRealError_

  
  !!* Returns the error of the reciprocal space summation for a certain cutoff
  !!* @param gg Cutoff radius
  !!* @param c6sum VdW coefficient sum.
  !!* @param eta Summation parameter
  !!* @return Estimated error of the summation.
  function getDispReciprocalError_(gg, c6sum, eta) result(err)
    real(dp), intent(in) :: gg, c6sum, eta
    real(dp) :: err

#ifdef EXTERNALERFC
    err = c6sum/(6.0_dp * sqrt(pi)) * (1.0_dp / eta**6) &
        &*(gg * eta * exp(-1.0_dp * (0.5_dp*gg*eta)**2) &
        &+ sqrt(pi) * extErfc(0.5_dp * gg * eta))
#else
    err = c6sum/(6.0_dp * sqrt(pi)) * (1.0_dp / eta**6) &
        &*(gg * eta * exp(-1.0_dp * (0.5_dp*gg*eta)**2) &
        &+ sqrt(pi) * erfc(0.5_dp * gg * eta))
#endif
    
  end function getDispReciprocalError_

end module DispCommon
