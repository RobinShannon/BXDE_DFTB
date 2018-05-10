!!* Contains routines related to finite electron temperature, including 
!!* Fermi, Gaussian and Methfessel-Paxton broadening functions.
!!* @todo Add other methods, including possibly Pederson and Jackson method
!!* PRB 43, 7312 (1991). Also fix exact occupation for electron numers, using
!!* interpolation instead of bisection.
module eTemp
# include "allocate.h"
# include "assert.h"    
  use accuracy, only : dp, elecTol, elecTolMax, mExpArg
  use external_erfc
  use message
  use Hermite
  use constants, only : pi
  use factorial, only : fact
  implicit none

  private

  integer, parameter :: Fermi = 0             !* Definition of a type of
                                              !* broadening function -
                                              !* Fermi-Dirac in this case
  integer, parameter :: Gaussian = 1          !* Definition of a type of
                                              !* broadening function -
                                              !* Gaussian in this case
  integer, parameter :: Methfessel = 1        !* Definition of a type of
                                              !* broadening function -
                                              !* Methfessel-Paxton,
                                              !* for higher orders use
                                              !* Methfessel + n as a value

  public :: Efilling, Fermi, Gaussian, Methfessel
  

  real(dp), parameter :: epsilon2 = 2.0_dp * epsilon(1.0_dp)


  
contains

  !!* Driver to calculate electron filling, the band-structure energy at T and
  !!* extrapolated to T=0K, and the entropy of the electron energy for the
  !!* Mermin free energy.
  !!* @param Ebs Band structure energy at T
  !!* @param Ef Fermi energy for given distribution
  !!* @param TS Entropy
  !!* @param E0 Band structure energy extrapolated to T=0K 
  !!* @param filling Electron occupancies
  !!* @param eigenvals The eigenvalues of the levels, 1st index is energy
  !!*   2nd index is k-point and 3nd index is spin
  !!* @param nelectrons Number of electrons
  !!* @param kT Thermal energy in atomic units
  !!* @param kWeight k-point weightings
  !!* @param distrib Choice of distribution functions, currently
  !!*   Fermi, Gaussian and Methfessle-Paxton supported. The flags is defined
  !!*   symbolically, so (Methfessel + 2) gives the 2nd order M-P scheme
  !!* @note use slices of eigenvalues and other input arrays if you want a
  !!*   different Fermi energy and filling for different k-points and/or spins.
  !!* @note If no electrons are present, the Fermi energy is set to zero per
  !!*   default.
  subroutine Efilling(Ebs,Ef,TS,E0,filling,eigenvals,nelectrons,kT,kWeight, &
      & distrib)
    real(dp), intent(out) :: Ebs
    real(dp), intent(out) :: Ef
    real(dp), intent(out) :: TS
    real(dp), intent(out) :: E0
    real(dp), intent(out) :: filling(:,:,:)
    real(dp), intent(in) :: eigenvals(:,:,:)
    real(dp), intent(in) :: nelectrons
    real(dp), intent(in) :: kT
    real(dp), intent(in) :: kWeight(:)
    integer, intent(in) :: distrib
    
    real(dp) :: upperEf, lowerEf
    real(dp) :: nElec
    real(dp) :: nElecMax, nElecMin, maxEig, minEig
    real(dp) :: EfOld
    
    ASSERT(size(eigenvals,dim=1) > 0)
    ASSERT(all(shape(filling) == shape(eigenvals)))
    ASSERT(size(eigenvals,dim=3) >= 1)
    ASSERT(size(eigenvals,dim=3) <= 2)
    ASSERT(nelectrons >= 0.0_dp)
    ! Not a tight enough bound ? :
    ASSERT(ceiling(nelectrons) <= size(eigenvals,dim=1))
    ASSERT(kT > 0.0_dp)
    ASSERT(size(kWeight) > 0)
    ASSERT(all(kWeight >= 0.0_dp))

    ASSERT(distrib >= Fermi)

    !! If no electrons there, we are ready
    if (nelectrons < epsilon(1.0_dp)) then
      filling(:,:,:) = 0.0_dp
      Ebs = 0.0_dp
      Ef = 0.0_dp
      TS = 0.0_dp
      E0 = 0.0_dp
      return
    end if
    
    !! find maximum and minimum possible value of Fermi Energy
    minEig = minval(eigenvals(1,:,:))
    maxEig = maxval(eigenvals(size(eigenvals, dim=1),:,:))
    !! Fermi level hopefully between highest and lowest eigenvalue
    upperEf = maxEig + 0.01_dp
    lowerEf = minEig - 0.01_dp

    !! but just to be on the safe side if the temperature is
    !! BIG compared to the bandwidth, or if the system has
    !! a fully filled band structure:
    nElecMax = electronCount(upperEf, eigenvals, kT, distrib, kWeight)
    nElecMin = electronCount(lowerEf, eigenvals, kT, distrib, kWeight)
    do while (nElecMin > nelectrons)
      lowerEf = 2.0_dp * (lowerEf-upperEf) + lowerEf
      nElecMin = electronCount(lowerEf, eigenvals, kT, distrib, kWeight)
    end do
    do while (nElecMax < nelectrons)
      upperEf = 2.0_dp * (upperEf-lowerEf) + lowerEf
      nElecMax = electronCount(upperEf, eigenvals, kT, distrib, kWeight)
    end do
    
    Ef = 0.5_dp * (upperEf + lowerEf)
    nElec = electronCount(Ef, eigenvals, kT, distrib, kWeight)

    !! Bisection as long as nr. electrons is not accurate enough or
    !! next change in the Fermi level would go below precision
    do while (abs(nElectrons - nElec) > elecTol &
        &.and. abs(upperEf - lowerEf) >= max(abs(Ef)*epsilon2, epsilon2))
      if ((nElecMax >= nElecMin) .eqv. (nelectrons >= nElec)) then
        lowerEf = Ef
        nElecMin = nElec
      else
        upperEf = Ef
        nElecMax = nElec
      end if
      Ef = 0.5_dp*(upperEf + lowerEf)
      nElec = electronCount(Ef, eigenvals, kT, distrib, kWeight)
    end do

    !! If number of electrons deviates from theoretical value too much: stop
    if (abs(nElectrons - nElec) > elecTolMax) then
      call error("Fermi level search did not converge.")
    end if

    ! Polish resulting root with Newton-Raphson type steps
    if (abs(nElectrons - nElec) > elecTol) then
      if (distrib == Fermi) then ! only derivs for Fermi so far
        if (abs(derivElectronCount(Ef,eigenvals,kT,distrib,kWeight))>= &
            & epsilon(1.0_dp))&
            & then
          EfOld = Ef
          Ef = Ef - (electronCount(Ef,eigenvals, kT,distrib,kWeight) - &
              & nElectrons)/derivElectronCount(Ef,eigenvals,kT,distrib,kWeight)
          do while( abs( &
              & electronCount(EfOld,eigenvals, kT,distrib,kWeight)-nElectrons)&
              & > &
              & abs(electronCount(Ef,eigenvals, kT,distrib,kWeight)-nElectrons))
            if (abs(derivElectronCount(Ef,eigenvals,kT,distrib,kWeight)) >= &
                & epsilon(1.0_dp)) then
              EfOld = Ef
              Ef = Ef - (electronCount(Ef,eigenvals, kT,distrib,kWeight) - &
                  & nElectrons) / &
                  & derivElectronCount(Ef,eigenvals,kT,distrib,kWeight)
            else
              exit
            end if
          end do
          Ef = EfOld
        end if
      end if
    end if
    
    nElec = electronCount(Ef, eigenvals, kT, distrib, kWeight)
    
    call electronFill(Ebs,filling,TS,E0,Ef,eigenvals,kT,distrib,kWeight)
    ! re-scale to give exact number of electrons, this is a temporay hack
    if (nElec > epsilon(1.0_dp)) then
      filling(:,:,:) = filling(:,:,:) * nelectrons/nElec
    end if
    
  end subroutine Efilling



!!* Calculates the number of electrons for a given Fermi energy and 
!!* distribution function
!!* @param Ef Fermi energy for given distribution
!!* @param eigenvals The eigenvalues of the levels, 1st index is energy
!!* 2nd index is k-point and 3nd index is spin
!!* @param kT Thermal energy in atomic units
!!* @param distrib Choice of distribution functions, currently
!!* Fermi, Gaussian and Methfessle-Paxton supported. The flags is defined
!!* sumbolically, so (Methfessel + 2) gives the 2nd order M-P scheme
!!* @param kWeight k-point weightings
  function electronCount(Ef,eigenvals,kT,distrib,kWeight)
    real(dp) :: electronCount
    real(dp), intent(in) :: Ef 
    real(dp), intent(in) :: eigenvals(:,:,:)
    real(dp), intent(in) :: kT
    integer, intent(in) :: distrib
    real(dp), intent(in) :: kWeight(:)

    integer :: MPorder
    real(dp) :: w
    real(dp), allocatable :: A(:)
    real(dp), allocatable :: hermites(:)
    integer i, j , k, l, ispin
    real(dp) :: occ, x    

    w = 1.0_dp/kT
    electronCount=0.0_dp
    if (distrib /= Fermi) then
      MPorder = distrib - 1
      ALLOCATE_(A,(0:MPorder))
      ALLOCATE_(hermites,(0:2*MPorder))
      call Aweights(A,MPorder)
      do ispin = 1, size(eigenvals,dim=3)
        do i = 1, size(kWeight)
          do j = 1, size(eigenvals,dim=1)
            if (eigenvals(j,i,ispin)>(Ef-3.0_dp*w)) then
              exit
            else
              electronCount=electronCount+kWeight(i)
            end if
          end do
          do k = j, size(eigenvals,dim=1)
            if (eigenvals(k,i,ispin)>(Ef+3.0_dp*w)) then
              exit
            end if            
            x = ( eigenvals(k,i,ispin) - Ef ) / kT
            call hX(hermites,MPorder*2,x)
#ifdef EXTERNALERFC
            occ = 0.5_dp*extErfc(x)
#else
            occ = 0.5_dp*erfc(x)
#endif
            do l=1, MPorder
              occ = occ + A(l) * hermites(2*l-1) * exp(-x**2)
            end do
            electronCount = electronCount + occ * kWeight(i)
          end do
        end do
      end do
      DEALLOCATE_(A)
      DEALLOCATE_(hermites)
    else
      do ispin = 1, size(eigenvals,dim=3)
        do i = 1, size(kWeight)
          do j = 1, size(eigenvals,dim=1)
            x = ( eigenvals(j,i,ispin) - Ef ) / kT
! Where the compiler does not handle inf gracefully, trap the exponential
! function for small input values
#ifdef EXPTRAP 
            if (x <= mExpArg) then            
              electronCount = electronCount + kWeight(i)/(1.0_dp + exp(x))
            endif
#else
            electronCount = electronCount + kWeight(i)/(1.0_dp + exp(x))
#endif
          end do
        end do
      end do
    end if
  end function electronCount
  
  !!* Calculates the derivative of the number of electrons for a given Fermi
  !!* energy and distribution function
  !!* @param Ef Fermi energy for given distribution
  !!* @param eigenvals The eigenvalues of the levels, 1st index is energy
  !!* 2nd index is k-point and 3nd index is spin
  !!* @param kT Thermal energy in atomic units
  !!* @param distrib Choice of distribution functions, currently
  !!* Fermi supported.
  !!* @param kWeight k-point weightings
  !!* @todo support MP
  function derivElectronCount(Ef,eigenvals,kT,distrib,kWeight)
    real(dp) :: derivElectronCount
    real(dp), intent(in) :: Ef
    real(dp), intent(in) :: eigenvals(:,:,:)
    real(dp), intent(in) :: kT
    integer, intent(in) :: distrib
    real(dp), intent(in) :: kWeight(:)
    
    real(dp) :: w
    integer i, j, ispin
    real(dp) :: x
    
    w = 1.0_dp/kT
    derivElectronCount=0.0_dp
    if (distrib /= Fermi) then
      call error("Fermi distribution only supported")
    else
      do ispin = 1, size(eigenvals,dim=3)
        do i = 1, size(kWeight)
          do j = 1, size(eigenvals,dim=1)
            x = ( eigenvals(j,i,ispin) - Ef ) * w
            if (x<10.0_dp) then
              ! Where the compiler does not handle inf gracefully, 
              ! trap the exponential function for small input values
#ifdef EXPTRAP
              if (x <= mExpArg) then
                derivElectronCount = derivElectronCount + &
                    & (w*kWeight(i)) * (exp(x)/((1.0_dp + exp(x))**2))
              endif
#else
              derivElectronCount = derivElectronCount + &
                  & (w*kWeight(i)) * (exp(x)/((1.0_dp + exp(x))**2))
#endif
            end if
          end do
        end do
      end do
    end if
  end function derivElectronCount
  
!!* Calculate filling and TS for the given eigenspectrum and distribution
!!* function and Fermi energy.
!!* @ref G. Kresse and J. Furthm&uuml;ller, Phys. Rev. B vol 54, pp 11169
!!* (1996).
!!* @ref M. Methfessel and A. T. Paxton,, Phys. Rev. B vol 40, pp 3616 (1989)
!!* @ref F. Wagner, Th.\ Laloyaux and M. Scheffler, Phys. Rev. B, vol 57 pp
!!* 2102 (1998)
!!* @param Eband Band structure energy at T
!!* @param filling Electron occupancies
!!* @param TS Entropy
!!* @param E0 Band structure energy extrapolated to T=0K 
!!* @param Ef Fermi energy for given distribution
!!* @param eigenvals The eigenvalues of the levels, 1st index is energy
!!* 2nd index is k-point and 3nd index is spin
!!* @param kT Thermal energy in atomic units
!!* @param distrib Choice of distribution functions, currently
!!* Fermi, Gaussian and Methfessle-Paxton supported. The flags is defined
!!* sumbolically, so (Methfessel + 2) gives the 2nd order M-P scheme
!!* @param kWeight k-point weightings
  subroutine electronFill(Eband,filling,TS,E0,Ef,eigenvals,kT,distrib,kWeights)
    real(dp), intent(out) :: Eband
    real(dp), intent(out) :: filling(:,:,:)
    real(dp), intent(out) :: TS
    real(dp), intent(out) :: E0
    real(dp), intent(in) :: Ef 
    real(dp), intent(in) :: eigenvals(:,:,:)
    real(dp), intent(in) :: kT
    integer, intent(in) :: distrib
    real(dp), intent(in) :: kWeights(:)

    integer :: MPorder
    integer :: kpts
    real(dp) :: w
    real(dp), allocatable :: A(:)
    real(dp), allocatable :: hermites(:)
    integer i, j , k, l, ispin
    real(dp) :: occ, x

    kpts = size(kWeights)
    
    Eband = 0.0_dp
    TS = 0.0_dp
    filling(:,:,:)=0.0_dp
    w = 1.0_dp/kT
    E0 = 0.0_dp
    
    ! The Gaussian and Methfessel-Paxton broadening functions first
    if (distrib /= Fermi) then
      MPorder = distrib - 1
      ALLOCATE_(A,(0:MPorder))
      ALLOCATE_(hermites,(0:2*MPorder))
      call Aweights(A,MPorder)
      do ispin = 1, size(eigenvals,dim=3)
        do i = 1, kpts
          do j = 1, size(eigenvals,dim=1)
            if (eigenvals(j,i,ispin)>(Ef-3.0_dp*w)) then
              exit
            else
              filling(j,i,ispin)=1.0_dp
              Eband = Eband + eigenvals(j,i,ispin)
            end if
          end do
          do k = j, size(eigenvals,dim=1)
            if (eigenvals(k,i,ispin)>(Ef+3.0_dp*w)) then
              exit
            end if
            x = ( eigenvals(k,i,ispin) - Ef ) / kT
            call hX(hermites,MPorder*2,x)
            ! Gauusian broadened occupancy
#ifdef EXTERNALERFC
            occ = 0.5_dp*extErfc(x)
#else
            occ = 0.5_dp*erfc(x)
#endif
            ! Gaussian broadening entropy
            TS = TS + kWeights(i)*0.5_dp*exp(-x**2)/sqrt(pi)
            ! Methfessel-Paxton occupation sum
            do l=1, MPorder
              occ = occ + A(l) * hermites(2*l-1) * exp(-x**2)
            end do
            filling(k,i,ispin) = occ
            ! Sum up the band-structure energy, including k-point weight where
            ! needed
            Eband = Eband + kWeights(i)*filling(k,i,ispin) * eigenvals(k,i&
                &,ispin)
            ! Methfessel-Paxton broadening entropy
            do l=1, MPorder
              TS = TS + kWeights(i)*0.5_dp*A(l)*hermites(2*l)*exp(-x**2)
            end do
          end do
        end do
      end do
      DEALLOCATE_(A)
      DEALLOCATE_(hermites)
      TS = TS * kT
      E0 = (real(MPorder+1,dp)*(Eband - TS) + Eband)/real(MPorder+2,dp)
    else
      do ispin = 1, size(eigenvals,dim=3)
        do i = 1, kpts
          do j = 1, size(eigenvals,dim=1)            
            x = ( eigenvals(j,i,ispin) - Ef ) / kT
! Where the compiler does not handle inf gracefully, trap the exponential
! function for small values
#ifdef EXPTRAP
            if (x > mExpArg) then
              filling(j,i,ispin) = 0.0_dp
            else
              filling(j,i,ispin) = 1.0_dp/(1.0_dp + exp(x))
            endif
#else
            filling(j,i,ispin) = 1.0_dp/(1.0_dp + exp(x))
#endif
            if (filling(j,i,ispin) <= elecTol) then
              exit
            end if
            if (filling(j,i,ispin) < (1.0_dp - elecTol) ) then
              ! Fermi-Dirac entropy :
              TS = TS - kWeights(i)*(filling(j,i,ispin)* &
                  & log(filling(j,i,ispin))+(1.0_dp-filling(j,i,ispin))&
                  & *log(1.0_dp-filling(j,i,ispin)))
            end if
            Eband = Eband + kWeights(i)*(filling(j,i,ispin) * eigenvals(j,i&
                &,ispin))
          end do
        end do
      end do
      TS = TS * kT
      E0 = Eband - 0.5_dp * TS
    end if
  end subroutine electronFill

!!* Calculate the weighting factors for the Methfessel-Paxton smearing scheme
!!* @param A returned weighting values for the scheme, given by 
!!* $A_n = \frac{(-1)^n}{n!4^n\sqrt{\pi]}$
!!* @param n the required order to calculate $A_n$ up to
!!* @ref M. Methfessel and A. T. Paxton, Phys. Rev. B Vol 40, pp 3616 (1989)
  subroutine Aweights(A,n)
    real(dp), intent(out) :: A(0:)
    integer, intent(in) :: n
    real(dp) :: nbang(0:n)
    integer i
    ASSERT(n>=0)
    ASSERT(size(A)>=n)
    A(:) = 0.0_dp
    call fact(nbang,n)
    do i = 0, n
       A(i) = real((-1)**i,dp)/(nbang(i)*real(4**i,dp)*sqrt(pi))
    end do
  end subroutine Aweights
  
end module eTemp
