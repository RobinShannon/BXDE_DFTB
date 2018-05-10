!!* Contains an DIIS mixer
!!* @description
!!*   The DIIS mixing is done by building a weighted combination over the
!!*   previous input charges to minimise the residue of the error.  Only 
!!*   a specified number of previous charge vectors are considered.
!!* @note In order to use the mixer you have to create and reset it.
module DIISMixer
# include "allocate.h"
# include "assert.h"  
  use accuracy
  use lapackroutines, only : gesv
  implicit none
  
  private
  
  !!* Contains the necessary data for an DIIS mixer
  type ODIISMixer
    private
    real(dp) :: initMixParam                !!* Initial mixing parameter
    integer :: mPrevVector                  !!* Max. nr. of stored prev. vectors
    integer :: nPrevVector                  !!* Nr. of stored previous vectors
    integer :: nElem                        !!* Nr. of elements in the vectors
    integer :: indx                         !!* Index for the storage
    real(dp), pointer :: prevQInput(:,:)    !!* Stored previous input charges
    real(dp), pointer :: prevQDiff(:,:)     !!* Stored prev. charge differences
    logical :: tFromStart                   !!* True if DIIS used from iteration
                                            !!* 2 as well as mixing
  end type ODIISMixer
  
  
  !!* Creates an DIISMixer instance
  interface create
    module procedure DIISMixer_create
  end interface
  
  !!* Destroys an DIISMixer instance
  interface destroy
    module procedure DIISMixer_destroy
  end interface
  
  !!* Resets the mixer
  interface reset
    module procedure DIISMixer_reset
  end interface
  
  !!* Does the mixing
  interface mix
    module procedure DIISMixer_mix
  end interface
  
  
  public :: ODIISMixer
  public :: create, destroy, reset, mix
  
contains
  
  !!* Creates an DIIS mixer instance.
  !!* @param self         Pointer to an initialized DIIS mixer on exit
  !!* @param nGeneration  Nr. of generations (including actual) to consider
  !!* @param initMixParam Mixing parameter for the first nGeneration cycles
  !!* @param tFromStart   True if using DIIS from iteration 2 as well as mixing
  subroutine DIISMixer_create(self, nGeneration, initMixParam,tFromStart)
    type(ODIISMixer), pointer :: self
    integer, intent(in) :: nGeneration
    real(dp), intent(in) :: initMixParam
    logical, intent(in), optional :: tFromStart
    
    ASSERT(nGeneration >= 2)
    
    INITALLOCATE_P(self)
    self%nElem = 0
    self%mPrevVector = nGeneration !- 1
    
    INITALLOCATE_PARR(self%prevQInput, (self%nElem, self%mPrevVector))
    INITALLOCATE_PARR(self%prevQDiff, (self%nElem, self%mPrevVector))
    
    self%initMixParam = initMixParam
    
    if (present(tFromStart)) then
      self%tFromStart = tFromStart
    else
      self%tFromStart = .false.
    end if
    
  end subroutine DIISMixer_create
  
  !!* Makes the mixer ready for a new SCC cycle
  !!* @param self  DIIS mixer instance
  !!* @param nElem Nr. of elements in the vectors to mix
  subroutine DIISMixer_reset(self, nElem)
    type(ODIISMixer), pointer :: self
    integer, intent(in) :: nElem
    

    ASSERT(nElem > 0)
    
    if (nElem /= self%nElem) then
      self%nElem = nElem
      DEALLOCATE_PARR(self%prevQInput)
      DEALLOCATE_PARR(self%prevQDiff)
      ALLOCATE_PARR(self%prevQInput, (self%nElem, self%mPrevVector))
      ALLOCATE_PARR(self%prevQDiff, (self%nElem, self%mPrevVector))
    end if
    self%nPrevVector = 0
    self%indx = 0
    
  end subroutine DIISMixer_reset
  
  !!* Destroys the DIIS mixer
  !!* @param self DIIS mixer instance.
  subroutine DIISMixer_destroy(self)
    type(ODIISMixer), pointer :: self
    
    if (associated(self)) then
      DEALLOCATE_PARR(self%prevQInput)
      DEALLOCATE_PARR(self%prevQDiff)
    end if
    DEALLOCATE_P(self)
    
  end subroutine DIISMixer_destroy
  
  !!* Mixes charges according to the DIIS method
  !!* @param self       Pointer to the diis mixer
  !!* @param qInpResult Input charges on entry, mixed charges on exit.
  !!* @param qDiff      Charge difference
  subroutine DIISMixer_mix(self, qInpResult, qDiff)
    type(ODIISMixer), pointer :: self
    real(dp), intent(inout) :: qInpResult(:)
    real(dp), intent(in)    :: qDiff(:)
    
    real(dp), allocatable :: aa(:,:), bb(:,:)
    integer :: ii, jj
    
    ASSERT(size(qInpResult) == self%nElem)
    ASSERT(size(qDiff) == self%nElem)
    
    if (self%nPrevVector < self%mPrevVector) then
      self%nPrevVector = self%nPrevVector + 1
    end if
    
    call storeVectors(self%prevQInput, self%prevQDiff, self%indx, &
        &qInpResult, qDiff, self%mPrevVector)
    
    qInpResult(:) = 0.0_dp
    
    if (self%tFromStart .or. self%nPrevVector == self%mPrevVector) then
      
      ALLOCATE_(aa, (self%nPrevVector+1, self%nPrevVector+1))
      ALLOCATE_(bb, (self%nPrevVector+1, 1))
      
      aa(:,:) = 0.0_dp
      bb(:,:) = 0.0_dp
      
      do ii = 1, self%nPrevVector
        do jj = 1, self%nPrevVector
          aa(ii, jj) = dot_product( self%prevQDiff(:, ii), self%prevQDiff(:, jj) )
        end do
      end do
      aa(self%nPrevVector+1, 1:self%nPrevVector) = -1.0_dp
      aa(1:self%nPrevVector, self%nPrevVector+1) = -1.0_dp
      
      bb(self%nPrevVector+1,1) = -1.0_dp
      
      !! Solve system of linear equations
      call gesv(aa, bb)
      
      qInpResult(:) = 0.0_dp
      do ii = 1, self%nPrevVector
        qInpResult(:) = qInpResult(:) + bb(ii,1) * ( &
            & self%prevQInput(:,ii) + self%prevQDiff(:,ii) )
      end do
      
      DEALLOCATE_(aa)
      DEALLOCATE_(bb)
      
    end if
    
    if (self%nPrevVector < self%mPrevVector) then
      !! First few iterations return simple mixed vector
      qInpResult(:) = qInpResult(:) + self%initMixParam * qDiff(:)
    end if
    
  end subroutine DIISMixer_mix
  
  !!* Stores a vector pair in a limited storage. If the stack is full, oldest
  !!* vector pair is overwritten.
  !!* @param prevQInp    Contains previous vectors of the first type
  !!* @param prevQDiff   Contains previous vectors of the second type
  !!* @param indx        Indexing of data
  !!* @param qInput      New first vector
  !!* @param qDiff       New second vector
  !!* @param mPrevVector Size of the stacks.
  subroutine storeVectors(prevQInp, prevQDiff, indx, qInput, qDiff, &
      &mPrevVector)
    real(dp), intent(inout) :: prevQInp(:,:)
    real(dp), intent(inout) :: prevQDiff(:,:)
    integer, intent(inout) :: indx
    real(dp), intent(in) :: qInput(:)
    real(dp), intent(in) :: qDiff(:)
    integer, intent(in) :: mPrevVector
    
!    indx = indx + 1
    indx = mod(indx, mPrevVector) + 1
    prevQInp(:,indx) = qInput(:)
    prevQDiff(:,indx) = qDiff(:)
    
  end subroutine storeVectors
  
end module DIISMixer
