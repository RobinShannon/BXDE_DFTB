!!* Contains a general mixer which calls the desired real mixers.
module Mixer
# include "allocate.h"
# include "assert.h"
  use accuracy
  use SimpleMixer
  use AndersonMixer
  use BroydenMixer
  use DIISMixer
  implicit none
  
  private
  
  !!* Interface type for various mixers.
  type OMixer
    private
    integer :: mixerType
    type(OSimpleMixer),   pointer :: pSimpleMixer
    type(OAndersonMixer), pointer :: pAndersonMixer
    type(OBroydenMixer),  pointer :: pBroydenMixer
    type(ODIISMixer),  pointer :: pDIISMixer
  end type OMixer

  !!* Creates mixer
  interface create
    module procedure Mixer_createSimple
    module procedure Mixer_createAnderson
    module procedure Mixer_createBroyden
    module procedure Mixer_createDIIS
  end interface

  !!* Destroys mixer
  interface destroy
    module procedure Mixer_destroy
  end interface

  !!* Resets mixer
  interface reset
    module procedure Mixer_reset
  end interface

  !!* Does actual mixing
  interface mix
    module procedure Mixer_mix
    module procedure Mixer_mixR2
    module procedure Mixer_mixR3
    module procedure Mixer_mixR4
  end interface

  
  public :: OMixer
  public :: create, destroy, reset, mix

  
  !! Constants for the different mixers.
  integer, parameter :: iSimpleMixer = 1
  integer, parameter :: iAndersonMixer = 2
  integer, parameter :: iBroydenMixer = 3
  integer, parameter :: iDIISMixer = 4

contains

  !!* Common stuff for initializing a Mixer.
  !!* @param self Initialized mixer on exit.
  subroutine Mixer_createCommon(self)
    type(OMixer), pointer :: self

    INITALLOCATE_P(self)
    self%pSimpleMixer => null()
    self%pAndersonMixer => null()
    self%pBroydenMixer => null()
    self%pDIISMixer => null()

  end subroutine Mixer_createCommon

  

  !!* Initializes a mixer as simple mixer.
  !!* @param self Mixer instance
  !!* @param pSimple Pointer to a valid simple mixer instance.
  subroutine Mixer_createSimple(self, pSimple)
    type(OMixer), pointer :: self
    type(OSimpleMixer), pointer :: pSimple

    call Mixer_createCommon(self)
    self%mixerType = iSimpleMixer
    self%pSimpleMixer => pSimple

  end subroutine Mixer_createSimple

  

  !!* Initializes a mixer as Anderson mixer.
  !!* @param self Mixer instance
  !!* @param pAnderson Pointer to a valid Anderson mixer instance.
  subroutine Mixer_createAnderson(self, pAnderson)
    type(OMixer), pointer :: self
    type(OAndersonMixer), pointer :: pAnderson

    call Mixer_createCommon(self)
    self%mixerType = iAndersonMixer
    self%pAndersonMixer => pAnderson

  end subroutine Mixer_createAnderson


  
  !!* Initializes a mixer as Broyden mixer
  !!* @param self Mixer instance
  !!* @param pBroyed Pointer to a valid Broyden mixer instance.
  subroutine Mixer_createBroyden(self, pBroyden)
    type(OMixer), pointer :: self
    type(OBroydenMixer), pointer :: pBroyden

    call Mixer_createCommon(self)
    self%mixerType = iBroydenMixer
    self%pBroydenMixer => pBroyden

  end subroutine Mixer_createBroyden
  
  !!* Initializes a mixer as DIIS mixer
  !!* @param self Mixer instance
  !!* @param pDIIS Pointer to a valid DIIS mixer instance.
  subroutine Mixer_createDIIS(self, pDIIS)
    type(OMixer), pointer :: self
    type(ODIISMixer), pointer :: pDIIS
    
    call Mixer_createCommon(self)
    self%mixerType = iDIISMixer
    self%pDIISMixer => pDIIS
    
  end subroutine Mixer_createDIIS
  
  !!* Resets the mixer
  !!* @param self  Mixer instance.
  !!* @param nElem Size of the vectors to mix.
  subroutine Mixer_reset(self, nElem)
    type(OMixer), pointer :: self
    integer, intent(in) :: nElem

    select case (self%mixerType)
    case(iSimpleMixer)
      call reset(self%pSimpleMixer, nElem)
    case (iAndersonMixer)
      call reset(self%pAndersonMixer, nElem)
    case (iBroydenMixer)
      call reset(self%pBroydenMixer, nElem)
    case (iDIISMixer)
      call reset(self%pDIISMixer, nElem)
    end select

  end subroutine Mixer_reset


  
  !!* Destroys the mixer
  !!* @param self Mixer instance
  subroutine Mixer_destroy(self)
    type(OMixer), pointer :: self
    
    if (associated(self)) then
      select case (self%mixerType)
      case (iSimpleMixer)
        call destroy(self%pSimpleMixer)
      case (iAndersonMixer)
        call destroy(self%pAndersonMixer)
      case (iBroydenMixer)
        call destroy(self%pBroydenMixer)
      case (iDIISMixer)
        call destroy(self%pDIISMixer)
      end select
    end if
    DEALLOCATE_P(self)
    
  end subroutine Mixer_destroy
  

  
  !!* Mixes to vectors together
  !!* @param self Mixer instance.
  !!* @param qInpRes Input vector on entry, result vector on exit.
  !!* @param qDiff   Difference vector
  subroutine Mixer_mix(self, qInpRes, qDiff)
    type(OMixer), pointer    :: self
    real(dp),      intent(inout) :: qInpRes(:)
    real(dp),      intent(in) :: qDiff(:)

    select case (self%mixerType)
    case (iSimpleMixer)
      call mix(self%pSimpleMixer, qInpRes, qDiff)
    case (iAndersonMixer)
      call mix(self%pAndersonMixer, qInpRes, qDiff)
    case (iBroydenMixer)
      call mix(self%pBroydenMixer, qInpRes, qDiff)
    case (iDIISMixer)
      call mix(self%pDIISMixer, qInpRes, qDiff)
    end select
    
  end subroutine Mixer_mix



  !!* Wrapper for mixing real rank 2 arrays (e.g. density matrix)
  !!* @param self Mixer instance.
  !!* @param qInpRes Input array on entry, result vector on exit.
  !!* @param qDiff   Difference array
  subroutine Mixer_mixR2(self, qInpRes, qDiff)
    type(OMixer), pointer :: self
    real(dp), intent(inout) :: qInpRes(:,:)
    real(dp), intent(in) :: qDiff(:,:)

    real(dp), allocatable :: qInpResR1(:), qDiffR1(:)

    ALLOCATE_(qInpResR1, (size(qInpRes)))
    ALLOCATE_(qDiffR1, (size(qDiff)))
    qInpResR1(:) = reshape(qInpRes, (/ size(qInpRes) /))
    qDiffR1(:) = reshape(qDiff, (/ size(qDiff) /))
    call mix(self, qInpResR1, qDiffR1)
    qInpRes(:,:) = reshape(qInpResR1, shape(qInpRes))
    DEALLOCATE_(qInpResR1)
    DEALLOCATE_(qDiffR1)

  end subroutine Mixer_mixR2



  !!* Wrapper for mixing real rank 3 arrays
  !!* @param self Mixer instance.
  !!* @param qInpRes Input array on entry, result vector on exit.
  !!* @param qDiff   Difference array
  subroutine Mixer_mixR3(self, qInpRes, qDiff)
    type(OMixer), pointer :: self
    real(dp), intent(inout) :: qInpRes(:,:,:)
    real(dp), intent(in) :: qDiff(:,:,:)

    real(dp), allocatable :: qInpResR1(:), qDiffR1(:)
    

    ALLOCATE_(qInpResR1, (size(qInpRes)))
    ALLOCATE_(qDiffR1, (size(qDiff)))
    qInpResR1(:) = reshape(qInpRes, (/ size(qInpRes) /))
    qDiffR1(:) = reshape(qDiff, (/ size(qDiff) /))
    call mix(self, qInpResR1, reshape(qDiff, (/size(qDiff)/)))
    qInpRes(:,:,:) = reshape(qInpResR1, shape(qInpRes))
    DEALLOCATE_(qInpResR1)
    DEALLOCATE_(qDiffR1)
    
  end subroutine Mixer_mixR3

  !!* Wrapper for mixing real rank 3 arrays
  !!* @param self Mixer instance.
  !!* @param qInpRes Input array on entry, result vector on exit.
  !!* @param qDiff   Difference array
  subroutine Mixer_mixR4(self, qInpRes, qDiff)
    type(OMixer), pointer :: self
    real(dp), intent(inout) :: qInpRes(:,:,:,:)
    real(dp), intent(in) :: qDiff(:,:,:,:)

    real(dp), allocatable :: qInpResR1(:), qDiffR1(:)
    
    
    ASSERT(all(shape(qDiff)==shape(qInpRes)))

    ALLOCATE_(qInpResR1, (size(qInpRes)))
    ALLOCATE_(qDiffR1, (size(qDiff)))
    qInpResR1(:) = reshape(qInpRes, (/ size(qInpRes) /))
    qDiffR1(:) = reshape(qDiff, (/ size(qDiff) /))
    call mix(self, qInpResR1, reshape(qDiff, (/size(qDiff)/)))
    qInpRes(:,:,:,:) = reshape(qInpResR1, shape(qInpRes))
    DEALLOCATE_(qInpResR1)
    DEALLOCATE_(qDiffR1)
    
  end subroutine Mixer_mixR4
  
end module Mixer
