!!* Wrapper modules collecting the common interface of various dispersion
!!* models.
module Dispersion
#include "assert.h"
#include "allocate.h"
  use Accuracy
  use DispSlaterKirkwood
  use DispUFF
  use Periodic, only : TNeighborList
  implicit none
  private

  public :: ODispersion, TDispersionInp
  public :: init, destruct, getRCutoff, updateCoords, updateLatVecs
  public :: getEnergies, addGradients, getStress
  public :: iSlaterKirkwood, iVdWUFF

  !!* Dispersion type
  type ODispersion
    private
    integer :: iDisp
    type(ODispSlaKirk), pointer :: pSlaterKirkwood
    type(ODispUFF), pointer :: pVdWUFF
    logical :: tInit = .false.
  end type ODispersion

  !!* Input for various dispersion models.
  type TDispersionInp
    integer :: iDisp
    type(TDispSlaKirkInp), pointer :: pSlaKirkInp
    type(TDispUFFInp), pointer :: pVdWUFFInp
  end type TDispersionInp

  !!* Initialization interface.
  interface init
    module procedure Disp_initSlaterKirkwood
    module procedure Disp_initVdWUFF
  end interface

  !!* Destruction interface.
  interface destruct
    module procedure Disp_destruct
    module procedure DispInp_destruct
  end interface

  !!* Coordinate update.
  interface updateCoords
    module procedure Disp_updateCoords
  end interface

  !!* Lattice vector update
  interface updateLatVecs
    module procedure Disp_updateLatVecs
  end interface

  !!* Real space cutoff.
  interface getRCutoff
    module procedure Disp_getRCutoff
  end interface

  !!* Delivering dispersion energy.
  interface getEnergies
    module procedure Disp_getEnergies
  end interface

  !!* Adding dispersion gradient.
  interface addGradients
    module procedure Disp_addGradients
  end interface

  !!* Delivering dispersion stress tensor.
  interface getStress
    module procedure Disp_getStress
  end interface
  
  !! Enumerating available dispersion methods.
  integer, parameter :: iSlaterKirkwood = 1
  integer, parameter :: iVdWUFF = 2


contains

  !!* Initializes disperson with the SlaterKirkwood model.
  !!* @param sf Self.
  !!* @param pSlaterKirkwood Initialized SlaterKirkwood dispersion model.
  subroutine Disp_initSlaterKirkwood(sf, pSlaterKirkwood)
    type(ODispersion), intent(inout) :: sf
    type(ODispSlaKirk), pointer :: pSlaterKirkwood

    ASSERT(.not. sf%tInit)
    
    sf%iDisp = iSlaterKirkwood
    sf%pSlaterKirkwood => pSlaterKirkwood
    sf%tInit = .true.

  end subroutine Disp_initSlaterKirkwood


  !!* Initializes disperson with the UFF model.
  !!* @param sf Self.
  !!* @param pVdWUFF Initialized UFF dispersion model.
  subroutine Disp_initVdWUFF(sf, pVdWUFF)
    type(ODispersion), intent(out) :: sf
    type(ODispUFF), pointer :: pVdWUFF

    ASSERT(.not. sf%tInit)
    
    sf%iDisp = iVdWUFF
    sf%pVdWUFF => pVdWUFF
    sf%tInit = .true.

  end subroutine Disp_initVdWUFF


  !!* Destructor.
  !!* @param sf Self.
  subroutine Disp_destruct(sf)
    type(ODispersion), intent(inout) :: sf

    ASSERT(sf%tInit)

    select case(sf%iDisp)
    case(iSlaterKirkwood)
      call destruct(sf%pSlaterKirkwood)
      DEALLOCATE_P(sf%pSlaterKirkwood)
    case(iVdWUFF)
      call destruct(sf%pVdWUFF)
      DEALLOCATE_P(sf%pVdWUFF)
    end select
    sf%tInit = .false.

  end subroutine Disp_destruct


  !!* Updating the coordinates.
  !!* @param sf Self.
  !!* @param neigh Neighborlist
  !!* @param img2CentCell Mapping of atoms to image in central cell.
  !!* @param coords Current coordinates.
  !!* @param species0 Type of atoms in the central cell.
  subroutine Disp_updateCoords(sf, neigh, img2CentCell, coords, species0)
    type(ODispersion), intent(inout) :: sf
    type(TNeighborList), intent(in) :: neigh
    integer, intent(in) :: img2CentCell(:)
    real(dp), intent(in) :: coords(:,:)
    integer, intent(in) ::  species0(:)

    ASSERT(sf%tInit)
    
    select case(sf%iDisp)
    case(iSlaterKirkwood)
      call updateCoords(sf%pSlaterKirkwood, neigh, img2CentCell, coords)
    case(iVdWUFF)
      call updateCoords(sf%pVdWUFF, neigh, img2CentCell, coords, species0)
    end select
    
  end subroutine Disp_updateCoords


  !!* Updating the lattice vectors
  !!* @param sf Self.
  !!* @param latVecs  New lattice vectors
  !!* @param recVecs  New reciprocal vectors
  !!* @param vol  New unit cell volume
  !!* @param species0  Species in the unit cell.
  subroutine Disp_updateLatVecs(sf, latVecs, recVecs, vol, species0)
    type(ODispersion), intent(inout) :: sf
    real(dp), intent(in) :: latVecs(:,:), recVecs(:,:)
    real(dp), intent(in) :: vol
    integer, intent(in) :: species0(:)

    ASSERT(sf%tInit)
    
    select case(sf%iDisp)
    case(iSlaterKirkwood)
      call updateLatVecs(sf%pSlaterKirkwood, latVecs, recVecs, vol)
    case(iVdWUFF)
      call updateLatVecs(sf%pVdWUFF, latVecs, recVecs, vol, species0)
    end select
    
  end subroutine Disp_updateLatVecs


  !!* Returns real space cutoff for Ewald-like summation.
  !!* @param sf Self.
  !!* @return Real space cutoff.
  function Disp_getRCutoff(sf) result(cutoff)
    type(ODispersion), intent(in) :: sf
    real(dp) :: cutoff

    ASSERT(sf%tInit)
    
    select case(sf%iDisp)
    case(iSlaterKirkwood)
      cutoff = getRCutoff(sf%pSlaterKirkwood)
    case(iVdWUFF)
      cutoff =  getRCutoff(sf%pVdWUFF)
    end select
    
  end function Disp_getRCutoff


  !!* Returns energy per atom due to dispersion.
  !!* @param sf Self.
  !!* @param energies Energy per atom on return.
  subroutine Disp_getEnergies(sf, energies)
    type(ODispersion), intent(in) :: sf
    real(dp), intent(out) :: energies(:)

    ASSERT(sf%tInit)
    
    select case(sf%iDisp)
    case(iSlaterKirkwood)
      call getEnergies(sf%pSlaterKirkwood, energies)
    case(iVdWUFF)
      call getEnergies(sf%pVdWUFF, energies)
    end select

  end subroutine Disp_getEnergies

  
  !!* Adds gradients due to dispersion.
  !!* @param sf Self.
  !!* @param gradients Adds contribution to the gradients due to the
  !!* dispersion model.
  subroutine Disp_addGradients(sf, gradients)
    type(ODispersion), intent(in) :: sf
    real(dp), intent(inout) :: gradients(:,:)

    ASSERT(sf%tInit)
    
    select case(sf%iDisp)
    case(iSlaterKirkwood)
      call addGradients(sf%pSlaterKirkwood, gradients)
    case(iVdWUFF)
      call addGradients(sf%pVdWUFF, gradients)
    end select

  end subroutine Disp_addGradients
  
  !!* Gets stress tensor due to dispersion.
  !!* @param sf Self.
  !!* @param stress tensor from the used dispersion model.
  subroutine Disp_getStress(sf, stress)
    type(ODispersion), intent(in) :: sf
    real(dp), intent(out) :: stress(3,3)
    
    ASSERT(sf%tInit)
    
    select case(sf%iDisp)
    case(iSlaterKirkwood)
!      call getStress(sf%pSlaterKirkwood, stress)
    case(iVdWUFF)
      call getStress(sf%pVdWUFF, stress)
    end select
    
  end subroutine Disp_getStress
  
  !!* Destructor for the dispersion input.
  !!* @param sf Self.
  subroutine DispInp_destruct(sf)
    type(TDispersionInp), intent(inout) :: sf

    select case(sf%iDisp)
    case(iSlaterKirkwood)
      call destruct(sf%pSlaKirkInp)
      DEALLOCATE_P(sf%pSlaKirkInp)
    case(iVdWUFF)
      call destruct(sf%pVdWUFFInp)
      DEALLOCATE_P(sf%pVdWUFFInp)
    end select

  end subroutine DispInp_destruct


end module Dispersion
