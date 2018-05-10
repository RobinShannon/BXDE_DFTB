module TypeGeometry
# include "allocate.h"
  use accuracy
  implicit none
  private

  public :: TGeometry, destruct, normalize
  
  !!* Type for containing geometrical information
  type TGeometry
    integer           :: nAtom
    logical           :: tPeriodic
    logical           :: tFracCoord
    integer,  pointer :: species(:)  => null()
    real(dp), pointer :: coords(:,:) => null()
    integer           :: nSpecie
    real(dp), pointer :: origin(:) => null()
    real(dp), pointer :: latVecs(:,:) => null()
    real(dp), pointer :: recVecs2p(:,:) => null()
    character(mc), pointer :: specieNames(:) => null()
  end type TGeometry

  !!* Interface for destroying a TGeometry type
  interface destruct
    module procedure destruct_Geometry
  end interface

  !!* Interface for cleaning up a geometry against non-existent atom species
  interface normalize
    module procedure Geometry_normalize
  end interface


contains

  !!* Destroys the TGeometry object
  !!* @param self Geometry object
  subroutine destruct_Geometry(self)
    type(TGeometry), intent(inout) :: self

    DEALLOCATE_PARR(self%species)
    DEALLOCATE_PARR(self%coords)
    DEALLOCATE_PARR(self%origin)
    DEALLOCATE_PARR(self%latVecs)
    DEALLOCATE_PARR(self%recVecs2p)
    DEALLOCATE_PARR(self%specieNames)
    
  end subroutine destruct_Geometry

  !!* Normalises a geometry object to be safe against the abscence of any atoms
  !!* of a species specified in the input file
  !!* @param self Geometry object
  subroutine Geometry_normalize(sf)
    type(TGeometry), intent(inout) :: sf

    logical, allocatable :: inUse(:)
    integer, pointer :: oldSpecies(:)
    character(mc), pointer :: oldSpecieNames(:)
    integer :: ind, iSp

    ALLOCATE_(inUse, (sf%nSpecie)) 
    do iSp = 1, sf%nSpecie
      inUse(iSp) = any(sf%species == iSp)
    end do
    if (.not. all(inUse)) then !some of the species are redundant, so re-index
      oldSpecies => sf%species
      oldSpecieNames => sf%specieNames
      sf%nSpecie = count(inUse)
      INITALLOCATE_PARR(sf%species, (size(oldSpecies)))
      INITALLOCATE_PARR(sf%specieNames, (sf%nSpecie))
      ind = 1
      do iSp = 1, size(oldSpecieNames)
        if (.not. inUse(iSp)) then
          continue
        end if
        sf%specieNames(ind) = oldSpecieNames(iSp)
        where (oldSpecies == iSp)
          sf%species = ind
        end where
        ind = ind + 1
      end do
    end if
    DEALLOCATE_(inUse)
    
  end subroutine Geometry_normalize
 

end module TypeGeometry
