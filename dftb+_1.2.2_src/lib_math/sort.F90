!!* Various types of sorting routines, and related stuff
!!* @todo add other algorithms, radix? definitely not quicksort though
!!* adaptive heap sorts?
module sorting
#include "assert.h"
#include "allocate.h"
  use accuracy, only : dp
  implicit none
  private

  public :: heap_sort, index_heap_sort, unique
  

  !!* Heap sort algorithm - O(N log(N)) performance
  interface heap_sort
    module procedure heap_sort_real
    module procedure heap_sort_int
  end interface

  !!* Heap sort algorithm - O(N log(N)) performance, provides an index value
  !!* instead of re-ordering values
  interface index_heap_sort
    module procedure index_heap_sort_real
    module procedure index_heap_sort_int
  end interface

  !!* Function to count number of unique elements in a sorted array of value
  !!* greater than 0 and place them at the start of the array in order
  interface unique
    module procedure unique_int
  end interface unique

  
contains

  !!* real case in-place heap sort
  !!* @param array Array of values to be sorted
  !!* @param tolerance Tolerance for equality of two elements
  !!* @ref based on Numerical Recipes Software 1986-92
  subroutine heap_sort_real(array, tolerance)
    real(dp), intent(inout)        :: array(:)
    real(dp), intent(in), optional :: tolerance
    
    integer  :: n, ir, ij, il, ii, ik
    real(dp) :: tmpReal
    real(dp) :: tol

    if (present(tolerance)) then
      tol = tolerance
    else
      tol = epsilon(1.0_dp)
    end if

    n = size(array)
    if (n <= 1) return
    il = n/2 + 1
    ir = n
    ik = 1
    do while (ik == 1)
      if (il > 1) then
        il = il - 1
        tmpReal = array(il)
      else
        tmpReal = array(ir)
        array(ir) = array(1)
        ir = ir - 1
        if(ir == 1)then
          array(1) = tmpReal
          return
        end if
      end if
      ii = il
      ij = 2 * il
      do while (ij <= ir)
        if (ij < ir) then
          if (array(ij) < array(ij+1) - tol) then
            ij = ij + 1
          end if
        end if
        if(tmpReal < array(ij) - tol) then
          array(ii) = array(ij)
          ii = ij
          ij = 2 * ij
        else
          ij = ir + 1
        end if
      end do
      array(ii) = tmpReal
    end do

  end subroutine heap_sort_real


  
  !!* integer case in-place heap sort
  !!* @param array Array of values to be sorted
  !!* @ref based on Numerical Recipes Software 1986-92
  subroutine heap_sort_int(array)
    integer, intent(inout) :: array(:)

    integer :: n, ii, ir, ij, il, ik
    integer :: tmpInt
    
    n = size(array)
    if (n <= 1) return

    il = n/2 + 1
    ir = n
    ik = 1
    do while (ik == 1)
      if (il > 1) then
        il = il - 1
        tmpInt = array(il)
      else
        tmpInt = array(ir)
        array(ir) = array(1)
        ir = ir - 1
        if(ir == 1)then
          array(1) = tmpInt
          return
        end if
      end if
      ii = il
      ij = 2 * il
      do while (ij <= ir)
        if (ij < ir) then
          if (array(ij) < array(ij+1)) then
            ij = ij + 1
          end if
        end if
        if(tmpInt < array(ij)) then
          array(ii) = array(ij)
          ii = ij
          ij = 2 * ij
        else
          ij = ir + 1
        end if
      end do
      array(ii) = tmpInt
    end do

  end subroutine  heap_sort_int


  
  !!* Real case heap sort returning an index.
  !!* @param indx Indexing array on return
  !!* @param array Array of values to be sorted
  !!* @param tolerance Tolerance for equality of two elements
  !!* @ref based on Numerical Recipes Software 1986-92
  subroutine index_heap_sort_real(indx,array, tolerance)
    integer,  intent(out)          :: indx(:)
    real(dp), intent(in)           :: array(:)
    real(dp), intent(in), optional :: tolerance

    integer :: n, ir, ij, il, ii, ik
    integer :: indxTmp
    real(dp) :: arrayTmp, tol

    ASSERT(size(array)==size(indx))

    if (present(tolerance)) then
      tol = tolerance
    else
      tol = epsilon(0.0_dp)
    end if

    do ii=1,size(indx)
      indx(ii) = ii
    end do
    n = size(array)
    if (n <= 1) return
    il=n/2+1
    ir=n
    ik = 1
    do while (ik == 1)
      if (il.gt.1) then
        il=il-1
        indxTmp=indx(il)
        arrayTmp=array(indxTmp)
      else
        indxTmp=indx(ir)
        arrayTmp=array(indxTmp)
        indx(ir)=indx(1)
        ir=ir-1
        if (ir.lt.1) then
          indx(1)=indxTmp
          return
        end if
      end if
      ii=il
      ij=2 * il
      do while (ij <= ir)
        if (ij < ir) then
          if (array(indx(ij)) < array(indx(ij+1)) - tol) then
            ij = ij + 1
          end if
        end if
        if(arrayTmp < array(indx(ij)) - tol) then
          indx(ii)=indx(ij)
          ii=ij
          ij=2*ij
        else
          ij = ir + 1
        end if
      end do
      indx(ii)=indxTmp        
    end do
    
  end subroutine index_heap_sort_real


  
  !!* real case heap sort returning an index
  !!* @param indx Indexing array on return
  !!* @param array Array of values to be sorted
  !!* @ref based on Numerical Recipes Software 1986-92
  subroutine index_heap_sort_int(indx, array)
    integer, intent(out) :: indx(:)
    integer, intent(in) :: array(:)
    
    integer :: n, ir, ij, il, ii, ik
    integer :: indxTmp, arrayTmp

    ASSERT(size(array)==size(indx))

    do ii=1,size(indx)
      indx(ii) = ii
    end do
    n = size(array)
    if (n <= 1) return
    il=n/2+1
    ir=n
    ik = 1
    do while (ik == 1)
      if (il.gt.1) then
        il=il-1
        indxTmp=indx(il)
        arrayTmp=array(indxTmp)
      else
        indxTmp=indx(ir)
        arrayTmp=array(indxTmp)
        indx(ir)=indx(1)
        ir=ir-1
        if (ir.lt.1) then
          indx(1)=indxTmp
          return
        end if
      end if
      ii=il
      ij=2 * il
      do while (ij <= ir)
        if (ij < ir) then
          if (array(indx(ij)) < array(indx(ij+1))) then
            ij = ij + 1
          end if
        end if
        if(arrayTmp < array(indx(ij))) then
          indx(ii)=indx(ij)
          ii=ij
          ij=2*ij
        else
          ij = ir + 1
        end if
      end do
      indx(ii)=indxTmp        
    end do
    
  end subroutine index_heap_sort_int



  !!* Function to count number of unique elements in a sorted array of value
  !!* greater than 0 and place them at the start of the array in order.
  !!* @param array Array to make unique.
  !!* @param arraySize Constraints the effect of the subroutine on the first
  !!*   n elements, where n is the value for arraySize. (default: size(array))
  !!* @return Number of unique elements.
  !!* @todo check that the elements are in sorted order, and generalise for
  !!* decreasing order as well as increasing
  function unique_int(array, arraySize) result(nUnique)
    integer, intent(inout) :: array(:)
    integer, intent(in), optional :: arraySize
    integer :: nUnique

    integer :: ii, ij, nn

    if (present(arraySize)) then
      nn = arraySize
    else
      nn = size(array)
    end if

    ASSERT(nn >= 1 )
    ASSERT(nn <= size(array))
    ASSERT(all(array(:nn) > 0))
    
    ii = 1
    do ij = 2, nn
      if (array(ij) /= array(ii)) then
        ii = ii + 1
        array(ii) = array(ij)
      end if
    end do
    nUnique = ii
    
  end function unique_int


end module sorting
