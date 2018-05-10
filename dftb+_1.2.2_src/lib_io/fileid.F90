module FileId
  use message
  
  private

  public :: getFileId
  
  integer, parameter :: minId = 20
  integer, parameter :: maxId = 65535
  
contains

  function getFileId()
    integer :: getFileId

    integer, save :: curId = minId

    if (curId > maxId) then
      call error("getFileId: No more free file identification numbers")
    end if
    getFileId = curId
    curId = curId + 1
    
  end function getFileId

end module FileId
