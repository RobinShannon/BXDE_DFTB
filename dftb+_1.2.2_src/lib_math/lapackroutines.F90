!!* Contains F90 wrapper functions for some commonly used lapack calls needed
!!* in the code
!!* @todo Proper documentation, and add more lapack calls
module lapackRoutines
#include "assert.h"
#include "allocate.h"
  use accuracy
  use message
  implicit none

  private
  
  character(len=100) :: error_string  !* Used to return runtime diagnostics

  !!* Computes the solution to a real system of linear equations
  !!* A * X = B,
  !!* where A is an N-by-N matrix and X and B are N-by-NRHS matrices
  !!* @param aa  Contains the coefficients on entry, the LU factorisation
  !!*   on exit.
  !!* @param bb  Right hand side(s) of the linear equation on entry, solution(s)
  !!*   on exit.
  !!* @param nEquation  The size of the problem (nr. of variables and
  !!*    equations). Must be only specified if different from size(aa, dim=1).
  !!* @param nSolution  Nr. of right hand sides (nr. of solutions). Must be only
  !!*   specified if different from size(b, dim=2).
  !!* @param iError  Error flag. If present, Lapack error flags are reported
  !!*   and noncritical errors (iError > 0) will not abort the program.
  interface gesv
    module procedure gesv_real
    module procedure gesv_dble
  end interface gesv

  !!* Computes the LU decomposition of a general rectangular matrix using
  !!* partial pivoting with row interchanges.
  !!* @param aa      Matrix to decompose on entry, L and U on exit. Unit
  !!*   diagonal elements of L are not stored.
  !!* @param ipiv    Pivot indices, row i of the matrix was interchanged with
  !!*   row ipiv(i).
  !!* @param nRow    Number of rows of the matrix to decompose. (Necessary if
  !!*   different from the number of rows of the passed matrix)
  !!* @param nColumn Number of rows of the matrix to decompose. (Necessary if
  !!*   different from the number of columns of the passed matrix)
  !!* @param iError  Error flag. Zero on successfull exit. If not present, any
  !!*   lapack error causes program termination. If passed only fatal lapack
  !!*   errors with error flag < 0 cause abort.
  !!* @desc
  !!*   The decomposition has the form: A = P*L*U, where P is a permutation
  !!*   matrix, L is a lower triangular matrix with unit diagonal elements and
  !!*   U is an upper triangular matrix.
  interface getrf
    module procedure getrf_real
    module procedure getrf_dble
  end interface

  !!* Computes the inverse of a matrix using LU factorization computed by getrf.
  !!* @param aa      Matrix to decompose on entry, L and U on exit. Unit
  !!*   diagonal elements of L are not stored.
  !!* @param ipiv    Pivot indices, as calculated by getri
  !!* @param nRow    Number of rows of the matrix to decompose. (Necessary if
  !!*   different from the number of rows of the passed matrix)
  !!* @param iError  Error flag. Zero on successfull exit. If not present, any
  !!*   lapack error causes program termination. If present, only fatal lapack
  !!*   errors with error flag < 0 cause abort.
  interface getri
    module procedure getri_real
    module procedure getri_dble
  end interface

  !!* Solves a system of linear equations A*X = B with a real
  !!* symmetric matrix A using the factorization A = U*D*U**T or
  !!* A = L*D*L**T computed by DSYTRF.
  !!* @param A on entry, the symmetric matrix A.  If UPLO = 'U', the leading
  !!*   N-by-N upper triangular part of A contains the upper
  !!*   triangular part of the matrix A, and the strictly lower
  !!*   triangular part of A is not referenced.  If UPLO = 'L', the
  !!*   leading N-by-N lower triangular part of A contains the lower
  !!*   triangular part of the matrix A, and the strictly upper
  !!*   triangular part of A is not referenced.
  !!*   On exit, the block diagonal matrix D and the multipliers used
  !!*   to obtain the factor U or L
  !!* @param B On entry, the right hand side matrix B. On exit, the solution
  !!*   matrix X.
  !!* @param nRow  Number of rows of the matrix to decompose. (Necessary if
  !!*   different from the number of rows of the passed matrix)
  !!* @param uplo upper or lower triangle of the matrix, defaults to lower
  !!* @param iError Error flag. Zero on successfull exit. If not present, any
  !!*   lapack error causes program termination. If present, only fatal lapack
  !!*   errors with error flag < 0 cause abort.
  interface sytrs
    module procedure sytrs_dble
    module procedure sytrs_real
  end interface
  
  !!* returns a vector of random numers, either from a uniform or normal
  !!* distribution
  !!* @param iDist choice of distribution (1: uniform (0,1), 2: uniform (-1,1),
  !!* 3: normal (0,1)
  !!* @param iSeed INTEGER array, dimension (4)
  !!* On entry, the seed of the random number generator; the array
  !!* elements must be between 0 and 4095, and ISEED(4) must be
  !!* odd.
  !!* On exit, the seed is updated.
  !!* @param x On exit, vector of random numbers
  interface larnv
    module procedure larnv_real
    module procedure larnv_dble
    module procedure larnv_cplx
    module procedure larnv_dblecplx
  end interface larnv
  
  public :: gesv, getri, getrf, matinv, sytrs, larnv


contains

  ! Single precision version of gesv
  subroutine gesv_real(aa, bb, nEquation, nSolution, iError)
    real(rsp), intent(inout) :: aa(:,:)
    real(rsp), intent(inout) :: bb(:,:)
    integer, intent(in), optional :: nEquation
    integer, intent(in), optional :: nSolution
    integer, intent(out), optional :: iError

    integer :: info
    integer :: nn, nrhs, lda, ldb
    integer, allocatable :: ipiv(:)

    lda = size(aa, dim=1)
    if (present(nEquation)) then
      ASSERT(nEquation >= 1 .and. nEquation <= lda)
      nn = nEquation
    else
      nn = lda
    end if
    ASSERT(size(aa, dim=2) >= nn)

    ldb = size(bb, dim=1)
    ASSERT(ldb >= nn)
    nrhs = size(bb, dim=2)
    if (present(nSolution)) then
      ASSERT(nSolution <= nrhs)
      nrhs = nSolution
    end if
    
    info = 0
    ALLOCATE_(ipiv, (nn))
    call sgesv(nn, nrhs, aa, lda, ipiv, bb, ldb, info)
    DEALLOCATE_(ipiv)
    
    if (info < 0) then
99000 format ('Failure in linear equation solver sgesv,', &
          & ' illegal argument at position ',i10)
      write (error_string, 99000) info
      call error(error_string)
    else
      if (present(iError)) then
        iError = info
      elseif (info > 0) then
99010 format ('Linear dependent system in linear equation solver sgetrf,', &
          & ' info flag: ',i10)
        write (error_string, 99010) info
        call error(error_string)
      end if
    end if

  end subroutine gesv_real



  ! Double precision version of gesv
  subroutine gesv_dble(aa, bb, nEquation, nSolution, iError)
    real(rdp), intent(inout) :: aa(:,:)
    real(rdp), intent(inout) :: bb(:,:)
    integer, intent(in), optional :: nEquation
    integer, intent(in), optional :: nSolution
    integer, intent(out), optional :: iError

    integer :: info
    integer :: nn, nrhs, lda, ldb
    integer, allocatable :: ipiv(:)

    lda = size(aa, dim=1)
    if (present(nEquation)) then
      ASSERT(nEquation >= 1 .and. nEquation <= lda)
      nn = nEquation
    else
      nn = lda
    end if
    ASSERT(size(aa, dim=2) >= nn)

    ldb = size(bb, dim=1)
    ASSERT(ldb >= nn)
    nrhs = size(bb, dim=2)
    if (present(nSolution)) then
      ASSERT(nSolution <= nrhs)
      nrhs = nSolution
    end if
    
    info = 0
    ALLOCATE_(ipiv, (nn))
    call dgesv(nn, nrhs, aa, lda, ipiv, bb, ldb, info)
    DEALLOCATE_(ipiv)
    
    if (info < 0) then
99020 format ('Failure in linear equation solver dgesv,', &
          & ' illegal argument at position ',i10)
      write (error_string, 99020) info
      call error(error_string)
    else
      if (present(iError)) then
        iError = info
      elseif (info > 0) then
99030 format ('Linear dependent system in linear equation solver dgesv,', &
          & ' info flag: ',i10)
        write (error_string, 99030) info
        call error(error_string)
      end if
    end if

  end subroutine gesv_dble



  !!* Single precision version of getrf.
  subroutine getrf_real(aa, ipiv, nRow, nColumn, iError)
    real(rsp), intent(inout) :: aa(:,:)
    integer, intent(out) :: ipiv(:)
    integer, optional, intent(in) :: nRow
    integer, optional, intent(in) :: nColumn
    integer, optional, intent(out) :: iError
    
    integer :: mm, nn, lda, info

    lda = size(aa, dim=1)
    nn = size(aa, dim=2)
    if (present(nRow)) then
      ASSERT(nRow >= 1 .and. nRow <= lda)
      mm = nRow
    else
      mm = lda
    end if
    if (present(nColumn)) then
      ASSERT(nColumn >= 1 .and. nColumn <= nn)
      nn = nColumn
    end if
    ASSERT(size(ipiv) == min(mm, nn))

    call sgetrf(mm, nn, aa, lda, ipiv, info)

    if (info < 0) then
99040 format ('Failure in LU factorisation sgetrf,', &
          & ' illegal argument at position ',i10)
      write (error_string, 99040) info
      call error(error_string)
    else
      if (present(iError)) then
        iError = info
      elseif (info > 0) then
99050 format ('Factor U is exactly zero in sgetrf,', &
          & ' info flag is ',i10)
        write (error_string, 99050) info
        call error(error_string)
      end if
    end if
    
  end subroutine getrf_real



  !!* Double precision version of getrf.
  subroutine getrf_dble(aa, ipiv, nRow, nColumn, iError)
    real(rdp), intent(inout) :: aa(:,:)
    integer, intent(out) :: ipiv(:)
    integer, optional, intent(in) :: nRow
    integer, optional, intent(in) :: nColumn
    integer, optional, intent(out) :: iError
    
    integer :: mm, nn, lda, info

    lda = size(aa, dim=1)
    nn = size(aa, dim=2)
    if (present(nRow)) then
      ASSERT(nRow >= 1 .and. nRow <= lda)
      mm = nRow
    else
      mm = lda
    end if
    if (present(nColumn)) then
      ASSERT(nColumn >= 1 .and. nColumn <= nn)
      nn = nColumn
    end if
    ASSERT(size(ipiv) == min(mm, nn))

    call dgetrf(mm, nn, aa, lda, ipiv, info)

    if (info < 0) then
99060 format ('Failure in LU factorisation dgetrf,', &
          & ' illegal argument at position ',i10)
      write (error_string, 99060) info
      call error(error_string)
    else
      if (present(iError)) then
        iError = info
      elseif (info > 0) then
99070 format ('Factor U is exactly zero in dgetrf,', &
          & ' info flag is ',i10)
        write (error_string, 99070) info
        call error(error_string)
      end if
    end if
    
  end subroutine getrf_dble



  !!* Single precision version of getri.
  subroutine getri_real(aa, ipiv, nRow, iError)
    real(rsp), intent(inout) :: aa(:,:)
    integer, intent(in) :: ipiv(:)
    integer, intent(in), optional :: nRow
    integer, intent(out), optional :: iError

    integer :: nn, lda, info, lwork
    real(rsp), allocatable :: work(:)
    real(rsp) :: work2(1)
    
    lda = size(aa, dim=1)
    if (present(nRow)) then
      ASSERT(nRow >= 1 .and. nRow <= lda)
      nn = nRow
    else
      nn = lda
    end if
    ASSERT(size(aa, dim=2) >= nn)
    ASSERT(size(ipiv) == nn)

    lwork = -1
    call sgetri(nn, aa, lda, ipiv, work2, lwork, info)
    lwork = int(work2(1))

    ALLOCATE_(work, (lwork))
    call sgetri(nn, aa, lda, ipiv, work, lwork, info)
    DEALLOCATE_(work)

    if (info < 0) then
99080 format ('Failure in LU factorisation (sgetri),', &
          & ' illegal argument at position ',i10)
      write (error_string, 99080) info
      call error(error_string)
    else
      if (present(iError)) then
        iError = info
      elseif (info > 0) then
99090 format ('Factor U is exactly zero in sgetri,', &
          & ' info flag is ',i10)
        write (error_string, 99090) info
        call error(error_string)
      end if
    end if

  end subroutine getri_real
  


  !!* Double precision version of getri.
  subroutine getri_dble(aa, ipiv, nRow, iError)
    real(rdp), intent(inout) :: aa(:,:)
    integer, intent(in) :: ipiv(:)
    integer, intent(in), optional :: nRow
    integer, intent(out), optional :: iError

    integer :: nn, lda, info, lwork
    real(rdp), allocatable :: work(:)
    real(rdp) :: work2(1)
    
    lda = size(aa, dim=1)
    if (present(nRow)) then
      ASSERT(nRow >= 1 .and. nRow <= lda)
      nn = nRow
    else
      nn = lda
    end if
    ASSERT(size(aa, dim=2) >= nn)
    ASSERT(size(ipiv) == nn)

    lwork = -1
    call dgetri(nn, aa, lda, ipiv, work2, lwork, info)
    lwork = int(work2(1))

    ALLOCATE_(work, (lwork))
    call dgetri(nn, aa, lda, ipiv, work, lwork, info)
    DEALLOCATE_(work)

    if (info < 0) then
99100 format ('Failure in LU factorisation (dgetri), illegal argument at&
          & position ', i10)
      write (error_string, 99100) info
      call error(error_string)
    else
      if (present(iError)) then
        iError = info
      elseif (info > 0) then
99110   format ('Factor U is exactly zero in dgetri,', &
            & ' info flag is ',i10)
        write (error_string, 99110) info
        call error(error_string)
      end if
    end if

  end subroutine getri_dble
  


  !!* Inverts a matrix.
  !!* @param aa     Matrix to invert on entry, inverted matrix on exit
  !!* @param nRow   Nr. of rows of the matrix (if different from size(aa, dim=1)
  !!* @param iError Error flag. Returns 0 on successfull operation. If this
  !!*   variable is not specified, any occuring error (e.g. singular matrix)
  !!*   stops the program.
  subroutine matinv(aa, nRow, iError)
    real(dp), intent(inout) :: aa(:,:)
    integer, intent(in), optional :: nRow
    integer, intent(out), optional :: iError

    integer :: nn, info
    integer, allocatable :: ipiv(:)

    nn = size(aa, dim=1)
    if (present(nRow)) then
      ASSERT(nRow >= 1 .and. nRow <= nn)
      nn = nRow
    end if
    ASSERT(size(aa, dim=2) >= nn)
    
    ALLOCATE_(ipiv, (nn))
    call getrf(aa, ipiv, nRow=nn, nColumn=nn, iError=info)
    if (info == 0) then
      call getri(aa, ipiv, nRow=nn, iError=info)
    end if
    DEALLOCATE_(ipiv)

    if (present(iError)) then
      iError = info
    elseif (info /= 0) then
99120 format ('Matrix inversion failed because of error in getrf or getri.', &
          & ' Info flag: ',i10)
      write (error_string, 99120) info
      call error(error_string)
    end if
    
  end subroutine matinv

  !!* Single precision version of sytrs
  subroutine sytrs_real(A,B, nRow, uplo,iError)
    real(rsp), intent(inout)        :: A(:,:)
    real(rsp), intent(inout)        :: B(:,:)
    integer, intent(in), optional   :: nRow
    character, intent(in), optional :: uplo
    integer, intent(out), optional  :: iError
    
    integer, allocatable :: ipiv(:)
    character :: iUplo
    integer :: nn, lda, ldb, info, lwork, nrhs
    real(rsp), allocatable :: work(:)
    real(rsp) :: work2(1)

    lda = size(A, dim=1)
    if (present(nRow)) then
      ASSERT(nRow >= 1 .and. nRow <= lda)
      nn = nRow
    else
      nn = lda
    end if
    
    ldb = size(b, dim=1)
    ASSERT(ldb >= nn)
    
    if (present(uplo)) then
       iUplo = uplo
    else
       iUplo = 'L'
    end if
    ASSERT(iUplo == 'u' .or. iUplo == 'U' .or. iUplo == 'l' .or. iUplo == 'L')

    ASSERT(size(A, dim=2) >= nn)
    nrhs = size(B, dim=2)
    
    ALLOCATE_(ipiv, (nn))
    
    lwork = -1
    call ssytrf(iUplo, nn, A, lda, ipiv, work2, lwork, info)
    lwork = int(work2(1))
    if (info == 0) then
      ALLOCATE_(work, (lwork))
      call ssytrf(iUplo, nn, A, lda, ipiv, work, lwork, info)
      DEALLOCATE_(work)
    end if

    if (info == 0) then
       call ssytrs(iUplo, nn, nrhs, A, lda, ipiv, B, ldb, info)
     end if
     
    DEALLOCATE_(ipiv)
    
    if (present(iError)) then
       iError = info
    elseif (info /= 0) then
99130  format ('Solution failed because of error in sytrf or sytrs.',&
           & ' Info flag: ',i10)
       write (error_string, 99130) info
       call error(error_string)
    end if
        
  end subroutine sytrs_real

  !!* Double precision version of sytrs
  subroutine sytrs_dble(A,B, nRow, uplo,iError)
    real(rdp), intent(inout)        :: A(:,:)
    real(rdp), intent(inout)        :: B(:,:)
    integer, intent(in), optional   :: nRow
    character, intent(in), optional :: uplo
    integer, intent(out), optional  :: iError
    
    integer, allocatable :: ipiv(:)
    character :: iUplo
    integer :: nn, lda, ldb, info, lwork, nrhs
    real(rdp), allocatable :: work(:)
    real(rdp) :: work2(1)

    lda = size(A, dim=1)
    if (present(nRow)) then
      ASSERT(nRow >= 1 .and. nRow <= lda)
      nn = nRow
    else
      nn = lda
    end if
    
    ldb = size(b, dim=1)
    ASSERT(ldb >= nn)
    
    if (present(uplo)) then
       iUplo = uplo
    else
       iUplo = 'L'
    end if
    ASSERT(iUplo == 'u' .or. iUplo == 'U' .or. iUplo == 'l' .or. iUplo == 'L')

    ASSERT(size(A, dim=2) >= nn)
    nrhs = size(B, dim=2)
    
    ALLOCATE_(ipiv, (nn))
    
    lwork = -1
    call dsytrf(iUplo, nn, A, lda, ipiv, work2, lwork, info)
    lwork = int(work2(1))
    if (info == 0) then
      ALLOCATE_(work, (lwork))
      call dsytrf(iUplo, nn, A, lda, ipiv, work, lwork, info)
      DEALLOCATE_(work)
    end if

    if (info == 0) then
       call dsytrs(iUplo, nn, nrhs, A, lda, ipiv, B, ldb, info)
     end if
     
    DEALLOCATE_(ipiv)
    
    if (present(iError)) then
       iError = info
    elseif (info /= 0) then
99130  format ('Solution failed because of error in sytrf or sytrs.',&
           & ' Info flag: ',i10)
       write (error_string, 99130) info
       call error(error_string)
    end if
        
  end subroutine sytrs_dble
  
  !!* single precision version of larnv
  subroutine larnv_real(iDist,iSeed,x)
    integer, intent(in)    :: iDist
    integer, intent(inout) :: iSeed(4)
    real(rsp), intent(out) :: x(:)

    integer :: n

    ASSERT(iDist > 0)
    ASSERT(iDist < 4)
    ASSERT(all(iSeed(:) >= 0))
    ASSERT(all(iSeed(:) <= 4095))
    ASSERT(mod(iSeed(4),2) == 1)
    ASSERT(size(x) > 0)    
    n = size(x)
    x(:) = 0.0
    call SLARNV( iDist, iSeed, n, x )
  end subroutine larnv_real

  !!* double precision version of larnv
  subroutine larnv_dble(iDist,iSeed,x)
    integer, intent(in)    :: iDist
    integer, intent(inout) :: iSeed(4)
    real(rdp), intent(out) :: x(:)

    integer :: n

    ASSERT(iDist > 0)
    ASSERT(iDist < 4)
    ASSERT(all(iSeed(:) >= 0))
    ASSERT(all(iSeed(:) <= 4095))
    ASSERT(mod(iSeed(4),2) == 1)
    ASSERT(size(x) > 0)    
    n = size(x)
    x(:) = 0.0d0
    call DLARNV( iDist, iSeed, n, x )
  end subroutine larnv_dble

  !!* complex version of larnv
  subroutine larnv_cplx(iDist,iSeed,x)
    integer, intent(in)       :: iDist
    integer, intent(inout)    :: iSeed(4)
    complex(rsp), intent(out) :: x(:)

    integer :: n

    ASSERT(iDist > 0)
    ASSERT(iDist < 4)
    ASSERT(all(iSeed(:) >= 0))
    ASSERT(all(iSeed(:) <= 4095))
    ASSERT(mod(iSeed(4),2) == 1)
    ASSERT(size(x) > 0)    
    n = size(x)
    x(:) = 0.0
    call CLARNV( iDist, iSeed, n, x )
  end subroutine larnv_cplx

  !!* double complex precision version of larnv
  subroutine larnv_dblecplx(iDist,iSeed,x)
    integer, intent(in)       :: iDist
    integer, intent(inout)    :: iSeed(4)
    complex(rdp), intent(out) :: x(:)

    integer :: n

    ASSERT(iDist > 0)
    ASSERT(iDist < 4)
    ASSERT(all(iSeed(:) >= 0))
    ASSERT(all(iSeed(:) <= 4095))
    ASSERT(mod(iSeed(4),2) == 1)
    ASSERT(size(x) > 0)
    n = size(x)
    x(:) = 0.0d0
    call ZLARNV( iDist, iSeed, n, x )
  end subroutine larnv_dblecplx
  
end module lapackRoutines
