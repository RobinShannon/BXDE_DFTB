!!* contains some miscellaneous QM related bits and pieces.
module QM
#include "assert.h"
  use accuracy, only : dp
  
  implicit none
  
  public

  !!* constructs a commutator  
  interface commutator
    module procedure C_
  end interface

  ! perform a unitary transformation on a matrix $X^\prime = U X U^\dag$
  interface unitary
    module procedure U_
  end interface
  
contains
  
  !!* constructs a commutator for given matrices C = [A,B]
  !!* @param C result of commutator
  !!* @param first matrix
  !!* @param second matrix
  subroutine C_(C,A,B)
    complex(dp), intent(out) :: C(:,:)
    complex(dp), intent(in)  :: A(:,:)
    complex(dp), intent(in)  :: B(:,:)

    ASSERT(all(shape(C)==shape(A)))
    ASSERT(all(shape(C)==shape(B)))
    ASSERT(size(C,dim=1)==size(C,dim=2))
    
    C = matmul(A,B) - matmul(B,A)
    
  end subroutine C_

  !!* unitary transformation on a matrix $X^\prime = U X U^\dag$
  !!* @param Xprime matrix in new basis
  !!* @param X matrix in original basis
  !!* @param U unitary matrix
  !!* @todo test that U is actually unitary in an assert_env block
  subroutine U_(Xprime,X,U)
    complex(dp), intent(out) :: Xprime(:,:)
    complex(dp), intent(in) :: X(:,:)
    complex(dp), intent(in) :: U(:,:)
    
    complex(dp) :: work(size(X,dim=1),size(X,dim=2))

    ASSERT(all(shape(Xprime)==shape(X)))
    ASSERT(all(shape(X)==shape(U)))
    ASSERT(size(X,dim=1)==size(X,dim=2))
    
    work = matmul(x,transpose(conjg(u)))
    Xprime = matmul(u,work)
    
  end subroutine U_
  
end module QM
