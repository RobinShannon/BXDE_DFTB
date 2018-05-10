!!* Where the compiler requires an external complementary error function to be
!!* called add in a wrapper interface for the needed functions, since both the
!!* real and the double precision cases then should have different function
!!* names. This external wrapper mechanism only becomes active if EXTERNALERFC
!!* is defined in the preprocessor.
module external_erfc
  implicit none
  private

  public :: extErfc, extErf
  
!!* wrapper for the two precisions of erfc
!!* @todo add the complex cases if needed - is the function defined for them?
  interface extErfc
     module procedure real_extErfc
     module procedure dble_extErfc
  end interface


  !!* wrapper for the two precisions of erfc
  !!* @todo add the complex cases if needed - is the function defined for them?
  interface extErf
     module procedure real_extErf
     module procedure dble_extErf
   end interface

   
contains

  !!* real precision external complementary error function routine for
  !!* compilers that get this from an external library
  !!* @param x value to calculate erfc(x) of
  function real_extErfc(x)
    use accuracy, only : rsp
    implicit none
    real(rsp) :: real_extErfc
    real(rsp) :: erfc
    real(rsp), intent(in) :: x
#ifdef EXTERNALERFC
    external :: erfc
#endif
    real_extErfc = erfc(x)
  end function real_extErfc


  
  !!* double precision external complementary error function routine for
  !!* compilers that get this from an external library
  !!* @param x value to calculate erfc(x) of
  function dble_extErfc(x)
    use accuracy, only : rdp
    implicit none
    real(rdp) :: dble_extErfc
    real(rdp), intent(in) :: x
    real(rdp) :: derfc
#ifdef EXTERNALERFC
    external :: derfc
#endif
    dble_extErfc = derfc(x)
  end function dble_extErfc


  !!* real precision external error function routine for
  !!* compilers that get this from an external library
  !!* @param x value to calculate erf(x) of
  function real_extErf(x)
    use accuracy, only : rsp
    implicit none
    real(rsp) :: real_extErf
    real(rsp) :: erf
    real(rsp), intent(in) :: x
#ifdef EXTERNALERFC
    external :: erf
#endif
    real_extErf = erf(x)
  end function real_extErf


  
  !!* double precision external error function routine for
  !!* compilers that get this from an external library
  !!* @param x value to calculate erf(x) of
  function dble_extErf(x)
    use accuracy, only : rdp
    implicit none
    real(rdp) :: dble_extErf
    real(rdp), intent(in) :: x
    real(rdp) :: derf
#ifdef EXTERNALERFC
    external :: derf
#endif
    dble_extErf = derf(x)
  end function dble_extErf

end module external_erfc
