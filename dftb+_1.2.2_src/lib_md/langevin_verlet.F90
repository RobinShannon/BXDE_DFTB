!!* Velocity Verlet intergrator.
module VelocityVerlet
#include "assert.h"
#include "allocate.h"
use Accuracy
use FileId
use message
implicit none
private

public :: OLangevinVerlet
public :: create, destroy, bxdInvert, bxdInvert2, next, state

character(len=120) :: error_string !* Used to return runtime diagnostics

!!* Data for the integrator.
type OLangevinVerlet
private
integer :: nAtom                     !* Nr. of atoms
integer :: temperature
real(dp) :: deltaT                   !* time step for the integrator
real(dp) :: friction
real(dp), pointer :: positions(:,:) !* list of particle positions
real(dp), pointer :: velocities(:,:) !* list of particle velocities
real(dp), pointer :: oldPositions(:,:) !* list of particle positions
real(dp), pointer :: oldVelocities(:,:) !* list of particle velocities
logical           :: vHalfPresent = .false. !* do we have the v(t-.5)
end type OLangevinVerlet

interface create
module procedure LangevinVerlet_velocities
end interface

interface destroy
module procedure LangevinVerlet_destroy
end interface

interface bxdInvert
module procedure LangevinVerlet_bxdInvert
end interface

interface bxdInvert2
module procedure LangevinVerlet_bxdInvert2
end interface

interface next
module procedure LangevinVerlet_next
end interface

interface rescale
module procedure LangevinVerlet_rescale
end interface

interface state
module procedure LangevinVerlet_state
end interface

contains

!!* Creates a VelocityVerlet object from given external velocities for the
!!* t-th time step, this means later we have to reconstruct the Vel. Verlet
!!* t+.5 velocities
!!* @param self Pointer to the initialised object on exit.
!!* @param deltaT Integration time step.
!!* @param positions Position of the atoms.
!!* @param pThermostat Pointer to a thermostat.
!!* @param velocities list of initial velocities
subroutine LangevinVerlet_velocities(self, deltaT, positions, velocities, temperature, friction)
type(OLangevinVerlet), pointer       :: self
real(dp), intent(in)                 :: deltaT
real(dp), intent(in)                 :: positions(:,:)
real(dp), intent(in)                 :: velocities(:,:)

ASSERT(size(positions, dim=1) == 3)

INITALLOCATE_P(self)
self%nAtom = size(positions, dim=2)
INITALLOCATE_PARR(self%velocities, (3, self%nAtom))
INITALLOCATE_PARR(self%positions, (3, self%nAtom))
INITALLOCATE_PARR(self%OldVelocities, (3, self%nAtom))
INITALLOCATE_PARR(self%oldPositions, (3, self%nAtom))

self%deltaT = deltaT
self%temperature = temperature
self%friction = friction
self%positions(:,:) = positions(:,:)

self%velocities(:,:) = velocities(:,:)

self%vHalfPresent = .false. ! assumes the V read in corresponds to the
! current coordinates, so we should reconstruct the t+.5 velocities when
! possible once forces are available for the coordinates

end subroutine LangevinVerlet_velocities


!!* removes an integrator example
!!* @param self the instance to deallocate
subroutine LangevinVerlet_destroy(self)
type(OLangevinVerlet), pointer :: self

if (.not. associated(self)) then
return
end if
DEALLOCATE_PARR(self%velocities)
DEALLOCATE_PARR(self%positions)
DEALLOCATE_PARR(self%oldVelocities)
DEALLOCATE_PARR(self%oldPositions)
end if
DEALLOCATE_P(self)

end subroutine LangevinVerlet_destroy

!!* Perform velocity inversion according to BXD constraints
!!* Currently set up for energy constraints but can be readily modified
!!* @param self Pointer to the initialised object on exit.
!!* @param masses Mass matrix for atoms in the system
!!* @param delPhi Derivative of the constraint with respect to positions

subroutine LangevinVerlet_bxdInvert(self, masses, delPhi)
type(OVelocityVerlet), pointer       :: self
real(dp), intent(in)                 :: masses(:,:)
real(dp), intent(in)                 :: delPhi(:,:)


integer :: ii
integer :: jj
real(dp) :: lambda
real(dp) :: numerator
real(dp) :: denominator

real(dp), allocatable :: invM_transPhi(:,:)
real(dp), allocatable :: invM(:,:)
ALLOCATE_(invM_transPhi,(3,self%nAtom))
ALLOCATE_(invM,(3,self%nAtom))

self%positions(:,:) = self%oldPositions(:,:)
self%vHalfPresent = .false.
!! Invert masses compentent wise to get invM
invM(:,:) = 1.0_dp / masses(:,:)
!! Get equivalent of mass multipled by transpose of del_Phi
invM_transPhi(:,:) = delPhi(:,:) * invM(:,:)
DO ii = 1, 3
print *, 'ii'
print *, ii
numerator = numerator + dot_product((-2.0_dp * delPhi(ii,:)),self%oldVelocities(ii,:))
END DO
DO jj = 1, 3
denominator = denominator + dot_product(delPhi(jj,:),invM_transPhi(jj,:))
print *, 'denom'
print *, denominator
END DO
lambda = numerator / denominator

!! calculate new velocities and store
self%velocities(:,:) = self%oldVelocities(:,:) + (lambda * invM_transPhi(:,:))
self%oldVelocities(:,:) = self%velocities(:,:)

end subroutine LangevinVerlet_bxdInvert



subroutine LangevinVerlet_bxdInvert2(self, masses, delPhi, delPhi2)
type(OLangevinVerlet), pointer       :: self
real(dp), intent(in)                 :: masses(:,:)
real(dp), intent(in)                 :: delPhi(:,:)
real(dp), intent(in)                 :: delPhi2(:,:)


integer :: ii
real(dp) :: lambda1
real(dp) :: lambda2
real(dp) :: a
real(dp) :: b
real(dp) :: c
real(dp) :: d
real(dp) :: e
real(dp) :: f

real(dp), allocatable :: invM_transPhi(:,:)
real(dp), allocatable :: invM_transPhi2(:,:)
real(dp), allocatable :: invM(:,:)
ALLOCATE_(invM_transPhi,(3,self%nAtom))
ALLOCATE_(invM_transPhi2,(3,self%nAtom))
ALLOCATE_(invM,(3,self%nAtom))

self%positions(:,:) = self%oldPositions(:,:)
self%vHalfPresent = .false.
!! Invert masses compentent wise to get invM
invM(:,:) = 1.0_dp / masses(:,:)
!! Get equivalent of mass multipled by transpose of del_Phi
invM_transPhi(:,:) = delPhi(:,:) * invM(:,:)
invM_transPhi2(:,:) = delPhi2(:,:) * invM(:,:)

DO ii = 1, 3
a = a + dot_product((delPhi(ii,:)),invM_transPhi(ii,:))
END DO

DO ii = 1, 3
b = b + dot_product((delPhi(ii,:)),invM_transPhi2(ii,:))
END DO

DO ii = 1, 3
c = c + dot_product((delPhi2(ii,:)),invM_transPhi(ii,:))
END DO

DO ii = 1, 3
d = d + dot_product((delPhi2(ii,:)),invM_transPhi2(ii,:))
END DO

DO ii = 1, 3
e = e + dot_product((delPhi(ii,:)),self%oldVelocities(ii,:))
END DO

DO ii = 1, 3
f = f + dot_product((delPhi2(ii,:)),self%oldVelocities(ii,:))
END DO


lambda2 = (( c * d) - ( a * f ))/(( b * c ) - ( a * d ))
lambda1 = (c - lambda2 * b) / a

!! calculate new velocities and store
self%velocities(:,:) = self%oldVelocities(:,:) + (lambda1 * invM_transPhi(:,:)) + (lambda2 * invM_transPhi2(:,:))
self%oldVelocities(:,:) = self%velocities(:,:)

end subroutine LangevinVerlet_bxdInvert2




!!* Takes a timestep for the MD integrator, optionally with a thermostat.
!!* @param self integrator to propogate
!!* @param accel Accelerations.
!!* @param newCoord displaced coordinates
!!* @param newVelocity velocity of displaced coords
!!* @caveat Due to the way the velocity Verlet is split, the returned
!!* velocity is for 1 complete MD step behind the returned positions
!!* so print positions, then call next and then print velocities
!!* to get agreement between the positions and velocities.
subroutine LangevinVerlet_next(self, accel, newCoord, newVelocity, idx, count, reversal)
type(OLangevinVerlet), pointer :: self
real(dp),intent(in) :: accel(:,:)
real(dp),intent(out) :: newCoord(:,:)
real(dp),intent(out) :: newVelocity(:,:)
integer, intent(in) :: idx(:)
integer, intent(in) :: count
logical, intent(in) :: reversal

integer  :: ii

real(dp) :: n
real(dp) :: a
real(dp) :: b
real(dp) :: c

CALL RANDOM_NUMBER(n)

a = (2.0_dp - self%friction * self%deltaT) / (2.0_dp + self%friction * self%deltaT)
b = sqrt( self%temperature * self%deltT * 0.5_dp * self%friction )
c = 2 * self%deltaT / (2 + self%friction * self%deltaT)

self%oldPositions = self%positions(:,:)


ASSERT(associated(self))

newCoord(:,:) = 0.0_dp
newVelocity(:,:) = 0.0_dp

! start from the usual ordering of velocity verlet method (two cycles
! shown):
! a.1 v(t+.5dt) = v(t)        + .5*a(t)*dt -- a(t) external
! a.2 r(t + dt) = r(t)        + v(t+.5dt)*dt
! a.3 v(t + dt) = v(t+.5dt)   + .5*a(t+dt)*dt -- a(t+dt) external call
! b.1 v(t+1.5dt) = v(t+dt)    + .5*a(t+dt)*dt -- a(t+dt) external
! b.2 r(t + 2dt) = r(t+dt)    + v(t+1.5dt)*dt
! b.3 v(t + 2dt) = v(t+1.5dt) + .5*a(t+2dt)*dt -- a(t+2dt) external call
!
! and cut out a.3 b.1 b.2 as the cycle :
! a.3 v(t)      = v(t-.5dt)+ .5*a(t)*dt -- a(t) input
! b.1 v(t+.5dt) = v(t)     + .5*a(t)*dt
! b.2 r(t+dt)   = r(t)     + v(t+.5dt)*dt
!
! so :
! v_out   = v_store + .5*a_input*dt
! v_store = v_out   + .5*a_input*dt
! r_out   = r_store + v_store*dt
! r_store = r_out
! where v_out is one MD step behind the positions returned.

if (self%vHalfPresent) then
newVelocity(:,:) = self%velocities(:,:) &
& + 0.5_dp * accel(:,:) * self%deltaT * b * n
else
newVelocity(:,:) = self%velocities(:,:)
self%vHalfPresent=.true.
end if

self%oldVelocities(:,:) = newVelocity(:,:)
self%velocities(:,:) = a * newVelocity(:,:) + 0.5_dp * accel(:,:) * self%deltaT + b * n
newCoord(:,:) = self%positions(:,:) + self%velocities(:,:) * c
self%positions(:,:) = newCoord(:,:)

if (associated(self%pThermostat)) then
call updateVelocities(self%pThermostat, self%velocities)
end if

end subroutine VelocityVerlet_next

!!* Rescale the cell parameters and coordinates according to the tensorial
!!* version of the Berensen barostat (allows cell shape changes if the
!!* external pressure/stress is non isotropic)
!!* @param self integrator to rescale
!!* @param coord atom coordinates to rescale
!!* @param latVecs lattice vectors to rescale
!!* @param pressureTensor system stress tensor
!!* @note the forms of the isotropic and anisotropic Beresdsen barostats in
!!* the literature are slightly incompatible in their definitions
subroutine VelocityVerlet_rescale(self,coord,latVecs,pressureTensor)
type(OVelocityVerlet), pointer :: self
real(dp),intent(inout)         :: coord(:,:)
real(dp),intent(inout)         :: latVecs(3,3)
real(dp),intent(in)            :: pressureTensor(3,3)

real(dp) :: scale(3,3)
real(dp) :: scaleIso, Pext, P
integer  :: ii

ASSERT(self%tBarostat)

if (self%tIsotropic) then ! isotropic Berendsen, not quite consistent
! with anisotropic but its in the literature...
Pext = 0.0_dp
P = 0.0_dp
do ii = 1, 3
Pext = self%Pressure(ii,ii) / 3.0_dp
P = P + pressureTensor(ii,ii) / 3.0_dp
end do
scaleIso = (1.0_dp - self%BarostatStrength*(Pext - P))**(1.0_dp/3.0_dp)
self%positions(:,:) = self%positions(:,:) * scaleIso
coord(:,:) = coord(:,:) * scaleIso
latVecs(:,:) = latVecs(:,:) * scaleIso
else
scale = 0.0_dp
do ii = 1, 3
scale(ii,ii) = 1.0_dp
end do
scale = scale - self%BarostatStrength*(self%Pressure-pressureTensor)
do ii = 1, self%nAtom
self%positions(:,ii) = matmul(self%positions(:,ii),scale)
coord(:,ii) = matmul(coord(:,ii),scale)
end do
latVecs(:,:) = matmul(latVecs(:,:),scale)
end if

end subroutine VelocityVerlet_rescale

subroutine VelocityVerlet_state(self,fd)
type(OVelocityVerlet), pointer :: self
integer,intent(in)             :: fd

if (associated(self%pThermostat)) then
call state(self%pThermostat,fd)
end if

end subroutine VelocityVerlet_state

end module VelocityVerlet
