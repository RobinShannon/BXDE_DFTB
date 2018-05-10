!! Note, this file is not read in any more. It had been kept since, it contains
!! a few draft ideas for the velocity Omelyan integrator.
!!* Velocity-Verlet and similar integrator routines
module Integrator
#include "allocate.h"
#include "assert.h"
  use accuracy
  use message
  use constants
  use ranlux
  implicit none

  private

  public :: OIntegrator
  public :: create, next, destroy
  public :: optimal_xi

  !!* Creates MD instance
  interface create
    module procedure Integrator_create
    module procedure Integrator_create_thermostat
  end interface

  !!* Destroys an MD instance
  interface destroy
    module procedure Integrator_destroy
  end interface
  
  !!* next MD step
  interface next
    module procedure Integrator_next
  end interface


  !!* Contains data for the integrator
  type OIntegrator
    private
    type(ORanlux),pointer :: pRanlux     !* random number generator
    integer               :: integrator  !* choice of integrator
    real(dp)              :: DeltaT      !* time step for the integrator
    real(dp), pointer     :: mass(:)     !* list of particle masses
    real(dp), pointer :: position(:,:,:) !* list of particle positions
                                         !* over time (last index if needed
                                         !* by integrator)
    real(dp)          :: kT              !* thermal energy
    integer           :: thermostat      !* choice of thermostat
    integer           :: step            !* step in multi-step integrators
    real(dp)          :: xi              !* Constant for Omelyan 2nd
                                         !* order integrator
    real(dp)          :: wvScale         !* Andersen thermostat probability
    real(dp), pointer :: velocity(:,:,:) !* list of particle velocities
                                         !* over time (last index if needed
                                         !* by integrator)
  end type OIntegrator

  !!* choices of integrator - velocity Verlet
  integer, parameter :: velocityVerlet = 1

  !!* choices of integrator - Omelyan, Mryglod and Folk 2nd order
  integer, parameter :: velocityOmelyan = 2


contains

  !!* Initialiser for a velocity-Verlet/Omelyan algorithm, with the velocity
  !!* at t = .5 step chosen either as a Maxwell-Boltzmann distribution or 
  !!* from supplied velocities
  !!* @param self The integrator to initialise
  !!* @param pRanlux Ranlux random generator
  !!* @param inChoice Which integrator to initialise as
  !!* @param DeltaT Time step
  !!* @param masses list of particle masses for the chemical species
  !!* @param species list of atomic species for all atom
  !!* @param positions list of initial positions
  !!* @param kT initial thermal energy, in Hartree, to choose a Boltzmann
  !!* distribution if velocities are not suplied
  !!* @param velocities initial velocities if supplied
  !!* @ref Omelyan, Mryglod and Folk, Phys. Rev. E, 65, 056706 (2002).
  !!* @author Ben Hourahine
  !!* @note either velocities or an initial temperature must be set.
  subroutine Integrator_create(self,pRanlux,intChoice,DeltaT,masses,species, &
      & positions, kT, velocities)
    type(OIntegrator), pointer :: self
    type(ORanlux), pointer     :: pRanlux
    integer,  intent(in)       :: intChoice
    real(dp), intent(in)       :: DeltaT
    real(dp), intent(in)       :: masses(:)
    integer, intent(in)        :: species(:)
    real(dp), intent(in)       :: positions(:,:)
    real(dp), intent(in), optional :: kT
    real(dp), intent(in), optional :: velocities(:,:)

    integer :: ii

    ASSERT(.not.associated(self))
    ASSERT(size(positions,dim=1) == 3)
    ASSERT(associated(pRanlux))
    ASSERT(present(kT).neqv.present(velocities))

    if (intChoice /= velocityOmelyan .and. intChoice /= velocityVerlet) then
      call error("Unknown MD integrator")
    end if
    if (intChoice == velocityOmelyan) then
      call error("Omelyan, Mryglod and Folk not working yet, sorry!")
    end if

    INITALLOCATE_P(self)
    INITALLOCATE_PARR(self%mass,(size(positions,dim=2)))

    if (self%integrator == velocityVerlet) then
      INITALLOCATE_PARR(self%velocity,(3,size(positions,dim=2),1))
      INITALLOCATE_PARR(self%position,(3,size(positions,dim=2),1))
    else
      INITALLOCATE_PARR(self%velocity,(3,size(positions,dim=2),2))
      INITALLOCATE_PARR(self%position,(3,size(positions,dim=2),2))
    end if

    self%pRanlux => pRanlux
    self%integrator = intChoice
    self%DeltaT = DeltaT
    do ii = 1, size(positions,dim=2)
      self%mass(ii) = masses(species(ii))
    end do
    self%position(:,:,1) = positions(:,:)
    self%thermostat = 0
    self%step = 0
    self%xi = 0.0_dp
    ! set the xi parameter if the Omelyan-type integrator
    if (intChoice == velocityOmelyan) then
      self%xi = optimal_xi()
    end if

    if (present(velocities)) then
      self%velocity(:,:,1) = velocities(:,:)
    else
      self%kT = kT
      call Integrator_init_velocity(self)
    end if
    
  end subroutine Integrator_create

  

  !!* Initialiser for a thermostated velocity-Verlet/Omelyan algorithm, 
  !!* with the velocity at t = .5 step chosen either as a Maxwell-Boltzmann 
  !!* distribution or from supplied velocities
  !!* @param self The integrator to initialise
  !!* @param pRanlux Ranlux random generator
  !!* @param inChoice Which integrator to initialise as
  !!* @param DeltaT Time step
  !!* @param masses list of particle masses for the chemical species
  !!* @param species list of atomic species for all atom
  !!* @param positions list of initial positions
  !!* @param kT initial thermal energy, in Hartree, to choose a Boltzmann
  !!* distribution if velocities are not suplied
  !!* @param thermostat choice of thermostating method
  !!* @parakm wvScale rescale probability
  !!* @param velocities initial velocities if supplied
  !!* @ref Omelyan, Mryglod and Folk, Phys. Rev. E, 65, 056706 (2002).
  !!* @author Ben Hourahine
  !!* @todo other thermostats, other integrators than v. Verlet and constraints
  subroutine Integrator_create_thermostat(self,pRanlux,intChoice,DeltaT,masses,&
      & species, positions, kT,thermostat,wvScale,velocities)
    type(OIntegrator), pointer :: self
    type(ORanlux), pointer     :: pRanlux
    integer,  intent(in)       :: intChoice
    real(dp), intent(in)       :: DeltaT
    real(dp), intent(in)       :: masses(:)
    integer, intent(in)        :: species(:)
    real(dp), intent(in)       :: positions(:,:)
    real(dp), intent(in)       :: kT
    integer, intent(in)        :: thermostat
    real(dp), intent(in)       :: wvScale
    real(dp), intent(in), optional :: velocities(:,:)

    integer :: ii

    ASSERT(.not.associated(self))
    ASSERT(size(positions,dim=1) == 3)
    ASSERT(associated(pRanlux))

    if (intChoice /= velocityOmelyan .and. intChoice /= velocityVerlet) then
      call error("Unknown MD integrator")
    end if
    if (intChoice == velocityOmelyan) then
      call error("Omelyan, Mryglod and Folk not working yet, sorry!")
    end if

    INITALLOCATE_P(self)
    INITALLOCATE_PARR(self%mass,(size(positions,dim=2)))

    if (self%integrator == velocityVerlet) then
      INITALLOCATE_PARR(self%velocity,(3,size(positions,dim=2),1))
      INITALLOCATE_PARR(self%position,(3,size(positions,dim=2),1))
    else
      INITALLOCATE_PARR(self%velocity,(3,size(positions,dim=2),2))
      INITALLOCATE_PARR(self%position,(3,size(positions,dim=2),2))
    end if

    self%pRanlux => pRanlux
    self%integrator = intChoice
    self%DeltaT = DeltaT
    do ii = 1, size(positions,dim=2)
      self%mass(ii) = masses(species(ii))
    end do
    self%position(:,:,1) = positions(:,:)
    self%kT = kT
    if (thermostat == 1 .or. thermostat == 2) then
      self%thermostat = thermostat
    else
      call error("Unknown thermostat!")
    end if
    self%wvScale = wvScale
    self%step = 0
    self%xi = 0.0_dp
    ! set the xi parameter if the Omelyan-type integrator
    if (intChoice == velocityOmelyan) then
      self%xi = optimal_xi()
    end if

    if (present(velocities)) then
      self%velocity(:,:,1) = velocities(:,:)
    else
      call Integrator_init_velocity(self)
    end if
    
  end subroutine Integrator_create_thermostat


  
  !!* Takes a timestep for the MD integrator, optionally with a thermostat.
  !!* Two versions of the Andersen thermostat are implemented, either the 
  !!* @param self integrator propogate
  !!* @param force forces acting on the integrator
  !!* @param newCoord displaced coordinates
  !!* @param newVelocity velocity of displaced coords
  !!* @param KE kinetic energy of particles
  !!* @param T temperature in K of the system
  !!* @caveat Due to the way the velocity Verlet is split, the returned
  !!* velocity is for 1 complete MD step behind the returned positions 
  !!* so print positions, then call next and then print velocities
  !!* to get agreement between the positions and velocities.
  !!* @todo other integrators than v. Verlet and constraints, other thermostats
  !!* than Andersen
  !!* @ref Andersen J. Chem. Phys. 72. 2384 (1980)
  subroutine Integrator_next(self,force,newCoord,newVelocity,KE,T)
    type(OIntegrator), pointer :: self
    real(dp),intent(in) :: force(:,:)
    real(dp),intent(out) :: newCoord(:,:)
    real(dp),intent(out) :: newVelocity(:,:)
    real(dp),intent(out) :: KE
    real(dp),intent(out), optional :: T
    
    real(dp) :: rescaleChance
    integer :: ii

    ASSERT(associated(self))

    newCoord(:,:) = 0.0_dp
    newVelocity(:,:) = 0.0_dp
    
    if (self%integrator == velocityVerlet) then
      
      if (self%thermostat == 1) then ! Andersen thermostat with a 
        ! probability that each atom's velocity can be drawn from a Maxwell
        ! Boltzmann distribution
        do ii = 1, size(self%mass(:))
          call getRandom(self%pRanlux, rescaleChance)
          if (rescaleChance <= self%wvScale) then
            call MaxwellBoltzmann( self%velocity(:,ii,1), &
                & self%mass(ii), self%kT, self%pRanlux )
          end if
        end do
        call rest_frame(self%velocity(:,:,1),self%mass)
      end if

      if (self%thermostat == 2) then ! Andersen thermostat again, but
        ! all atoms re-set at random
        call getRandom(self%pRanlux, rescaleChance)
        if (rescaleChance <= self%wvScale) then
          call Integrator_init_velocity(self)
        end if
      end if

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

      newVelocity(:,:) = self%velocity(:,:,1) + &
          & 0.5_dp *  * &
          & self%DeltaT
      self%velocity(:,:,1) = newVelocity(:,:) + &
          & 0.5_dp * force(:,:)/spread(self%mass(:),1,3) * &
          & self%DeltaT
      newCoord(:,:) = self%position(:,:,1) + &
          & self%velocity(:,:,1) * self%DeltaT
      self%position(:,:,1) = newCoord(:,:)
    else
      call error("This integrator not implemented yet, sorry")
    end if
    
    call eval_KE(KE,newVelocity(:,:),self%mass)
    if (present(T)) then
      call eval_T(T,newVelocity(:,:),self%mass,(3*size(self%mass)-3))
    end if
    
    self%step = self%step + 1

  end subroutine Integrator_next


  
  !!* removes an integrator example
  !!* @param self the instance to deallocate
  subroutine Integrator_destroy(self)
    type(OIntegrator), pointer :: self

    if (associated(self)) then
      DEALLOCATE_PARR(self%mass)
      DEALLOCATE_PARR(self%velocity)
      DEALLOCATE_PARR(self%position)
    end if
    DEALLOCATE_P(self)

  end subroutine Integrator_destroy

  
  !!* Set initial velocities from a random Maxwell-Boltzmann distribution
  !!* @parameter self an MD integrator instance
  !!* @param pRanlux Ranlux random generator
  !!* @todo use atom resolved energies/forces to pre-condition the velocities
  !!* towards a canonical distribution, i.e. attempt to distribute the total
  !!* energy at random correctly instead of just the kinetic part.
  subroutine Integrator_init_velocity(self)
    type(OIntegrator), pointer :: self
    
    integer  :: ii, jj
    real(dp) :: KE
    real(dp) :: ranvals(7)
    
    character(len=11), parameter :: formatVelocity = '(i4,3f12.5)'
    
    ASSERT(associated(self%pRanlux))
    
    do ii = 1, size(self%mass(:))
      call MaxwellBoltzmann( self%velocity(:,ii,1),self%mass(ii), &
          & self%kT, self%pRanlux )
    end do

    call rest_frame(self%velocity(:,:,1),self%mass)
    call eval_KE(KE,self%velocity(:,:,1),self%mass)
    !write(*,*)'Initial KE',KE
    call rescaleToT(self%velocity(:,:,1),self%mass,self%kT,KE)
    !write(*,*)'Target KE per degree of Freedom',self%kT
    call eval_KE(KE,self%velocity(:,:,1),self%mass)
    !write(*,*)'Scaled KE',KE

    
#if DEBUG >= 2 
    write(*,*)'Initial velocities for MD'
    do ii=1,size(self%velocity,dim=2)
      write(*,formatVelocity)ii,self%velocity(:,ii,1)
    end do
    write(*,*)
#endif
    
  end subroutine Integrator_init_velocity



  !!* return the optimal splitting value for the optimized Verlet-like
  !!* integrator of Omelyan, Mryglod and Folk
  function optimal_xi()
    real(dp) :: optimal_xi

    optimal_xi = 0.5_dp &
        & - ((2.0_dp*sqrt(326.0_dp)+36.0_dp)**(1.0_dp/3.0_dp))/12.0_dp &
        & + 1.0_dp/(6.0_dp*((2.0_dp*sqrt(326.0_dp)+36.0_dp)**(1.0_dp/3.0_dp)))

  end function optimal_xi


  
  !!* Shift velocities so the average momentum is 0
  !!* @param velocity particle velocities
  !!* @param mass particle masses
  subroutine  rest_frame(velocity,mass)
    real(dp), intent(inout) :: velocity(:,:)
    real(dp), intent(in) :: mass(:)
    
    real(dp) :: mv(3)
    
    ! calculate total momentum of the system
    mv(:) = sum( spread(mass(:),1,3) * velocity(:,:), dim=2)

    ! per atom
    mv(:) = mv(:)/real(size(mass),dp)
    
    ! and shift so that it is 0
    velocity(:,:) = velocity(:,:) - &
        & spread( mv(:), 2, size(mass) ) / spread( mass, 1, 3)

  end subroutine rest_frame


  
  !!* Calculate the kinetic energy of an integrator
  !!* @parameter KE resulting energy
  !!* @param velocity particle velocities
  !!* @param mass particle masses
  subroutine eval_KE(KE,velocity,mass)
    real(dp),intent(out) :: KE
    real(dp), intent(in) :: velocity(:,:)
    real(dp), intent(in) :: mass(:)
    
    ASSERT(size(velocity,dim=2)==size(mass))
    ASSERT(size(velocity,dim=1)==3)

    KE = 0.5_dp * sum(spread(mass(:),1,3) * velocity(:,:)**2 )
    
  end subroutine eval_KE

  
  
  !!* Calculate the kinetic temperature of an integrator
  !!* @parameter T resulting temperature in Kelvin
  !!* @param velocity particle velocities
  !!* @param mass particle masses
  !!* @param Nf degrees of freedom
  subroutine eval_T(T,velocity,mass,Nf)
    real(dp),intent(out) :: T
    real(dp), intent(in) :: velocity(:,:)
    real(dp), intent(in) :: mass(:)
    integer, intent(in)  :: Nf
    
    ASSERT(size(velocity,dim=2)==size(mass))
    ASSERT(size(velocity,dim=1)==3)
    
    T = sum(spread(mass(:),1,3) * velocity(:,:)**2 ) / (real(Nf,dp) * Boltzmann)
    
  end subroutine eval_T

  

  !!* Rescales the velocities of a system to match the target thermal energy
  !!* @param velocity particle velocities
  !!* @param mass particle masses
  !!* @param kTtarget intended kinetic energy per degree of freedom
  !!* @param KE current kinetic energy of the particles
  subroutine rescaleToT(velocity,mass,kTtarget,KE)
    real(dp), intent(inout) :: velocity(:,:)
    real(dp), intent(in)    :: mass(:)
    real(dp), intent(in)    :: kTtarget
    real(dp), intent(in)    :: KE
    
    ASSERT(size(velocity,dim=2)==size(mass))
    ASSERT(size(velocity,dim=1)==3)
    
    ! should the scale be 3*natoms-3 as the translational degrees of
    ! freedom are 0 in this frame?
    velocity(:,:) = velocity(:,:) * sqrt(kTtarget*real(3*size(mass)-3,dp) / &
        & (2.0_dp * KE))
    
  end subroutine rescaleToT


  
  !!* Converts a uniform distribution into a Gaussian distribution
  !!* @parameter eta1 number with Gaussian distribution
  !!* @parameter eta2 number with Gaussian distribution
  !!* @parameter u1 number with uniform distribution
  !!* @parameter u2 number from uniform distribution
  subroutine BoxMueller(eta1,eta2,u1,u2)
    real(dp), intent(out) :: eta1
    real(dp), intent(out) :: eta2
    real(dp), intent(in)  :: u1
    real(dp), intent(in)  :: u2
    
    real(dp) :: theta, a

    a = sqrt( -2.0_dp*log(u1) )
    theta = 2.0_dp*pi*u2

    eta1 = a * cos(theta)
    eta2 = a * sin(theta)

  end subroutine BoxMueller


  
  !!* draws an atom velocity from a Maxwell-Boltzmann distribution
  !!* @param velocity resulting velocity
  !!* @param mass atomic mass in a.u.
  !!* @param kT thermal energy in H for each degree of freedom
  !!* @param pRanlux pointer to a random number generator
  subroutine MaxwellBoltzmann(velocity,mass,kT,pRanlux)
    real(dp), intent(out)  :: velocity(:)
    real(dp), intent(in)   :: mass
    real(dp), intent(in)   :: kT
    type(ORanlux), pointer :: pRanlux
    
    real(dp) :: ranvals(7)
    real(dp) :: junk
    
    ASSERT(size(velocity)==3)
    
    ! use the uniform distribution to get a normal (Gaussian) distribution
    ! and then scale to Maxwell-Boltzmann 
    
    call getRandom(pRanlux, ranvals)
    
    call BoxMueller(velocity(1),velocity(2),ranvals(1),ranvals(2))
    call BoxMueller(velocity(3),junk,ranvals(3),ranvals(4))
    
    if (ranvals(5) < 0.5_dp) then
      velocity(1) = - velocity(1)
    end if
    if (ranvals(6) < 0.5_dp) then
      velocity(2) = - velocity(2)
    end if
    if (ranvals(7) < 0.5_dp) then
      velocity(3) = - velocity(3)
    end if
    
    velocity(:) = velocity(:) * sqrt(kT/mass)    
   
  end subroutine MaxwellBoltzmann

end module Integrator
