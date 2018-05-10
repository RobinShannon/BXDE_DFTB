# 1 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
!!* Global variables and initialization for the main program
!!* @todo Assignment (copy) operator for TNeighbors!!!
module initprogram

# 1 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/../includes/assert.h" 1
!! -*- f90 -*-
!! vim:syntax=fortran:
!!
!! Provides C-style assertions for Fortran
!! (compile the source with 0>=1 to use it)
!!









!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Debug turned on
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# 35 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/../includes/assert.h"









!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# 5 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90" 2

# 1 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/../includes/allocate.h" 1
!! -*- f90 -*-
!! Provides trapped allocate and deallocate for Fortran
!! (compile the source with 0>=2 to use the verbose version)
!!

use allocation




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Low-level macros.
!! Do NOT use them, use the 0-level specific ones instead
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




























































!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Debug turned on and its level is high enough
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# 115 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/../includes/allocate.h"



























!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Macros for initialization and allocation in one step
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!









!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# 6 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90" 2
  use inputdata_
  use constants
  use periodic
  use accuracy
  use intrinsicprinter
  use short_gamma
  use coulomb
  use message

  use Mixer
  use SimpleMixer
  use AndersonMixer
  use BroydenMixer
  use DIISMixer

  use GeoOpt
  use ConjugateGradient
  use SteepestDescent

  use Ranlux
  use MDCommon
  use MDintegrator
  use Connectivity
  use VelocityVerlet
  use Thermostat
  use DummyThermostat
  use AndersenThermostat
  use BerendsenThermostat
  use NHCThermostat
  use TempProfile
  use numDerivs2
  use LapackRoutines
  use SimpleAlgebra

  use SCC
  use SCC_init
  use SlakoCont
  use RepCont

  use Spin, only: Spin_getOrbitalEquiv, ud2qm, qm2ud
  use DFTBplsU

  use Dispersion
  use DispSlaterKirkwood
  use DispUFF

  use ThirdOrder

  use OrbitalEquiv
  use CommonTypes
  use sorting, only : heap_sort
  use Fifo
  use LinkedList
  implicit none


  character(*), parameter :: fChargeIn = "charges.bin"
  character(*), parameter :: fStopDriver = "stop_driver"
  character(*), parameter :: fStopSCC = "stop_scc"

  logical               :: tSCC            !* Is the calculation SCC?
  integer,  parameter   :: nCutoff = 1     !* Nr. of different cutoffs

  integer               :: nAtom           !* nr. of atoms
  integer               :: nAllAtom        !* nr. of all (image and orig) atoms
  integer,  pointer     :: Img2CentCell(:) !* nr. of original atom in centre
  integer               :: nType           !* nr of different types (nAtom)

  type(TOrbitals), pointer :: orb

  integer               :: nOrb            !* nr. of orbitals in the system
  integer               :: nAllOrb         !* nr. of orbitals for all atoms
  integer,  pointer     :: specie(:)       !* types of the atoms (nAllAtom)
  integer,  allocatable :: specie0(:)      !* type of the atoms (nAtom)
  real(dp), pointer     :: coord(:,:)      !* Coords of the atoms (3, nAllAtom)
  real(dp), allocatable, target :: coord0(:,:)   !* Coords (3, nAtom)
  real(dp), allocatable :: tmpCoords(:)    !* temporary array of coords
  real(dp), allocatable :: tmpWeight(:)    !* temporary weights
  real(dp), allocatable :: tmp3Coords(:,:) !* temporary array of coords (3,:)
  logical               :: tPeriodic       !* if calculation is periodic
  logical               :: tShowFoldedCoord!* Should central cell coordinates
  !* be output?
  logical               :: tFracCoord      !* are atomic coordinates fractional?
  real(dp)              :: sccTol          !* Tollerance for SCC cycle

  real(dp), allocatable, target :: latVec(:,:)  !* lattice vectors as columns
  real(dp), allocatable, target :: recVec(:,:)  !* reciprocal vecs as columns
  real(dp)              :: origLatVec(3,3) !* original lattice vectors if
  !* optimizing, normalized

  real(dp), allocatable :: recVec2p(:,:)   !* reciprocal vectors in 2pi units
  real(dp)              :: volume          !* cell volume
  real(dp)              :: recVolume       !* reciprocal cell volume
  real(dp),  pointer    :: cellVec(:,:)    !* translation vecs for interacting
                                           !* image cells (3, nImgCell + 1)
  real(dp), pointer     :: rCellVec(:,:)   !* cell vectors in absolute units
  integer,  pointer     :: iCellVec(:)     !* index in cellVec for each atom

  !!* ADT for neighbor parameters
  type(TNeighborList), save :: neighborList
  integer,  allocatable :: nNeighbor(:)    !* nr. of neighbors for SK + rep
  integer,  pointer     :: iPair(:,:)      !* H/S indexing array

  integer,  allocatable :: iAtomStart(:)   !* atom start pos for squared H/S

  real(dp), allocatable :: hubbU(:,:)      !* Hubbard Us (orbital, atom)
  real(dp), allocatable :: atomEigVal(:,:) !* self energy (orbital, atom)
  real(dp), allocatable :: referenceN0(:,:)!* reference n_0 charges for each
                                           !* atom
  real(dp), allocatable :: mass(:)         !* list of atomic masses

  type(OSlakoCont), pointer :: skHamCont
  type(OSlakoCont), pointer :: skOverCont
  type(ORepCont), pointer :: pRepCont
  real(dp) :: skCutoff
  real(dp) :: skRepCutoff
  logical :: tUseBuggyRepSum   !* if old (buggy) rep. summation should be used

  real(dp)              :: mCutoff        !* longest pair interaction

  real(dp), pointer     :: ham(:,:)       !* Hamiltonian
  real(dp), pointer     :: iHam(:,:)      !* imaginary part of the Hamiltonian
  real(dp), allocatable :: chargePerShell(:,:,:)
  real(dp), pointer     :: over(:)        !* Overlap


  integer               :: nKPoint        !* nr. of K-points
  real(dp), allocatable :: kPoint(:,:)    !* K-points
  real(dp), allocatable :: KWeight(:)     !* weight of the K-Points

  real(dp)              :: pressure        !* external pressure if periodic
  logical               :: tBarostat
  real(dp)              :: BarostatStrength

  logical               :: tRealHS        !* H and S are real

  real(dp), allocatable :: nEl(:)         !* nr. of electrons
  real(dp)              :: nEl0           !* Nr. of all electrons if neutral
  real(dp), allocatable :: W(:,:,:)       !* Spin W's	  !'
  real(dp), allocatable :: xi(:,:)        !* Spin orbit constants

  logical               :: tDFTBU         !* is this a DFTB+U calculation?
  integer               :: nDFTBUfunc     !* Choice of orbital functional
  real(dp), allocatable :: UJ(:,:)        !* list of U-J for species
  integer, allocatable  :: nUJ(:)         !* How many U-J for each species
  integer, allocatable  :: niUJ(:,:)      !* number of l-values of U-J for each
  !* block
  integer, allocatable  :: iUJ(:,:,:)     !* l-values of U-J for each block

  real(dp) :: tempElec                    !* electron temperature
  logical :: tFillKSep                    !* If K points should filled separat.
  logical :: tSetFillingTemp              !* Filling temp updated by MD.
  integer  :: iDistribFn = 0              !* Choice of electron distribution
                                          !* function, defaults to Fermi
  real(dp) :: tempAtom                    !* atomic temperature
  real(dp) :: deltaT                      !* MD stepsize
  integer :: solver                       !* eigensolver
  integer :: nSCCIter                     !* number of SCC iterations
  integer :: nSpin                        !* Number of spin components, 1
                                          !* is unpolarised, 2 is polarised, 4
  !* is noncolinear / spin-orbit
  logical :: tSpin                        !* is this a spin polarized
  !* calculation?
  logical :: tSpinOrbit                   !* is there spin-orbit coupling
  logical :: tDualSpinOrbit               !* Use block like dual representation
  logical :: tImHam                       !* complex hamiltonian in real space
  logical :: t2Component                  !* is this a two component
  !* calculation (spin orbit or non-collinear spin)

  real(dp) :: almix

  logical               :: tGeoOpt          !* Geometry optimization needed?
  logical               :: tCoordOpt        !* optimize internal coordinates?
  logical               :: tLatOpt          !* optimize lattice constants?
  logical               :: tLatOptFixAng    !* Fix angles between lattice
  !!* vectors when optimizing?
  logical               :: tLatOptFixLen(3)
  logical               :: tLatOptIsotropic
  logical               :: tMD              !* Is this a MD calculation?
  logical               :: tDerivs          !* Is this a derivatives calc?
  logical               :: tMulliken        !* Do we need Mulliken charges?
  logical               :: tDipole          ! calculate an electric dipole?
  logical               :: tAtomicEnergy    !* Do we need atom resolved E?
  logical               :: tEigenvecs       !* Print out eigenvectors?
  logical               :: tProjEigenvecs   !* Print eigenvector projections?
  logical               :: tForces          !* Do we need forces?
  integer               :: nMovedAtom       !* Number of moved atoms
  integer, allocatable  :: indMovedAtom(:)  !* Index of the moved atoms
  integer               :: nMovedCoord      !* Nr. of moved coordinates
  integer               :: nGeoSteps        !* Nr. of geo movements to do
  integer               :: nGeoConstr       !* Nr. of geometry constraints
  integer,  allocatable :: conAtom(:)       !* Index of constrained atoms
  real(dp), allocatable :: conVec(:,:)      !* Constraint vectors

  character(lc) :: geoOutFile     !* File containing end geometry

  logical               :: tAppendGeo       !* Append geometries in the output?
  logical :: tConvrgForces
  character(mc), allocatable :: specieName(:)

  real(dp)              :: random_pool(10)   !* pool of initial random numbers
  ! for future use. See comment in code at create(pRanlux, in this routine.

  logical            :: tInitialized = .false.
  private :: tInitialized

  type(OGeoOpt), pointer :: pGeoCoordOpt  !* General geometry optimizer
  type(OGeoOpt), pointer :: pGeoLatOpt    !* Geometry optimizer for lattice
                                          !* consts

  type(OMixer), pointer :: pChrgMixer    !* Charge mixer
  integer               :: iMixer        !* mixer number

  type(ORanlux), pointer :: pRanlux !* Random number generator

  type(OMDCommon), pointer :: pMDFrame  !* MD Framework
  type(OMDIntegrator), pointer :: pMDIntegrator !* MD integrator
  type(OConnect), pointer :: pConnect
  type(OTempProfile), pointer :: pTempProfile

  type(OnumDerivs), pointer :: pDerivDriver

  !! BXD input
  logical                  :: BXD
  logical                  :: COMbxd
  integer                  :: bxdBoxes
  integer                  :: bxdEquil
  integer                  :: bxdEvents
  real(dp)                 :: bxdLower
  real(dp)                 :: bxdUpper
  logical                  :: DOS
  integer                  :: DOSiter

  !! Connectivity
  logical                  :: revAll
  logical                  :: revDis
  integer                  :: rxnWindow
  logical                  :: endOnReact
  logical                  :: restrain
  integer                  :: COMindex
  real(dp)                 :: minCOM
  logical                  :: bxdWait

  !! Charge related variables
  real(dp), allocatable    :: q0(:, :, :)
  real(dp), allocatable    :: qShell0(:,:)
  real(dp), allocatable    :: qInput(:, :, :)
  real(dp), allocatable    :: qOutput(:, :, :)
  real(dp), allocatable    :: qBlockIn(:, :, :, :)
  real(dp), allocatable    :: qBlockOut(:, :, :, :)
  real(dp), allocatable    :: qiBlockIn(:, :, :, :)
  real(dp), allocatable    :: qiBlockOut(:, :, :, :)
  real(dp), allocatable    :: qInpRed(:), qOutRed(:), qDiffRed(:)
  integer, allocatable     :: iEqOrbitals(:,:,:)  !* Orbital equiv. relations
  integer :: nIneqOrb !* nr. of inequivalent orbitals
  integer :: nMixElements !* nr. of elements to go through the mixer - may
  ! contain reduced orbitals and also orbital blocks (if tDFTBU)
  integer, allocatable :: iEqBlockDFTBU(:,:,:,:) !* Orbital equivalency for
  !* orbital blocks
  integer, allocatable :: iEqBlockDFTBULS(:,:,:,:) !* Orbital equivalency for
  !* orbital blocks with spin-orbit


  !! External charges
  integer :: nExtChrg   !* Nr. of external charges
  logical :: tExtChrg   !* If external charges must be considered

  logical  :: tEField = .false. ! external electric field
  real(dp) :: EFieldStrength = 0.0_dp ! field strength
  real(dp) :: EfieldVector(3) = 0.0_dp ! field direction
  logical  :: tTDEfield = .false. ! time dependent
  real(dp) :: EfieldOmega = 0.0_dp ! angular frequency
  integer  :: EfieldPhase = 0 ! phase of field at step 0

  ! PDOS projection
  type(listIntR1), save :: iOrbRegion
  type(listCharLc), save :: regionLabels

  !! Third order
  logical :: t3rdFull
  type(OThirdOrder) :: thirdOrd

  !! Other stuff
  logical :: tReadChrg    !* If initial charges/dens mtx. from external file.
  logical :: tWriteTagged !* produce tagged output?
  logical :: tWriteDetailedXML !* Produce detailed.xml
  logical :: tWriteResultsTag !* Produce detailed.tag
  logical :: tWriteDetailedOut !* Produce detailed.out
  logical :: tWriteBandDat !* Produce band.dat
  logical :: tWriteHS, tWriteRealHS  !* Should HS (square and real) be printed?
  logical :: tMinMemory, tStoreEigvecs


  integer :: runId !* Program run id

  integer :: restartFreq    !* Frequency for saving restart info

  logical :: tDispersion  !* If dispersion should be calculated
  logical :: tStress = .true. !* Can stress be calculated? - start by
  !* assuming it can
  type(ODispersion), save :: myDispersion

  type(OFifoRealR2), allocatable :: storeEigvecsReal(:)
  type(OFifoCplxR2), allocatable :: storeEigvecsCplx(:)


  integer, parameter :: nInitNeighbor = 40  !* First guess for nr. of neighbors.
  private :: nInitNeighbor



contains


  !!* Initializes the variables in the module based on the parsed input
  !!* @param input Holds the parsed input data.
  subroutine initProgramVariables(input)
    type(inputData), intent(inout), target :: input

    !! Mixer related local variables
    integer  :: nGeneration
    real(dp) :: mixParam
    type(OSimpleMixer), pointer :: pSimpleMixer
    type(OAndersonMixer), pointer :: pAndersonMixer
    type(OBroydenMixer), pointer :: pBroydenMixer
    type(ODIISMixer), pointer :: pDIISMixer

    !! Geometry optimizer related local variables
    type(OConjGrad), pointer :: pConjGrad    !* Conjugate gradient driver
    type(OSteepDesc), pointer :: pSteepDesc  !* Steepest descent driver
    type(OConjGrad), pointer :: pConjGradLat   !* Conjugate gradient driver
    type(OSteepDesc), pointer :: pSteepDescLat !* Steepest descent driver

    !! MD related local variables
    type(OThermostat), pointer :: pThermostat
    type(ODummyThermostat), pointer :: pDummyTherm
    type(OAndersenThermostat), pointer :: pAndersenTherm
    type(OBerendsenThermostat), pointer :: pBerendsenTherm
    type(ONHCThermostat), pointer :: pNHCTherm

    type(OVelocityVerlet), pointer :: pVelocityVerlet

    integer :: ind, ii, jj, kk, iS, iAt, iSp, iSh
    integer :: iStart, iEnd

    !! Dispersion
    type(ODispSlaKirk), pointer :: pSlaKirk
    type(ODispUFF), pointer :: pVdWUFF

    character(lc) :: strTmp, strTmp2
    logical :: tFirst
    integer  :: iSeed, timeValues(8)
    real(dp) :: rTmp

    logical :: tExist ! Flag if some files do exist or not

    !! Orbital equivalency for SCC and Spin
    integer, allocatable :: iEqOrbSCC(:,:,:), iEqOrbSpin(:,:,:)
    !! Orbital equivalency for orbital potentials
    integer, allocatable :: iEqOrbDFTBU(:,:,:)

    !! Damped interactions
    logical, allocatable, target :: tDampedShort(:)
    type(TThirdOrderInp) :: thirdInp

    !! PDOS stuff
    integer :: iReg, nAtomRegion, nOrbRegion, iTmp
    integer, allocatable :: iAtomRegion(:)
    integer :: valshape(1)
    character(lc) :: tmpStr
    integer, allocatable :: tmpir1(:)

    type(TSCCInit), save :: sccInit



    !! Basic variables
    tSCC = input%ctrl%tScc
    tDFTBU = input%ctrl%tDFTBU
    tSpin = input%ctrl%tSpin
    if (tSpin) then
      nSpin = 2
    else
      nSpin = 1
    end if
    tSpinOrbit = input%ctrl%tSpinOrbit
    tDualSpinOrbit = input%ctrl%tDualSpinOrbit
    t2Component = input%ctrl%t2Component

    if (t2Component) then
      nSpin = 4
    end if

   !! BXD input
    BXD = input%ctrl%BXD
    COMbxd = input%ctrl%comBXD
    bxdBoxes = input%ctrl%bxdBoxes
    bxdEquil = input%ctrl%bxdEquil
    bxdEvents = input%ctrl%bxdEvents
    bxdLower = input%ctrl%bxdLower
    bxdUpper = input%ctrl%bxdUpper
    DOS = input%ctrl%DOS
    DOSiter = input%ctrl%DOSiter

  !! Connectivity
    revAll = input%ctrl%revAll
    revDis = input%ctrl%revDis
    rxnWindow = input%ctrl%rxnWindow
    endOnReact = input%ctrl%endOnReact
    restrain = input%ctrl%restrain
    COMindex = input%ctrl%COMindex
    minCOM = input%ctrl%minCOM
    bxdWait = input%ctrl%bxdWait

    sccTol = input%ctrl%sccTol
    nAtom = input%geom%nAtom
    nType = input%geom%nSpecie
    tPeriodic = input%geom%tPeriodic
    tShowFoldedCoord = input%ctrl%tShowFoldedCoord
    if (tShowFoldedCoord .and. .not. tPeriodic) then
      call error("Folding coordinates back into the central cell is&
          & meaningless for molecular boundary conditions!")
    end if
    tFracCoord = input%geom%tFracCoord
    solver = input%ctrl%iSolver
    if (tSCC) then
      nSCCIter = input%ctrl%maxIter
    else
      nSCCIter = 1
    end if

    if (tPeriodic) then
if (.not. allocated(latVec)) then;    allocate(latVec   (3, 3), stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas&
&01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 436, ioerr, 1);  else;    call allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 436&
&, "Array already allocated", 1);  endif
# 437 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"

      latVec(:,:) = input%geom%latVecs(:,:)
if (.not. allocated(recVec)) then;    allocate(recVec   (3, 3), stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas&
&01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 439, ioerr, 1);  else;    call allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 439&
&, "Array already allocated", 1);  endif
# 440 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (.not. allocated(recVec2p)) then;    allocate(recVec2p   (3, 3), stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/pan&
&asas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 440, ioerr, 1);  else;    call allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 440&
&, "Array already allocated", 1);  endif
# 441 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
      recVec2p = latVec(:,:)
      call matinv(recVec2p)
      recVec2p = reshape(recVec2p, (/3, 3/), order=(/2, 1/))
      recVec = 2.0_dp * pi * recVec2p
      volume = determinant33(latVec)
      recVolume = determinant33(recVec)
    else
if (.not. allocated(latVec)) then;    allocate(latVec   (0, 0), stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas&
&01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 448, ioerr, 1);  else;    call allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 448&
&, "Array already allocated", 1);  endif
# 449 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (.not. allocated(recVec)) then;    allocate(recVec   (0, 0), stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas&
&01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 449, ioerr, 1);  else;    call allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 449&
&, "Array already allocated", 1);  endif
# 450 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (.not. allocated(recVec2p)) then;    allocate(recVec2p   (0, 0), stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/pan&
&asas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 450, ioerr, 1);  else;    call allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 450&
&, "Array already allocated", 1);  endif
# 451 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
      volume = 0.0_dp
      recVolume = 0.0_dp
    end if

    orb => input%slako%orb

    !! Slater-Koster tables
    skHamCont => input%slako%skHamCont
    skOverCont => input%slako%skOverCont
    pRepCont => input%slako%repCont
    tUseBuggyRepSum = input%ctrl%useBuggyRepSum

if (.not. allocated(atomEigVal)) then;    allocate(atomEigVal   (orb%mShell, nType), stat=ioerr);    if (ioerr /= 0) call allocate&
&Error("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 463, ioerr, 1);  else;    call allocateError&
&("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 463, "Array already allocated", 1);  endif
# 464 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"


    atomEigVal(:,:) = input%slako%skSelf(1:orb%mShell, :)


if (.not. allocated(referenceN0)) then;    allocate(referenceN0  (orb%mShell, nType), stat=ioerr);    if (ioerr /= 0) call allocat&
&eError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 469, ioerr, 1);  else;    call allocateErro&
&r("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 469, "Array already allocated", 1);  endif
# 470 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
    referenceN0(:,:) = input%slako%skOcc(1:orb%mShell, :)

if (.not. allocated(mass)) then;    allocate(mass  (nType), stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas01/c&
&hem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 472, ioerr, 1);  else;    call allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 472&
&, "Array already allocated", 1);  endif
# 473 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
    mass(:) = input%slako%mass(:)
    ! Spin W's	!'
    if (tSpin) then
if (.not. allocated(W)) then;    allocate(W  (orb%mShell,orb%mShell,nType), stat=ioerr);    if (ioerr /= 0) call allocateError("/p&
&anfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 476, ioerr, 1);  else;    call allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 476&
&, "Array already allocated", 1);  endif
# 477 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
      W(:,:,:) = 0.0_dp
      do ii=1,nType
        do jj=1, orb%nShell(ii)
          do kk=1, orb%nShell(ii)
            W(jj,kk,ii)=input%ctrl%W(jj,kk,ii)
          end do
        end do
      end do
    else
if (.not. allocated(W)) then;    allocate(W  (0,0,0), stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas01/chem/rs&
&15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 486, ioerr, 1);  else;    call allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 486&
&, "Array already allocated", 1);  endif
# 487 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
    end if

    if (tSpinOrbit) then
if (.not. allocated(xi)) then;    allocate(xi  (orb%mShell,nType), stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/pana&
&sas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 490, ioerr, 1);  else;    call allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 490&
&, "Array already allocated", 1);  endif
# 491 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
      xi(:,:) = 0.0_dp
      do ii=1,nType
        do jj=1, orb%nShell(ii)
          xi(jj,ii)=input%ctrl%xi(jj,ii)
        end do
      end do
    else
if (.not. allocated(xi)) then;    allocate(xi  (0,0), stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas01/chem/rs&
&15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 498, ioerr, 1);  else;    call allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 498&
&, "Array already allocated", 1);  endif
# 499 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
    end if

    ! DFTB+U parameters
    if (tDFTBU) then
      nDFTBUfunc = input%ctrl%DFTBUfunc
if (.not. allocated(UJ)) then;    allocate(UJ  (size(input%ctrl%UJ,dim=1),size(input%ctrl%UJ,dim=2)), stat=ioerr);    if (ioerr /=&
& 0) call allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 504, ioerr, 1);  else;    c&
&all allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 504, "Array already allocated", &
&1);  endif
# 505 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (.not. allocated(nUJ)) then;    allocate(nUJ  (size(input%ctrl%nUJ)), stat=ioerr);    if (ioerr /= 0) call allocateError("/panf&
&s/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 505, ioerr, 1);  else;    call allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 505&
&, "Array already allocated", 1);  endif
# 506 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (.not. allocated(niUJ)) then;    allocate(niUJ  (size(input%ctrl%niUJ,dim=1),size(input%ctrl%niUJ,dim=2)), stat=ioerr);    if (&
&ioerr /= 0) call allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 506, ioerr, 1);  el&
&se;    call allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 506, "Array already allo&
&cated"&
&, 1);  endif
# 507 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (.not. allocated(iUJ)) then;    allocate(iUJ   (size(input%ctrl%iUJ,dim=1),         size(input%ctrl%iUJ,dim=2),size(input%ctrl%&
&iUJ,dim=3)), stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/init&
&program.F90"&
&, 508, ioerr, 1);  else;    call allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 508&
&, "Array already allocated", 1);  endif
# 508 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"


      UJ(:,:) = input%ctrl%UJ(:,:)
      nUJ(:) = input%ctrl%nUJ(:)
      niUJ(:,:) = input%ctrl%niUJ(:,:)
      iUJ(:,:,:) = input%ctrl%iUJ(:,:,:)
      do ii = 1, nType
        do jj = 1, nUJ(ii)
          if (niUJ(jj,ii)>1) then
            call heap_sort(iUJ(1:niUJ(jj,ii),jj,ii))
          end if
        end do
      end do
    else
if (.not. allocated(UJ)) then;    allocate(UJ  (0,0), stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas01/chem/rs&
&15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 522, ioerr, 1);  else;    call allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 522&
&, "Array already allocated", 1);  endif
# 523 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (.not. allocated(nUJ)) then;    allocate(nUJ  (0), stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas01/chem/rs&
&15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 523, ioerr, 1);  else;    call allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 523&
&, "Array already allocated", 1);  endif
# 524 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (.not. allocated(niUJ)) then;    allocate(niUJ  (0,0), stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas01/che&
&m/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 524, ioerr, 1);  else;    call allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 524&
&, "Array already allocated", 1);  endif
# 525 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (.not. allocated(iUJ)) then;    allocate(iUJ  (0,0,0), stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas01/che&
&m/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 525, ioerr, 1);  else;    call allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 525&
&, "Array already allocated", 1);  endif
# 526 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
    end if

    !! Cutoffs from SlaKo and repulsive
    skCutoff = max(getCutoff(skHamCont), getCutoff(skOverCont))
    skRepCutoff = max(skCutoff, getCutoff(pRepCont))
    mCutoff = skRepCutoff

    !! Get specie names and output file
    geoOutFile = input%ctrl%outFile
if (.not. allocated(specieName)) then;    allocate(specieName   (size(input%geom%specieNames)), stat=ioerr);    if (ioerr /= 0) ca&
&ll allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 535, ioerr, 1);  else;    call al&
&locateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 535, "Array already allocated", 1);  e&
&ndif
# 536 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
    specieName(:) = input%geom%specieNames(:)

    do ii = 1, nType
      do jj = ii+1, nType
        if (specieName(ii) == specieName(jj)) then
          write(tmpStr,"('Duplicate identical species labels in the geometry: ',&
              &A)")specieName(ii)
              call error(tmpStr)
        end if
      end do
    end do


    !! Initialise the SCC module (the two copies of the Hubbard Us are rather
    !! artifical, since the copy for the main program is only used for dumping
    !! into the tagged format for autotest)
if (.not. allocated(hubbU)) then;    allocate(hubbU   (orb%mShell, nType), stat=ioerr);    if (ioerr /= 0) call allocateError("/pa&
&nfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 552, ioerr, 1);  else;    call allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 552&
&, "Array already allocated", 1);  endif
# 553 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"


    hubbU(:,:) = input%slako%skHubbU(1:orb%mShell, :)
    if (tSCC) then
      sccInit%orb => orb
      if (tPeriodic) then
        sccInit%latVecs => latVec
        sccInit%recVecs => recVec
        sccInit%volume = volume
      end if
      sccInit%hubbU => input%slako%skHubbU
if (.not. allocated(tDampedShort)) then;    allocate(tDampedShort   (nType), stat=ioerr);    if (ioerr /= 0) call allocateError("/&
&panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 564, ioerr, 1);  else;    call allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 564&
&, "Array already allocated", 1);  endif
# 565 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
      if (input%ctrl%tDampH) then
        tDampedShort = (mass < 3.5_dp * amu__au)
        !tDampedShort(:) = (specieName == "H" .or. specieName == "h")
      else
        tDampedShort(:) = .false.
      end if
      sccInit%tDampedShort => tDampedShort
      sccInit%dampExp = input%ctrl%dampExp
      nExtChrg = input%ctrl%nExtChrg
      tExtChrg = (nExtChrg > 0)
      if (tExtChrg) then
        tStress = .false. ! Stress calculations not allowed


        sccInit%extCharges => input%ctrl%extChrg
        sccInit%blurWidths => input%ctrl%extChrgblurWidth
      end if
      if (associated(input%ctrl%chrgConstr)) then

        if (any(abs(input%ctrl%chrgConstr(:,2)) > epsilon(1.0_dp))) then
          sccInit%chrgConstraints => input%ctrl%chrgConstr
        end if
      end if

      if (associated(input%ctrl%thirdOrderOn)) then

        sccInit%thirdOrderOn => input%ctrl%thirdOrderOn
      end if

      sccInit%ewaldAlpha = input%ctrl%ewaldAlpha
      call init_SCC(sccInit)
      mCutoff = max(mCutoff, getSCCCutoff())

      ! Initialize full 3rd order module
      t3rdFull = input%ctrl%t3rdFull
      if (t3rdFull) then
        thirdInp%nAtom = nAtom
if (.not. allocated(thirdInp%hubbU)) then;    allocate(thirdInp%hubbU   (nType), stat=ioerr);    if (ioerr /= 0) call allocateErro&
&r("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 602, ioerr, 1);  else;    call allocateError("/p&
&anfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 602, "Array already allocated", 1);  endif
# 603 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
        thirdInp%hubbU(:) = hubbU(1,:)
if (.not. allocated(thirdInp%hubbUDeriv)) then;    allocate(thirdInp%hubbUDeriv   (nType), stat=ioerr);    if (ioerr /= 0) call al&
&locateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 604, ioerr, 1);  else;    call allocat&
&eError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 604, "Array already allocated", 1);  endif
# 605 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"

        thirdInp%hubbUDeriv(:) = input%ctrl%hubDerivs(:)
if (.not. allocated(thirdInp%damped)) then;    allocate(thirdInp%damped   (nType), stat=ioerr);    if (ioerr /= 0) call allocateEr&
&ror("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 607, ioerr, 1);  else;    call allocateError("&
&/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 607, "Array already allocated", 1);  endif
# 608 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
        thirdInp%damped = tDampedShort
        thirdInp%dampExp = input%ctrl%dampExp
        call init(thirdOrd, thirdInp)
        mCutoff = max(mCutoff, getCutoff(thirdord))
      end if
    end if


    !! Initial coordinates
if (.not. allocated(coord0)) then;    allocate(coord0   (3, nAtom), stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/pan&
&asas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 617, ioerr, 1);  else;    call allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 617&
&, "Array already allocated", 1);  endif
# 618 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"

    coord0(:,:) = input%geom%coords(:,:)
if (.not. allocated(specie0)) then;    allocate(specie0   (nAtom), stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/pana&
&sas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 620, ioerr, 1);  else;    call allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 620&
&, "Array already allocated", 1);  endif
# 621 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"

    specie0(:) = input%geom%species(:)
    if (tPeriodic) then
      !! Make some guess for the nr. of all interacting atoms
      nAllAtom = int((real(nAtom, dp)**(1.0_dp/3.0_dp) + 3.0_dp)**3)
    else
      nAllAtom = nAtom
    end if
nullify(coord);  if (.not. associated(coord)) then;    allocate(coord    (3, nAllAtom), stat=ioerr);    if (ioerr /= 0) call alloc&
&ateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 629, ioerr, 5);  else;    call allocateEr&
&ror("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 629, "Ptr. array already associated!",5);  end&
&if
# 630 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
nullify(specie);  if (.not. associated(specie)) then;    allocate(specie    (nAllAtom), stat=ioerr);    if (ioerr /= 0) call alloc&
&ateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 630, ioerr, 5);  else;    call allocateEr&
&ror("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 630, "Ptr. array already associated!",5);  end&
&if
# 631 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
nullify(img2CentCell);  if (.not. associated(img2CentCell)) then;    allocate(img2CentCell    (nAllAtom), stat=ioerr);    if (ioer&
&r /= 0) call allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 631, ioerr, 5);  else; 
   call allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 631, "Ptr. array already ass&
&ociated!"&
&,5);  endif
# 632 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
nullify(iCellVec);  if (.not. associated(iCellVec)) then;    allocate(iCellVec    (nAllAtom), stat=ioerr);    if (ioerr /= 0) call&
& allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 632, ioerr, 5);  else;    call allo&
&cateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 632, "Ptr. array already associated!",5)&
&;  endif
# 633 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (.not. allocated(iAtomStart)) then;    allocate(iAtomStart   (nAtom + 1), stat=ioerr);    if (ioerr /= 0) call allocateError("/&
&panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 633, ioerr, 1);  else;    call allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 633&
&, "Array already allocated", 1);  endif
# 634 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
    call buildSquaredAtomIndex(iAtomStart, orb)


    !! Initalize images (translations)
    if (tPeriodic) then
      nullify(cellVec)
      nullify(rCellVec)
      call getCellTranslations(cellVec, rCellVec, latVec, recVec2p, mCutoff)
    else
nullify(cellVec);  if (.not. associated(cellVec)) then;    allocate(cellVec    (3, 1), stat=ioerr);    if (ioerr /= 0) call alloca&
&teError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 643, ioerr, 5);  else;    call allocateErr&
&or("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 643, "Ptr. array already associated!",5);  endi&
&f
# 644 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
nullify(rCellVec);  if (.not. associated(rCellVec)) then;    allocate(rCellVec    (3, 1), stat=ioerr);    if (ioerr /= 0) call all&
&ocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 644, ioerr, 5);  else;    call allocate&
&Error("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 644, "Ptr. array already associated!",5);  e&
&ndif
# 645 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
      cellVec(:, 1) = (/ 0.0_dp, 0.0_dp, 0.0_dp /)
      rCellVec(:, 1) = (/ 0.0_dp, 0.0_dp, 0.0_dp /)
    end if

    !! Initialize neighborlist.
    call init(neighborList, nAtom, nInitNeighbor)
if (.not. allocated(nNeighbor)) then;    allocate(nNeighbor   (nAtom), stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/&
&panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 651, ioerr, 1);  else;    call allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 651&
&, "Array already allocated", 1);  endif
# 652 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"

    !! Intialize Hamilton and overlap
    if (tSCC) then
if (.not. allocated(chargePerShell)) then;    allocate(chargePerShell  (orb%mShell,nAtom,nSpin), stat=ioerr);    if (ioerr /= 0) c&
&all allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 655, ioerr, 1);  else;    call a&
&llocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 655, "Array already allocated", 1);  
endif
# 656 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
    else
if (.not. allocated(chargePerShell)) then;    allocate(chargePerShell  (0,0,0), stat=ioerr);    if (ioerr /= 0) call allocateError&
&("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 657, ioerr, 1);  else;    call allocateError("/pa&
&nfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 657, "Array already allocated", 1);  endif
# 658 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
    end if
nullify(ham);  if (.not. associated(ham)) then;    allocate(ham    (0, nSpin), stat=ioerr);    if (ioerr /= 0) call allocateError(&
&"/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 659, ioerr, 5);  else;    call allocateError("/pan&
&fs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 659, "Ptr. array already associated!",5);  endif
# 660 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
nullify(iHam);  if (.not. associated(iHam)) then;    allocate(iHam    (0, nSpin), stat=ioerr);    if (ioerr /= 0) call allocateErr&
&or("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 660, ioerr, 5);  else;    call allocateError("/&
&panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 660, "Ptr. array already associated!",5);  endif
# 661 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
nullify(over);  if (.not. associated(over)) then;    allocate(over    (0), stat=ioerr);    if (ioerr /= 0) call allocateError("/pa&
&nfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 661, ioerr, 5);  else;    call allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 661&
&, "Ptr. array already associated!",5);  endif
# 662 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
nullify(iPair);  if (.not. associated(iPair)) then;    allocate(iPair    (0, nAtom), stat=ioerr);    if (ioerr /= 0) call allocate&
&Error("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 662, ioerr, 5);  else;    call allocateError&
&("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 662, "Ptr. array already associated!",5);  endif
# 663 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"

    !! Brillouin zone sampling
    if (tPeriodic) then
      nKPoint = input%ctrl%nKPoint
if (.not. allocated(kPoint)) then;    allocate(kPoint   (3, nKPoint), stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/p&
&anasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 667, ioerr, 1);  else;    call allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 667&
&, "Array already allocated", 1);  endif
# 668 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (.not. allocated(kWeight)) then;    allocate(kWeight   (nKPoint), stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/pa&
&nasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 668, ioerr, 1);  else;    call allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 668&
&, "Array already allocated", 1);  endif
# 669 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"


      kPoint(:,:) = input%ctrl%KPoint(:,:)
      if (sum(input%ctrl%kWeight(:)) < epsilon(1.0_dp)) then
        call error("Sum of k-point weights should be greater than zero!")
      end if
      kWeight(:) = input%ctrl%kWeight(:) / sum(input%ctrl%kWeight(:))
    else
      nKPoint = 1
if (.not. allocated(kPoint)) then;    allocate(kPoint   (3, nKPoint), stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/p&
&anasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 678, ioerr, 1);  else;    call allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 678&
&, "Array already allocated", 1);  endif
# 679 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (.not. allocated(kWeight)) then;    allocate(kWeight   (nKpoint), stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/pa&
&nasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 679, ioerr, 1);  else;    call allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 679&
&, "Array already allocated", 1);  endif
# 680 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
      kPoint(:,1) = 0.0_dp
      kWeight(1) = 1.0_dp
    end if

    if ((.not. tPeriodic) .or. (nKPoint == 1 .and. &
        &all(kPoint(:, 1) == (/ 0.0_dp, 0.0_dp, 0.0_dp /)))) then
      tRealHS = .true.
    else
      tRealHS = .false.
    end if

    !! Other usefull quantities
    nOrb = orb%nOrb

    if (nSpin == 4) then
if (.not. allocated(nEl)) then;    allocate(nEl  (1), stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas01/chem/rs&
&15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 695, ioerr, 1);  else;    call allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 695&
&, "Array already allocated", 1);  endif
# 696 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
    else
if (.not. allocated(nEl)) then;    allocate(nEl  (nSpin), stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas01/che&
&m/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 697, ioerr, 1);  else;    call allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 697&
&, "Array already allocated", 1);  endif
# 698 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
    end if

    nEl0 = 0.0_dp
    do ii = 1, nAtom
      nEl0 = nEl0 + sum(input%slako%skOcc(1:orb%nShell(specie0(ii)), &
          & specie0(ii)))
    end do
    nEl(:) = 0.0_dp
    if (nSpin == 1 .or. nSpin == 4) then
      nEl(1) = nEl0 - input%ctrl%nrChrg
      if(ceiling(nEl(1)) > 2.0_dp*nOrb) then
        call error("More electrons than basis functions!")
      end if
    else
      nEl(1) = 0.5_dp * (nEl0 - input%ctrl%nrChrg + input%ctrl%nrSpinPol)
      nEl(2) = 0.5_dp * (nEl0 - input%ctrl%nrChrg - input%ctrl%nrSpinPol)
      if (any(ceiling(nEl(:)) > nOrb)) then
        call error("More electrons than basis functions!")
      end if
    end if

    if (.not.all(nEl(:) >= 0.0_dp)) then
      call error("Less than 0 electrons!")
    end if

    iDistribFn = input%ctrl%iDistribFn
    tempElec = input%ctrl%tempElec
    tSetFillingTemp = input%ctrl%tSetFillingTemp
    tFillKSep = input%ctrl%tFillKSep
    tempAtom = input%ctrl%tempAtom
    deltaT = input%ctrl%deltaT



    tImHam = tDualSpinOrbit ! .or. tBField

    !! Create equivalency relations
    if (tSCC) then
if (.not. allocated(iEqOrbitals)) then;    allocate(iEqOrbitals   (orb%mOrb, nAtom, nSpin), stat=ioerr);    if (ioerr /= 0) call a&
&llocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 736, ioerr, 1);  else;    call alloca&
&teError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 736, "Array already allocated", 1);  endif
# 737 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (.not. allocated(iEqOrbSCC)) then;    allocate(iEqOrbSCC   (orb%mOrb, nAtom, nSpin), stat=ioerr);    if (ioerr /= 0) call alloc&
&ateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 737, ioerr, 1);  else;    call allocateEr&
&ror("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 737, "Array already allocated", 1);  endif
# 738 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
      call SCC_getOrbitalEquiv(orb, specie0, iEqOrbSCC)
      if (nSpin == 1) then
        iEqOrbitals(:,:,:) = iEqOrbSCC(:,:,:)
      else
if (.not. allocated(iEqOrbSpin)) then;    allocate(iEqOrbSpin   (orb%mOrb, nAtom, nSpin), stat=ioerr);    if (ioerr /= 0) call all&
&ocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 742, ioerr, 1);  else;    call allocate&
&Error("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 742, "Array already allocated", 1);  endif
# 743 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
        call Spin_getOrbitalEquiv(orb, specie0, iEqOrbSpin)
        call OrbitalEquiv_merge(iEqOrbSCC, iEqOrbSpin, orb, iEqOrbitals)
if (allocated(iEqOrbSpin)) then;    deallocate(iEqOrbSpin, stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas01/ch&
&em/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&,745,ioerr, 2);  endif
# 746 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
      end if
if (allocated(iEqOrbSCC)) then;    deallocate(iEqOrbSCC, stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas01/chem&
&/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&,747,ioerr, 2);  endif
# 748 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
      nIneqOrb = maxval(iEqOrbitals)
      nMixElements = nIneqOrb
      if (tDFTBU) then
if (.not. allocated(iEqOrbSpin)) then;    allocate(iEqOrbSpin   (orb%mOrb, nAtom, nSpin), stat=ioerr);    if (ioerr /= 0) call all&
&ocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 751, ioerr, 1);  else;    call allocate&
&Error("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 751, "Array already allocated", 1);  endif
# 752 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (.not. allocated(iEqOrbDFTBU)) then;    allocate(iEqOrbDFTBU   (orb%mOrb, nAtom, nSpin), stat=ioerr);    if (ioerr /= 0) call a&
&llocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 752, ioerr, 1);  else;    call alloca&
&teError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 752, "Array already allocated", 1);  endif
# 753 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
        call DFTBplsU_getOrbitalEquiv(iEqOrbDFTBU,orb, specie0, nUJ, niUJ, iUJ)
        call OrbitalEquiv_merge(iEqOrbitals, iEqOrbDFTBU, orb, iEqOrbSpin)
        iEqOrbitals(:,:,:) = iEqOrbSpin(:,:,:)
        nIneqOrb = maxval(iEqOrbitals)
if (allocated(iEqOrbSpin)) then;    deallocate(iEqOrbSpin, stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas01/ch&
&em/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&,757,ioerr, 2);  endif
# 758 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (allocated(iEqOrbDFTBU)) then;    deallocate(iEqOrbDFTBU, stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas01/&
&chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&,758,ioerr, 2);  endif
# 759 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (.not. allocated(iEqBlockDFTBU)) then;    allocate(iEqBlockDFTBU  (orb%mOrb, orb%mOrb, nAtom, nSpin), stat=ioerr);    if (ioerr&
& /= 0) call allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 759, ioerr, 1);  else;  
  call allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 759, "Array already allocated&
&"&
&, 1);  endif
# 760 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
        call DFTBU_blockIndx(iEqBlockDFTBU, nIneqOrb, orb, specie0, &
            & nUJ, niUJ, iUJ)
        nMixElements = max(nMixElements,maxval(iEqBlockDFTBU)) ! as
        !  iEqBlockDFTBU does not include diagonal elements, so in the case of
        !  a purely s-block DFTB+U calculation, maxval(iEqBlockDFTBU) would
        !  return 0
        if (tImHam) then
if (.not. allocated(iEqBlockDFTBULS)) then;    allocate(iEqBlockDFTBULS  (orb%mOrb, orb%mOrb, nAtom, nSpin), stat=ioerr);    if (i&
&oerr /= 0) call allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 767, ioerr, 1);  els&
&e;    call allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 767, "Array already alloc&
&ated"&
&, 1);  endif
# 768 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
          call DFTBU_blockIndx(iEqBlockDFTBULS,nMixElements , orb, specie0, &
            & nUJ, niUJ, iUJ)
          nMixElements = max(nMixElements,maxval(iEqBlockDFTBULS))
        end if
      end if
    else
      nIneqOrb = nOrb
      nMixElements = 0
    end if

    if (.not.tDFTBU) then
if (.not. allocated(iEqBlockDFTBU)) then;    allocate(iEqBlockDFTBU  (0, 0, 0, 0), stat=ioerr);    if (ioerr /= 0) call allocateEr&
&ror("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 779, ioerr, 1);  else;    call allocateError("&
&/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 779, "Array already allocated", 1);  endif
# 780 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
    end if
    if (.not.(tDFTBU.and.tImHam)) then
if (.not. allocated(iEqBlockDFTBULS)) then;    allocate(iEqBlockDFTBULS  (0, 0, 0, 0), stat=ioerr);    if (ioerr /= 0) call alloca&
&teError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 782, ioerr, 1);  else;    call allocateErr&
&or("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 782, "Array already allocated", 1);  endif
# 783 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
    end if


    !! Initialize mixer
    !! (at the moment, the mixer doesn't need to know about the size of the
    !! vector to mix.)
    if (tSCC) then
      iMixer = input%ctrl%iMixSwitch
      nGeneration = input%ctrl%iGenerations
      mixParam = input%ctrl%almix
      select case (iMixer)
      case (1)
        call create(pSimpleMixer, mixParam)
        call create(pChrgMixer, pSimpleMixer)
      case (2)
        if (input%ctrl%andersonNrDynMix > 0) then
          call create(pAndersonMixer, nGeneration, mixParam, &
              &input%ctrl%andersonInitMixing, input%ctrl%andersonDynMixParams, &
              &input%ctrl%andersonOmega0)
        else
          call create(pAndersonMixer, nGeneration, mixParam, &
              &input%ctrl%andersonInitMixing, omega0=input%ctrl%andersonOmega0)
        end if
        call create(pChrgMixer, pAndersonMixer)
      case (3)
        call create(pBroydenMixer, nSCCIter, mixParam, &
            &input%ctrl%broydenOmega0, input%ctrl%broydenMinWeight, &
            &input%ctrl%broydenMaxWeight, input%ctrl%broydenWeightFac, &
            &nGeneration)
        call create(pChrgMixer, pBroydenMixer)
      case(4)
        call create(pDIISMixer,nGeneration, mixParam, input%ctrl%tFromStart)
        call create(pChrgMixer, pDIISMixer)
      case default
        call error("Unknown charge mixer type.")
      end select
    else
      nullify(pChrgMixer)
    end if

    !! initialise in cases where atoms move
    tGeoOpt = input%ctrl%tGeoOpt
    tCoordOpt = input%ctrl%tCoordOpt
    tLatOpt = (input%ctrl%tLatOpt .and. tPeriodic)
    if (tLatOpt) then
      if (tExtChrg) then
        ! Stop as not sure, what to do with the coordinates of the
        ! external charges, when the lattice changes.
        call error("External charges and lattice optimisation can not be used &
            &together.")
      end if
    end if
    if (tLatOpt) then
      tLatOptFixAng = input%ctrl%tLatOptFixAng
      tLatOptFixLen = input%ctrl%tLatOptFixLen
      tLatOptIsotropic = input%ctrl%tLatOptIsotropic
      if (tLatOptFixAng .or. any(tLatOptFixLen) .or. tLatOptIsotropic) then
        origLatVec(:,:) = latVec(:,:)
      end if
    end if
    pressure = input%ctrl%pressure
    tBarostat = input%ctrl%tBarostat
    BarostatStrength = input%ctrl%BarostatStrength

    tAppendGeo = input%ctrl%tAppendGeo
    tConvrgForces = (input%ctrl%tConvrgForces .and. tSCC) ! no point if not SCC
    tMD = input%ctrl%tMD
    tDerivs = input%ctrl%tDerivs
    tMulliken = input%ctrl%tMulliken
    tAtomicEnergy = input%ctrl%tAtomicEnergy
    tEigenvecs = input%ctrl%tEigenvecs

    ! Projection of eigenstates onto specific regions of the system
    tProjEigenvecs = input%ctrl%tProjEigenvecs
    if (tProjEigenvecs) then
      call init(iOrbRegion)
      call init(regionLabels)
      do iReg = 1, size(input%ctrl%tShellResInRegion)
        call elemShape(input%ctrl%iAtInRegion, valshape, iReg)
        nAtomRegion = valshape(1)
if (.not. allocated(iAtomRegion)) then;    allocate(iAtomRegion   (nAtomRegion), stat=ioerr);    if (ioerr /= 0) call allocateErro&
&r("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 863, ioerr, 1);  else;    call allocateError("/p&
&anfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 863, "Array already allocated", 1);  endif
# 864 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
        call intoArray(input%ctrl%iAtInRegion, iAtomRegion, iTmp, iReg)
        if (input%ctrl%tShellResInRegion(iReg)) then
          iSp = specie0(iAtomRegion(1))

          ! Create a separate region for each shell. It will contain
          ! the orbitals of that given shell for each atom in the region.
          do iSh = 1, orb%nShell(iSp)
            nOrbRegion = nAtomRegion &
                &* (orb%posShell(iSh + 1, iSp) - orb%posShell(iSh, iSp))
            ind = 1
            ! Create orbital index.
if (.not. allocated(tmpir1)) then;    allocate(tmpir1   (nOrbRegion), stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/p&
&anasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 875, ioerr, 1);  else;    call allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 875&
&, "Array already allocated", 1);  endif
# 876 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
            do ii = 1, nAtomRegion
              iAt = iAtomRegion(ii)
              do jj = orb%posShell(iSh, iSp), orb%posShell(iSh + 1, iSp) - 1
                tmpir1(ind) = iAtomStart(iAt) + jj - 1
                ind = ind + 1
              end do
            end do
            call append(iOrbRegion, tmpir1)
if (allocated(tmpir1)) then;    deallocate(tmpir1, stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas01/chem/rs151&
&64/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&,884,ioerr, 2);  endif
# 885 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
            write(tmpStr, "(A,'.',I0,'.out')") &
                &trim(input%ctrl%RegionLabel(iReg)), iSh
            call append(regionLabels, tmpStr)
          end do
        else
          ! We take all orbitals from all atoms.
          nOrbRegion = 0
          do ii = 1, nAtomRegion
            nOrbRegion = nOrbRegion + orb%nOrbAtom(iAtomRegion(ii))
          end do
          ind = 1
if (.not. allocated(tmpir1)) then;    allocate(tmpir1   (nOrbRegion), stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/p&
&anasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 896, ioerr, 1);  else;    call allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 896&
&, "Array already allocated", 1);  endif
# 897 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
          ! Create an index of the orbitals
          do ii = 1, nAtomRegion
            iAt = iAtomRegion(ii)
            do jj = 1, orb%nOrbAtom(iAt)
              tmpir1(ind) = iAtomStart(iAt) + jj - 1
              ind = ind + 1
            end do
          end do
          call append(iOrbRegion, tmpir1)
if (allocated(tmpir1)) then;    deallocate(tmpir1, stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas01/chem/rs151&
&64/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&,906,ioerr, 2);  endif
# 907 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
          write(tmpStr, "(A,'.out')") trim(input%ctrl%RegionLabel(iReg))
          call append(regionLabels, tmpStr)
        end if
if (allocated(iAtomRegion)) then;    deallocate(iAtomRegion, stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas01/&
&chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&,910,ioerr, 2);  endif
# 911 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
      end do
    end if

    tForces = input%ctrl%tForces

    tStress = ((tPeriodic .and. tForces).and.tStress) ! requires stress to
    !  already be possible and it being a periodic calculation with forces

    if (input%ctrl%nrChrg == 0.0_dp .and. (.not.tPeriodic) .and. tMulliken) then
      tDipole = .true.
    else
      tDipole = .false.
    end if

    nMovedAtom = input%ctrl%nrMoved
    nMovedCoord = 3 * nMovedAtom
    nGeoSteps = input%ctrl%maxRun
    if (nMovedAtom > 0) then
if (.not. allocated(indMovedAtom)) then;    allocate(indMovedAtom   (size(input%ctrl%indMovedAtom)), stat=ioerr);    if (ioerr /= &
&0) call allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 929, ioerr, 1);  else;    ca&
&ll allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 929, "Array already allocated", 1&
&);  endif
# 930 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
      indMovedAtom(:) = input%ctrl%indMovedAtom(:)
    else
if (.not. allocated(indMovedAtom)) then;    allocate(indMovedAtom   (0), stat=ioerr);    if (ioerr /= 0) call allocateError("/panf&
&s/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 932, ioerr, 1);  else;    call allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 932&
&, "Array already allocated", 1);  endif
# 933 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
    end if

    nullify(pGeoCoordOpt)
    if (tCoordOpt) then
if (.not. allocated(tmpCoords)) then;    allocate(tmpCoords  (nMovedCoord), stat=ioerr);    if (ioerr /= 0) call allocateError("/p&
&anfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 937, ioerr, 1);  else;    call allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 937&
&, "Array already allocated", 1);  endif
# 938 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
      tmpCoords(1:nMovedCoord) = reshape(coord0(:, indMovedAtom), &
          & (/ nMovedCoord /))
      select case (input%ctrl%iGeoOpt)
      case(1)
if (.not. allocated(tmpWeight)) then;    allocate(tmpWeight  (nMovedCoord), stat=ioerr);    if (ioerr /= 0) call allocateError("/p&
&anfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 942, ioerr, 1);  else;    call allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 942&
&, "Array already allocated", 1);  endif
# 943 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
        tmpWeight(1:nMovedCoord) = 0.5_dp * deltaT**2 &
            & / reshape(spread(mass(specie0(indMovedAtom(:))), 1, 3), &
            & (/nMovedCoord/))
        call create(pSteepDesc, size(tmpCoords), input%ctrl%maxForce, 0.2_dp,&
            & tmpWeight )
if (allocated(tmpWeight)) then;    deallocate(tmpWeight, stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas01/chem&
&/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&,948,ioerr, 2);  endif
# 949 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
        call create(pGeoCoordOpt, pSteepDesc)
      case (2)
        call create(pConjGrad, size(tmpCoords), input%ctrl%maxForce, 0.2_dp)
        call create(pGeoCoordOpt, pConjGrad)
      end select
      call reset(pGeoCoordOpt, tmpCoords)
    end if

    nullify(pGeoLatOpt)
    if (tLatOpt) then
      select case (input%ctrl%iGeoOpt)
      case(1)
if (.not. allocated(tmpWeight)) then;    allocate(tmpWeight  (9), stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panas&
&as01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 961, ioerr, 1);  else;    call allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 961&
&, "Array already allocated", 1);  endif
# 962 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
        tmpWeight = 1.0_dp
        call create(pSteepDescLat, 9, input%ctrl%maxForce, 0.2_dp,&
            & tmpWeight )
if (allocated(tmpWeight)) then;    deallocate(tmpWeight, stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas01/chem&
&/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&,965,ioerr, 2);  endif
# 966 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
        call create(pGeoLatOpt, pSteepDescLat)
      case(2)
        call create(pConjGradLat, 9, input%ctrl%maxForce, 0.2_dp)
        call create(pGeoLatOpt, pConjGradLat)
      end select
      if (tLatOptIsotropic ) then ! optimization uses scaling factor of lattice
        !  vectors
        call reset( pGeoLatOpt, &
            &(/1.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp/))
      else if (tLatOptFixAng) then
        call reset( pGeoLatOpt, & ! optimization uses scaling factor of unit
            !  cell
            &(/1.0_dp,1.0_dp,1.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp/))
      else
        call reset( pGeoLatOpt, reshape(latVec, (/ 9 /)) )
      end if
    end if

    if (.not.(tGeoOpt.or.tMD)) then
      nGeoSteps = 0
    end if

    !! Initialize constraints
    nGeoConstr = input%ctrl%nrConstr
    if (nGeoConstr > 0) then
if (.not. allocated(conAtom)) then;    allocate(conAtom   (input%ctrl%nrConstr), stat=ioerr);    if (ioerr /= 0) call allocateErro&
&r("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 991, ioerr, 1);  else;    call allocateError("/p&
&anfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 991, "Array already allocated", 1);  endif
# 992 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (.not. allocated(conVec)) then;    allocate(conVec   (3, input%ctrl%nrConstr), stat=ioerr);    if (ioerr /= 0) call allocateErr&
&or("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 992, ioerr, 1);  else;    call allocateError("/&
&panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 992, "Array already allocated", 1);  endif
# 993 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"


      conAtom(:) = input%ctrl%conAtom(:)
      conVec(:,:) = input%ctrl%conVec(:,:)
      do ii = 1, nGeoConstr
        conVec(:,ii) = conVec(:,ii) / sqrt(sum(conVec(:,ii)**2))
      end do
    else
if (.not. allocated(conAtom)) then;    allocate(conAtom   (0), stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas0&
&1/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 1001, ioerr, 1);  else;    call allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 10&
&01, "Array already allocated", 1);  endif
# 1002 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (.not. allocated(conVec)) then;    allocate(conVec   (3, 0), stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas&
&01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 1002, ioerr, 1);  else;    call allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 10&
&02, "Array already allocated", 1);  endif
# 1003 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
    end if

    !! Dispersion
    tDispersion = associated(input%ctrl%pDispInp)
    if (tDispersion) then
      select case (input%ctrl%pDispInp%iDisp)
      case (iSlaterKirkwood)
        tStress = .false.
        if (tLatOpt) then
          call error("Sorry, lattice optimisation and Slater-Kirkwood type &
              &dispersion can not be used together")
        end if
        if (tBarostat) then
          call error("Sorry, barostatic MD and Slater-Kirkwood type &
              &dispersion can not be used together")
        end if
nullify(pSlaKirk);  if (.not. associated(pSlaKirk)) then;    allocate(pSlaKirk, stat=ioerr);    if (ioerr /= 0) call allocateError&
&("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 1019, ioerr, 3);  else;    call allocateError("/p&
&anfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 1019, "Ptr. already associated!", 3);  endif
# 1020 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
        if (tPeriodic) then
          call init(pSlaKirk, input%ctrl%pDispInp%pSlaKirkInp, latVec, &
              &recVec, volume)
        else
          call init(pSlaKirk, input%ctrl%pDispInp%pSlaKirkInp)
        end if
        call init(myDispersion, pSlaKirk)
      case (iVdWUFF)
        ! tStress = tStress ! as this is consistent with stress if already
        ! allowed
nullify(pVdWUFF);  if (.not. associated(pVdWUFF)) then;    allocate(pVdWUFF, stat=ioerr);    if (ioerr /= 0) call allocateError("/&
&panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 1030, ioerr, 3);  else;    call allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 10&
&30, "Ptr. already associated!", 3);  endif
# 1031 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
        if (tPeriodic) then
          call init(pVdWUFF, input%ctrl%pDispInp%pVdWUFFInp, nAtom, specie0, &
              &latVec, recVec, volume)
        else
          call init(pVdWUFF, input%ctrl%pDispInp%pVdWUFFInp, nAtom)
        end if
        call init(myDispersion, pVdWUFF)
      end select
      mCutoff = max(mCutoff, getRCutoff(myDispersion))

    end if

    !! Generate a random id for the run. Seed with system time if possible.
    call system_clock(iSeed)
    !! Try date_and_time if system_clock does not work properly.
    if (iSeed < 1) then
      call date_and_time(values=timeValues)
      if (timeValues(5) >= 0) then
        iSeed = 1000 * (60 * (60 * timeValues(5) + timeValues(6)) &
            &+ timeValues(7)) + timeValues(8)
      end if
    end if
    !! Try  default seeder if attempts to use the clock failed.
    if (iSeed < 1) then
      call random_seed
      call random_number(rTmp)
      ! Make sure seed > 0.
      iSeed = int(real(huge(iSeed) - 1, dp) * rTmp) + 1
    end if
    call create(pRanlux, initSeed=iSeed)
    call getRandom(pRanlux, rTmp)
    runId = int(real(huge(runId) - 1, dp) * rTmp) + 1
    call destroy(pRanlux)

    !! Create random generator and pull off first 10
    !! random numbers to avoid disturbing the subsequent sequence.
    !! random_pool can then be used to seed other generators if needed,
    !! with a supply of random numbers controlled from the initial seed.
    !! If the size of random_pool is changed then reproducibility of
    !! the random numbers if initialised from a seed is lost.
    iSeed = input%ctrl%iSeed
    if (iSeed < 1) then
      iSeed = runId     ! No seed specified, use random runId
    end if
    call create(pRanlux, 3, iSeed)
    call getRandom(pRanlux,random_pool(:))


    !! MD stuff
    if (tMD) then
      !! Create MD framework.
      call create(pMDFrame, nMovedAtom, nAtom, input%ctrl%tMDstill)

      !! Create temperature profile, if thermostat is not the dummy one
      if (input%ctrl%iThermostat /= 0) then
        call create(pTempProfile, input%ctrl%tempMethods, &
            &input%ctrl%tempSteps, input%ctrl%tempValues)
      else
        nullify(pTempProfile)
      end if

      !! Create thermostat
      select case (input%ctrl%iThermostat)
      case (0) ! No thermostat
        call create(pDummyTherm, tempAtom, mass(specie0(indMovedAtom)), &
            &pRanlux, pMDFrame)
        call create(pThermostat, pDummyTherm)
      case (1) ! Andersen thermostat
        call create(pAndersenTherm, pRanlux, mass(specie0(indMovedAtom)), &
            &pTempProfile, input%ctrl%tRescale, input%ctrl%wvScale, pMDFrame)
        call create(pThermostat, pAndersenTherm)
      case (2) ! Berendsen thermostat
        call create(pBerendsenTherm, pRanlux, mass(specie0(indMovedAtom)), &
            &pTempProfile, input%ctrl%wvScale, pMDFrame)
        call create(pThermostat, pBerendsenTherm)
      case (3) ! Nose-Hoover-Chain thermostat
        if (input%ctrl%tInitNHC) then
          call create(pNHCTherm, pRanlux, mass(specie0(indMovedAtom)), &
              & pTempProfile, input%ctrl%wvScale, pMDFrame, input%ctrl%deltaT, &
              & input%ctrl%nh_npart, input%ctrl%nh_nys, input%ctrl%nh_nc, &
              & input%ctrl%xnose, input%ctrl%vnose, input%ctrl%gnose)
        else
          call create(pNHCTherm, pRanlux, mass(specie0(indMovedAtom)), &
              &pTempProfile, input%ctrl%wvScale, pMDFrame, input%ctrl%deltaT, &
              & input%ctrl%nh_npart, input%ctrl%nh_nys, input%ctrl%nh_nc)
        end if
        call create(pThermostat, pNHCTherm)
      end select

      !! Create MD integrator
      if (input%ctrl%tReadMDVelocities) then
        if (tBarostat) then
          call create(pVelocityVerlet, deltaT, coord0(:,indMovedAtom),&
              & pThermostat,input%ctrl%initialVelocities, &
              & BarostatStrength,pressure,input%ctrl%tIsotropic)
        else
          call create(pVelocityVerlet, deltaT, coord0(:,indMovedAtom),&
              & pThermostat,input%ctrl%initialVelocities)
        end if
      else
        if (tBarostat) then
          call create(pVelocityVerlet, deltaT, coord0(:,indMovedAtom),&
              & pThermostat, BarostatStrength,pressure,input%ctrl%tIsotropic)
        else
          call create(pVelocityVerlet, deltaT, coord0(:,indMovedAtom),&
              & pThermostat)
        end if
      end if

      call create(pMDIntegrator, pVelocityVerlet)
    else
      nullify(pMDIntegrator)
    end if

    if (tDerivs) then
if (.not. allocated(tmp3Coords)) then;    allocate(tmp3Coords   (3,nMovedAtom), stat=ioerr);    if (ioerr /= 0) call allocateError&
&("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 1146, ioerr, 1);  else;    call allocateError("/p&
&anfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 1146, "Array already allocated", 1);  endif
# 1147 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
      tmp3Coords = coord0(:,indMovedAtom)
      call create(pDerivDriver,tmp3Coords, &
          & input%ctrl%derivDelta)
      coord0(:,indMovedAtom) = tmp3Coords
if (allocated(tmp3Coords)) then;    deallocate(tmp3Coords, stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas01/ch&
&em/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&,1151,ioerr, 2);  endif
# 1152 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
      nGeoSteps = 2 * 3 * nMovedAtom - 1
    end if

    if (input%ctrl%tEField) then
      tEField = .true.
      EFieldStrength = input%ctrl%EFieldStrength
      EfieldVector(:) = input%ctrl%EfieldVector(:)
      tTDEfield = input%ctrl%tTDEfield
      EfieldOmega = input%ctrl%EfieldOmega
      EfieldPhase = input%ctrl%EfieldPhase
      if (tTDEfield .and. .not. tMD) then
        call error ("Time dependent electric fields only possible for MD!")
      end if
      ! parser should catch all of these:

    else
      tEField = .false.
      EFieldStrength = 0.0_dp
      EfieldVector(:) = 0.0_dp
      tTDEfield = .false.
      EfieldOmega = 0.0_dp
      EfieldPhase = 0
    end if

    !! Allocate charge arrays
    if (tMulliken) then ! automatically true if tSCC
if (.not. allocated(q0)) then;    allocate(q0   (orb%mOrb, nAtom, nSpin), stat=ioerr);    if (ioerr /= 0) call allocateError("/pan&
&fs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 1178, ioerr, 1);  else;    call allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 11&
&78, "Array already allocated", 1);  endif
# 1179 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
      q0(:,:,:) = 0.0_dp

if (.not. allocated(qShell0)) then;    allocate(qShell0   (orb%mShell, nAtom), stat=ioerr);    if (ioerr /= 0) call allocateError(&
&"/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 1181, ioerr, 1);  else;    call allocateError("/pa&
&nfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 1181, "Array already allocated", 1);  endif
# 1182 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
      qShell0(:,:) = 0.0_dp
    else
if (.not. allocated(q0)) then;    allocate(q0   (0,0,0), stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas01/chem&
&/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 1184, ioerr, 1);  else;    call allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 11&
&84, "Array already allocated", 1);  endif
# 1185 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (.not. allocated(qShell0)) then;    allocate(qShell0   (0,0), stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasa&
&s01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 1185, ioerr, 1);  else;    call allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 11&
&85, "Array already allocated", 1);  endif
# 1186 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
    end if

if (.not. allocated(qInput)) then;    allocate(qInput   (orb%mOrb, nAtom, nSpin), stat=ioerr);    if (ioerr /= 0) call allocateErr&
&or("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 1188, ioerr, 1);  else;    call allocateError("&
&/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 1188, "Array already allocated", 1);  endif
# 1189 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (.not. allocated(qOutput)) then;    allocate(qOutput   (orb%mOrb, nAtom, nSpin), stat=ioerr);    if (ioerr /= 0) call allocateE&
&rror("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 1189, ioerr, 1);  else;    call allocateError&
&("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 1189, "Array already allocated", 1);  endif
# 1190 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
    qInput(:,:,:) = 0.0_dp
    qOutput(:,:,:) = 0.0_dp

    if (tDFTBU) then
if (.not. allocated(qBlockIn)) then;    allocate(qBlockIn   (orb%mOrb, orb%mOrb, nAtom, nSpin), stat=ioerr);    if (ioerr /= 0) ca&
&ll allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 1194, ioerr, 1);  else;    call a&
&llocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 1194, "Array already allocated", 1); 
 endif
# 1195 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (.not. allocated(qBlockOut)) then;    allocate(qBlockOut   (orb%mOrb, orb%mOrb, nAtom, nSpin), stat=ioerr);    if (ioerr /= 0) &
&call allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 1195, ioerr, 1);  else;    call&
& allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 1195, "Array already allocated", 1)&
&;  endif
# 1196 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
      qBlockIn = 0.0_dp
      qBlockOut = 0.0_dp
      if (tImHam) then
if (.not. allocated(qiBlockIn)) then;    allocate(qiBlockIn   (orb%mOrb, orb%mOrb, nAtom, nSpin), stat=ioerr);    if (ioerr /= 0) &
&call allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 1199, ioerr, 1);  else;    call&
& allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 1199, "Array already allocated", 1)&
&;  endif
# 1200 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
        qiBlockIn = 0.0_dp
      else
if (.not. allocated(qiBlockIn)) then;    allocate(qiBlockIn   (0, 0, 0, 0), stat=ioerr);    if (ioerr /= 0) call allocateError("/p&
&anfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 1202, ioerr, 1);  else;    call allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 12&
&02, "Array already allocated", 1);  endif
# 1203 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
        qiBlockIn = 0.0_dp
      end if
    else
if (.not. allocated(qBlockIn)) then;    allocate(qBlockIn   (0, 0, 0, 0), stat=ioerr);    if (ioerr /= 0) call allocateError("/pan&
&fs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 1206, ioerr, 1);  else;    call allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 12&
&06, "Array already allocated", 1);  endif
# 1207 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (.not. allocated(qBlockOut)) then;    allocate(qBlockOut   (0, 0, 0, 0), stat=ioerr);    if (ioerr /= 0) call allocateError("/p&
&anfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 1207, ioerr, 1);  else;    call allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 12&
&07, "Array already allocated", 1);  endif
# 1208 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (.not. allocated(qiBlockIn)) then;    allocate(qiBlockIn   (0, 0, 0, 0), stat=ioerr);    if (ioerr /= 0) call allocateError("/p&
&anfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 1208, ioerr, 1);  else;    call allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 12&
&08, "Array already allocated", 1);  endif
# 1209 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
      qiBlockIn = 0.0_dp
      qBlockIn = 0.0_dp
      qBlockOut = 0.0_dp
    end if

    if (tImHam) then
if (.not. allocated(qiBlockOut)) then;    allocate(qiBlockOut   (orb%mOrb, orb%mOrb, nAtom, nSpin), stat=ioerr);    if (ioerr /= 0&
&) call allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 1215, ioerr, 1);  else;    ca&
&ll allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 1215, "Array already allocated", &
&1);  endif
# 1216 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
      qiBlockOut = 0.0_dp
    end if

    if (tSCC) then
if (.not. allocated(qDiffRed)) then;    allocate(qDiffRed   (nMixElements), stat=ioerr);    if (ioerr /= 0) call allocateError("/p&
&anfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 1220, ioerr, 1);  else;    call allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 12&
&20, "Array already allocated", 1);  endif
# 1221 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (.not. allocated(qInpRed)) then;    allocate(qInpRed   (nMixElements), stat=ioerr);    if (ioerr /= 0) call allocateError("/pan&
&fs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 1221, ioerr, 1);  else;    call allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 12&
&21, "Array already allocated", 1);  endif
# 1222 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (.not. allocated(qOutRed)) then;    allocate(qOutRed   (nMixElements), stat=ioerr);    if (ioerr /= 0) call allocateError("/pan&
&fs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 1222, ioerr, 1);  else;    call allocateError("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 12&
&22, "Array already allocated", 1);  endif
# 1223 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
      qDiffRed = 0.0_dp
      qInpRed = 0.0_dp
      qOutRed = 0.0_dp
    end if

    !! Initialize Mulliken charges
    if (tMulliken) then
      call initQFromShellChrg(q0, referenceN0, specie0, orb)
    end if


    tReadChrg = input%ctrl%tReadChrg
    if (tSCC) then
      do iAt = 1, nAtom
        iSp = specie0(iAt)
        do iSh = 1, orb%nShell(iSp)
          qShell0 (iSh,iAt) = sum(q0(orb%posShell(iSh,iSp): &
              & orb%posShell(iSh+1,iSp)-1,iAt,1))
        end do
      end do
      if (tReadChrg) then
        if (tDFTBU) then
          if (nSpin == 2) then
            call initQFromFile(qInput, fChargeIn, sum(nEl), orb, &
                & magnetisation=nEl(1)-nEl(2), qBlock=qBlockIn)
          else
            if (tImHam) then
              call initQFromFile(qInput, fChargeIn, nEl(1), orb, &
                  & qBlock=qBlockIn,qiBlock=qiBlockIn)
            else
              call initQFromFile(qInput, fChargeIn, nEl(1), orb, &
                  & qBlock=qBlockIn)
            end if
          end if
        else
          ! hack again caused by going from up/down to q and M
          if (nSpin == 2) then
            call initQFromFile(qInput, fChargeIn, sum(nEl), orb, &
                & magnetisation=nEl(1)-nEl(2))
          else
            call initQFromFile(qInput, fChargeIn, nEl(1), orb)
          end if
        end if
      else
        if (associated(input%ctrl%initialCharges)) then
          if (abs(sum(input%ctrl%initialCharges) - input%ctrl%nrChrg) &
              &> 1e-4_dp) then
            write(strTmp, "(A,G13.6,A,G13.6,A)") "Sum of initial charges does &
                &not match specified total charge. (", &
                &sum(input%ctrl%initialCharges), " vs. ", &
                &input%ctrl%nrChrg, ")"
            call warning(strTmp)
          end if
          call initQFromAtomChrg(qInput, input%ctrl%initialCharges, &
              &referenceN0, specie0, specieName, orb)
        else
          qInput(:,:,:) = q0(:,:,:)
        end if

        select case (nSpin)
        case (1)
          ! nothing to do
        case (2)
          if (associated(input%ctrl%initialSpins)) then
            do ii = 1, nAtom
              !! doesn't actually matter if additional spin
              !! polarization pushes charges to <0 as the initial charges aren't
              !! mixed in to later iterations
              qInput(1:orb%nOrbAtom(ii),ii,2) = &
                  & qInput(1:orb%nOrbAtom(ii),ii,1) * &
                  & input%ctrl%initialSpins(1,ii) / &
                  & sum(qInput(1:orb%nOrbAtom(ii),ii,1))
            end do
          else
            do ii = 1, nAtom
              qInput(1:orb%nOrbAtom(ii),ii,2) = &
                  & qInput(1:orb%nOrbAtom(ii),ii,1) * &
                  & (nEl(1)-nEl(2))/sum(qInput(:,:,1))
            end do
          end if
        case (4)
          if (tSpin) then
            if (.not.associated(input%ctrl%initialSpins)) then
              call error("Missing initial spins!")
            end if
            if (any(shape(input%ctrl%initialSpins)/=(/3,nAtom/))) then
              call error("Incorrect shape initialSpins array!")
            end if
            do ii = 1, nAtom
              do jj = 1, 3
                qInput(1:orb%nOrbAtom(ii),ii,jj+1) = &
                    & qInput(1:orb%nOrbAtom(ii),ii,1) * &
                    & input%ctrl%initialSpins(jj,ii) &
                    & / sum(qInput(1:orb%nOrbAtom(ii),ii,1))
              end do
            end do
          end if
        end select
        if (tDFTBU) then
          qBlockIn = 0.0_dp
          do iS = 1, nSpin
            do iAt = 1, nAtom
              iSp = specie0(iAt)
              do iSh = 1, orb%nShell(iSp)
                iStart = orb%posShell(iSh,iSp)
                iEnd = orb%posShell(iSh+1,iSp)-1
                rTmp = sum(qInput(iStart:iEnd,iAt,iS))
                rTmp = rTmp / real(iEnd+1-iStart,dp)
                do ii = iStart, iEnd
                  qBlockIn(ii,ii,iAt,iS) = rTmp
                end do
              end do
            end do
          end do
          if (tImHam) then
            qiBlockIn = 0.0_dp
          end if
        end if
      end if

      qInpRed = 0.0_dp
      if (nSpin == 2) then
        call qm2ud(qInput)
        if (tDFTBU) then
          call qm2ud(qBlockIn)
        end if
      end if

      call OrbitalEquiv_reduce(qInput, iEqOrbitals, orb, qInpRed(1:nIneqOrb))
      if (tDFTBU) then
        call AppendBlock_reduce( qBlockIn,iEqBlockDFTBU, orb, &
            & qInpRed )
        if (tImHam) then
          call AppendBlock_reduce( qiBlockIn,iEqBlockDFTBULS, orb, &
              & qInpRed, skew=.true. )
        end if
      end if

      if (nSpin == 2) then
        call ud2qm(qInput)
        if (tDFTBU) then
          call ud2qm(qBlockIn)
        end if
      end if
    end if

    !! Set various options
    tWriteTagged = input%ctrl%tWriteTagged
    tWriteDetailedXML = input%ctrl%tWriteDetailedXML
    tWriteResultsTag = input%ctrl%tWriteResultsTag
    tWriteDetailedOut = input%ctrl%tWriteDetailedOut
    tWriteBandDat = input%ctrl%tWriteBandDat
    tWriteHS = input%ctrl%tWriteHS
    tWriteRealHS = input%ctrl%tWriteRealHS

    !! Minimize memory usage?
    tMinMemory = input%ctrl%tMinMemory
    tStoreEigvecs = tMinMemory .and. (nKPoint > 1 .or. nSpin == 2 )
    if (tStoreEigvecs) then
      if (tRealHS.and.(.not.t2Component)) then
if (.not. allocated(storeEigvecsReal)) then;    allocate(storeEigvecsReal   (nSpin), stat=ioerr);    if (ioerr /= 0) call allocate&
&Error("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 1383, ioerr, 1);  else;    call allocateErro&
&r("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 1383, "Array already allocated", 1);  endif
# 1384 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (.not. allocated(storeEigvecsCplx)) then;    allocate(storeEigvecsCplx  (0), stat=ioerr);    if (ioerr /= 0) call allocateError&
&("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 1384, ioerr, 1);  else;    call allocateError("/p&
&anfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 1384, "Array already allocated", 1);  endif
# 1385 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
        do iS = 1, nSpin
          call init(storeEigvecsReal(iS), 0, "tmp_eigvr_")
        end do
      else
        if (t2Component) then
if (.not. allocated(storeEigvecsCplx)) then;    allocate(storeEigvecsCplx   (1), stat=ioerr);    if (ioerr /= 0) call allocateErro&
&r("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 1390, ioerr, 1);  else;    call allocateError("/&
&panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 1390, "Array already allocated", 1);  endif
# 1391 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (.not. allocated(storeEigvecsReal)) then;    allocate(storeEigvecsReal  (0), stat=ioerr);    if (ioerr /= 0) call allocateError&
&("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 1391, ioerr, 1);  else;    call allocateError("/p&
&anfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 1391, "Array already allocated", 1);  endif
# 1392 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
          call init(storeEigvecsCplx(1), 0, "tmp_eigvc_")
        else
          write(*,*)'Here'
if (.not. allocated(storeEigvecsCplx)) then;    allocate(storeEigvecsCplx   (nSpin), stat=ioerr);    if (ioerr /= 0) call allocate&
&Error("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 1395, ioerr, 1);  else;    call allocateErro&
&r("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 1395, "Array already allocated", 1);  endif
# 1396 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (.not. allocated(storeEigvecsReal)) then;    allocate(storeEigvecsReal  (0), stat=ioerr);    if (ioerr /= 0) call allocateError&
&("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 1396, ioerr, 1);  else;    call allocateError("/p&
&anfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 1396, "Array already allocated", 1);  endif
# 1397 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
          do iS = 1, nSpin
            call init(storeEigvecsCplx(iS), 0, "tmp_eigvc_")
          end do
        end if
      end if
    else
if (.not. allocated(storeEigvecsReal)) then;    allocate(storeEigvecsReal  (0), stat=ioerr);    if (ioerr /= 0) call allocateError&
&("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 1403, ioerr, 1);  else;    call allocateError("/p&
&anfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 1403, "Array already allocated", 1);  endif
# 1404 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (.not. allocated(storeEigvecsCplx)) then;    allocate(storeEigvecsCplx  (0), stat=ioerr);    if (ioerr /= 0) call allocateError&
&("/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90", 1404, ioerr, 1);  else;    call allocateError("/p&
&anfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 1404, "Array already allocated", 1);  endif
# 1405 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
    end if

    !! Check if stopfiles already exist and quit if yes
    inquire(file=fStopSCC, exist=tExist)
    if (tExist) then
      call error("Stop file '" // fStopSCC // "' already present at startup")
    end if
    inquire(file=fStopDriver, exist=tExist)
    if (tExist) then
      call error("Stop file '" // fStopDriver // "' already present at startup")
    end if

    restartFreq = input%ctrl%restartFreq

    tInitialized = .true.


    if (input%ctrl%tMD) then
      select case(input%ctrl%iThermostat)
      case (0)
        write(*,"('Mode:',T30,A,/,T30,A)") 'MD without scaling of velocities',&
            & '(a.k.a. NVE ensemble)'
      case (1)
        write(*,"('Mode:',T30,A,/,T30,A)") &
            & "MD with re-selection of velocities according to temperature", &
            & "(a.k.a. NVT ensemble using Andersen thermostating)"
      case(2)
        write(*,"('Mode:',T30,A,/,T30,A)") &
            & "MD with scaling of velocities according to temperature", &
            & "(a.k.a. not NVT ensemble using Berendsen thermostating)"
      case(3)
        write(*,"('Mode:',T30,A,/,T30,A)") &
            & "MD with scaling of velocities according to", &
            & "Nose-Hoover-Chain thermostat"

      case default
        call error("Unknown thermostat mode")
      end select
    elseif (tGeoOpt) then
      if (nGeoConstr > 0) then
        strTmp = "with constraints"
      else
        strTmp = ""
      end if
      select case (input%ctrl%iGeoOpt)
      case (1)
        write(*,"('Mode:',T30,A)")'Steepest descent' // trim(strTmp)
      case (2)
        write(*,"('Mode:',T30,A)") 'Conjugate gradient relaxation' &
            &// trim(strTmp)
      case default
        call error("Unknown optimisation mode")
      end select
    elseif (tDerivs) then
      write (*,"('Mode:',T30,A)") "2nd derivatives calculation"
      write (*,"('Mode:',T30,A)") "Calculated for atoms:"
      write (*,*) indMovedAtom
    else
      write (*,"('Mode:',T30,A)") "Static calculation"
    end if

    if (tSCC) then
      write (*, "(A,':',T30,A)") "Self consistent charges", "Yes"
      write (*, "(A,':',T30,E14.6)") "SCC-tolerance", sccTol
      write (*, "(A,':',T30,I14)") "Max. scc iterations", nSCCIter
      write (*, "(A,':',T30,E14.6)") "Ewald alpha parameter", getSCCEwaldPar()
      if (tDFTBU) then
        write (*, "(A,':',T35,A)") "Orbitally dependant functional", "Yes"
        write (*, "(A,':',T30,I14)") "Orbital functional number",nDFTBUfunc !
        !  use module to reverse look up name
      end if
    else
      write (*, "(A,':',T30,A)") "Self consistent charges", "No"
    end if

    select case (nSpin)
    case(1)
      write (*,"(A,':',T30,A)") "Spin polarisation", "No"
      write (*, "(A,':',T30,F12.6,/,A,':',T30,F12.6)") "Nr. of up electrons", &
          &0.5_dp*nEl(1), "Nr. of down electrons", 0.5_dp*nEl(1)
    case(2)
      write (*,"(A,':',T30,A)") "Spin polarisation", "Yes"
      write (*, "(A,':',T30,F12.6,/,A,':',T30,F12.6)") "Nr. of up electrons", &
          &nEl(1), "Nr. of down electrons", nEl(2)
    case(4)
      write (*,"(A,':',T30,A)") "Non-collinear calculation", "Yes"
      write (*, "(A,':',T30,F12.6)") "Nr. of electrons", nEl(1)
    end select

    if (tPeriodic) then
      write (*, "(A,':',T30,A)") "Periodic boundaries", "Yes"
      if (tLatOpt) then
        write (*, "(A,':',T30,A)") "Lattice optimisation", "Yes"
        write (*, "(A,':',T30,f12.6)") "Pressure", pressure
      end if
    else
      write (*, "(A,':',T30,A)") "Periodic boundaries", "No"
    end if

    select case (solver)
    case(1)
      write (strTmp, "(A)") "Standard"
    case(2)
      write (strTmp, "(A)") "Divide and Conquer"
    case(3)
      write (strTmp, "(A)") "Relatively robust (version 1)"
    case(4)
      write (strTmp, "(A)") "Relativel robust (version 2)"
    case default
      call error("Unknown eigensolver!")
    end select
    write (*, "(A,':',T30,A)") "Diagonalizer", trim(strTmp)

    if (tSCC) then
      select case (iMixer)
      case(1)
        write (strTmp, "(A)") "Simple"
      case(2)
        write (strTmp, "(A)") "Anderson"
      case(3)
        write (strTmp, "(A)") "Broyden"
      case(4)
        write (strTmp, "(A)") "DIIS"
      end select
      write (*, "(A,':',T30,A,' ',A)") "Mixer", trim(strTmp), "mixer"
      write (*, "(A,':',T30,F14.6)") "Mixing parameter", mixParam
      write (*, "(A,':',T30,I14)") "Maximal SCC-cycles", nSCCIter
      select case (iMixer)
      case(2)
        write (*, "(A,':',T30,I14)") "Nr. of chrg. vectors to mix", nGeneration
      case(3)
        write (*, "(A,':',T30,I14)") "Nr. of chrg. vec. in memory", nGeneration
      case(4)
        write (*, "(A,':',T30,I14)") "Nr. of chrg. vectors to mix", nGeneration
      end select
    end if

    if (tCoordOpt) then
      write (*, "(A,':',T30,I14)") "Nr. of moved atoms", nMovedAtom
    end if
    if (tGeoOpt) then
      write (*, "(A,':',T30,I14)") "Max. nr. of geometry steps", nGeoSteps
      write (*, "(A,':',T30,E14.6)") "Force tolerance", input%ctrl%maxForce
      if (input%ctrl%iGeoOpt == 1) then
        write (*, "(A,':',T30,E14.6)") "Step size", deltaT
      end if
    end if


    tFirst = .true.
    if (nGeoConstr > 0) then
      do ii = 1, nAtom
        do jj = 1, nGeoConstr
          if (conAtom(jj) == ii) then
            if (tFirst) then
              write(strTmp, "(A,':')") "Geometry constraints"
              tFirst = .false.
            else
              write(strTmp, "(A)") ""
            end if
            write(*, "(A,T30,'At',I4,': ',3F10.6)") trim(strTmp), &
                &ii, (conVec(kk,jj), kk=1,3)
          end if
        end do
      end do
    end if

    if (.not.input%ctrl%tSetFillingTemp) then
      write (*, "(A,':',T30,E14.6)") "Electronic temperature", tempElec
    end if
    if (tMD) then
      write (*, "(A,':',T30,E14.6)") "Time step", deltaT
      if (input%ctrl%iThermostat == 0 .and. .not.input%ctrl%tReadMDVelocities)&
          & then
        write (*, "(A,':',T30,E14.6)") "Temperature", tempAtom
      end if
      write (*, "(A,':',T30,I14)") "Random seed", iSeed
      if (input%ctrl%tRescale) then
        write (*, "(A,':',T30,E14.6)") "Rescaling probability", &
            &input%ctrl%wvScale
      end if
    end if

    if (tSCC) then
      if (input%ctrl%tReadChrg) then
        write (strTmp, "(A,A,A)") "Read in from '", trim(fChargeIn), "'"
      else
        write (strTmp, "(A,E11.3,A)") "Set automatically (system chrg: ", &
            &input%ctrl%nrChrg, ")"
      end if
      write (*, "(A,':',T30,A)") "Initial charges", trim(strTmp)
    end if

    do ii = 1, nType
      if (ii == 1) then
        write (strTmp, "(A,':')") "Included shells"
      else
        write (strTmp, "(A)") ""
      end if
      do jj = 1, orb%nShell(ii)
        if (jj == 1) then
          strTmp2 = trim(orbitalNames(orb%angShell(jj, ii) + 1))
        else
          strTmp2 = trim(strTmp2) // ", " &
              &// trim(orbitalNames(orb%angShell(jj, ii) + 1))
        end if
      end do
      write (*, "(A,T30,A2,':  ',A)") trim(strTmp), trim(specieName(ii)), &
          &trim(strTmp2)
    end do

    if (tPeriodic) then
      do ii = 1, nKPoint
        if (ii == 1) then
          write(strTmp, "(A,':')") "K-points and weights"
        else
          write(strTmp, "(A)") ""
        end if
        write(*,"(A,T28,I6,':',3F10.6,3X,F10.6)") trim(strTmp), ii, &
            & (kPoint(jj, ii), jj=1, 3), kWeight(ii)
      end do
    end if

    if (tDispersion) then
      select case(input%ctrl%pDispInp%iDisp)
      case (iSlaterKirkwood)
        write(*,"(A)") "Using Slater-Kirkwood dispersion corrections"
      case (iVdWUFF)
        write(*,"(A)") "Using Lennard-Jones dispersion corrections"
      case default
        call error("Unknown dispersion model - this shouldn't happen!")
      end select
    end if

    tFirst = .true.
    if (tSpin) then
      do ii = 1, nType
        do jj = 1, orb%nShell(ii)
          do kk = 1, orb%nShell(ii)
            if (tFirst) then
              write(strTmp, "(A,':')") "Spin coupling constants"
              tFirst = .false.
            else
              write(strTmp, "(A)") ""
            end if
            write(*, "(A,T30,A2,2X,I1,'(',A1,')-',I1,'(',A1,'): ',E14.6)") &
                &trim(strTmp), specieName(ii), &
                jj, orbitalNames(orb%angShell(jj, ii)+1), &
                &kk, orbitalNames(orb%angShell(kk, ii)+1), &
                &W(kk, jj, ii)
          end do
        end do
      end do
    end if

    tFirst = .true.
    if (tSpinOrbit) then
      if (tDualSpinOrbit) then
        write(*,"(A)")"Dual representation spin orbit"
      end if
      do ii = 1, nType
        do jj = 1, orb%nShell(ii)
          if (tFirst) then
            write(strTmp, "(A,':')") "Spin orbit constants"
            tFirst = .false.
          else
            write(strTmp, "(A)") ""
          end if
          write(*, "(A,T30,A2,2X,I1,'(',A1,'): ',E14.6)") &
                &trim(strTmp), specieName(ii), &
                jj, orbitalNames(orb%angShell(jj, ii)+1), &
                &xi(jj, ii)
          if (xi(jj, ii) /= 0.0_dp .and. orb%angShell(jj, ii) == 0) then
            call error("Program halt due to non-zero s-orbital spin-orbit &
                &coupling constant!")
          end if
        end do
      end do
    end if

    if (tSCC) then
      if (t3rdFull) then
        write (*, "(A,T30,A)") "Full 3rd order correction", "Yes"
      end if
      if (any(tDampedShort)) then
        write(*, "(A,T30,A)") "Damped SCC", "Yes"
        ii = count(tDampedShort)
        write(strTmp, "(A,I0,A)") "(A,T30,", ii, "(A,1X))"
        write(*, strTmp) "Damped specie(s):", pack(specieName, tDampedShort)
if (allocated(tDampedShort)) then;    deallocate(tDampedShort, stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas0&
&1/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&,1694,ioerr, 2);  endif
# 1695 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
      end if
    end if

    write (*, "(A,':')") "Extra options"
    if (tMulliken .and. .not. tSCC) then
      write (*, "(T30,A)") "Mulliken analysis"
    end if
    if (tForces .and. .not. (tMD .or. tGeoOpt .or. tDerivs)) then
      write (*, "(T30,A)") "Force calculation"
    end if
    if (tEigenvecs) then
      write (*, "(T30,A)") "Eigenvector printing"
    end if
    if (tExtChrg) then
      write (*, "(T30,A)") "External charges specified"
    end if

    if (tEField) then
      if (.not.tSCC) then
        call error("External electric field currently unsuported for non-SCC.")
      end if
      if (tTDEfield) then
        write (*, "(T30,A)") "External electric field specified"
        write (*, "(A,':',T30,E14.6)") "Angular frequency", EfieldOmega
      else
        write (*, "(T30,A)") "External static electric field specified"
      end if
      write (*, "(A,':',T30,E14.6)") "Field strength", EFieldStrength
      write (*, "(A,':',T30,3F9.6)") "Direction", EfieldVector
      if (tPeriodic) then
        call warning("Saw tooth potential used for periodic geometry &
            &- make sure there is a vacuum region!")
      end if
    end if

    if (tDFTBU) then
      do ii = 1, nType
        if (nUJ(ii)>0) then
          write(strTmp, "(A,':')") "U-J coupling constants"
          write(*, "(A,T25,A2)")trim(strTmp), specieName(ii)
          do jj = 1, nUJ(ii)
            write(strTmp, "(A,I1,A)")'(A,',niUJ(jj,ii),'I2,T25,A,F6.4)'
            write(*,trim(strTmp))'Shells:',iUJ(1:niUJ(jj,ii),jj,ii),'UJ:', &
                & UJ(jj,ii)
          end do
        end if
      end do

    end if

    if (.not.tStress) then
      if (tBarostat) then
        call error("Sorry, MD with a barost requires stress stress evaluation")
      end if
      if (tLatOpt) then
        call error("Sorry, lattice optimization requires stress tensor evaluation")
      end if
    end if

    if (tSpinOrbit .and. (tWriteHS .or. tWriteRealHS)) then
      call error("Writing of Hamiltonian currently not possible with spin orbit&
          & coupling enabled.")
    end if

    tInitialized = .true.

  end subroutine initProgramVariables





  !!* Destroys the program variables
  subroutine destroyProgramVariables

    integer :: iS



    tInitialized = .false.

if (allocated(latVec)) then;    deallocate(latVec, stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas01/chem/rs151&
&64/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&,1776,ioerr, 2);  endif
# 1777 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (allocated(recVec)) then;    deallocate(recVec, stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas01/chem/rs151&
&64/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&,1777,ioerr, 2);  endif
# 1778 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (allocated(recVec2p)) then;    deallocate(recVec2p, stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas01/chem/r&
&s15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&,1778,ioerr, 2);  endif
# 1779 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (allocated(hubbU)) then;    deallocate(hubbU, stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas01/chem/rs15164&
&/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&,1779,ioerr, 2);  endif
# 1780 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (allocated(atomEigVal)) then;    deallocate(atomEigVal, stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas01/ch&
&em/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&,1780,ioerr, 2);  endif
# 1781 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (allocated(referenceN0)) then;    deallocate(referenceN0, stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas01/&
&chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&,1781,ioerr, 2);  endif
# 1782 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (allocated(mass)) then;    deallocate(mass, stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas01/chem/rs15164/D&
&FTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&,1782,ioerr, 2);  endif
# 1783 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (associated(Img2CentCell)) then;    deallocate(Img2CentCell, stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas&
&01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 1783, ioerr, 6);    nullify(Img2CentCell);  endif
# 1784 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (associated(cellVec)) then;    deallocate(cellVec, stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas01/chem/rs&
&15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 1784, ioerr, 6);    nullify(cellVec);  endif
# 1785 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (associated(rCellVec)) then;    deallocate(rCellVec, stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas01/chem/&
&rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 1785, ioerr, 6);    nullify(rCellVec);  endif
# 1786 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (associated(iCellVec)) then;    deallocate(iCellVec, stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas01/chem/&
&rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 1786, ioerr, 6);    nullify(iCellVec);  endif
# 1787 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (associated(coord)) then;    deallocate(coord, stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas01/chem/rs1516&
&4/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 1787, ioerr, 6);    nullify(coord);  endif
# 1788 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (allocated(coord0)) then;    deallocate(coord0, stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas01/chem/rs151&
&64/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&,1788,ioerr, 2);  endif
# 1789 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (associated(specie)) then;    deallocate(specie, stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas01/chem/rs15&
&164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 1789, ioerr, 6);    nullify(specie);  endif
# 1790 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (allocated(specie0)) then;    deallocate(specie0, stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas01/chem/rs1&
&5164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&,1790,ioerr, 2);  endif
# 1791 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (associated(iPair)) then;    deallocate(iPair, stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas01/chem/rs1516&
&4/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 1791, ioerr, 6);    nullify(iPair);  endif
# 1792 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (allocated(iAtomStart)) then;    deallocate(iAtomStart, stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas01/ch&
&em/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&,1792,ioerr, 2);  endif
# 1793 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (associated(ham)) then;    deallocate(ham, stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas01/chem/rs15164/DF&
&TB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 1793, ioerr, 6);    nullify(ham);  endif
# 1794 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (associated(iHam)) then;    deallocate(iHam, stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas01/chem/rs15164/&
&DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 1794, ioerr, 6);    nullify(iHam);  endif
# 1795 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (associated(over)) then;    deallocate(over, stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas01/chem/rs15164/&
&DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&, 1795, ioerr, 6);    nullify(over);  endif
# 1796 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (allocated(W)) then;    deallocate(W, stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas01/chem/rs15164/DFTB+/d&
&ftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&,1796,ioerr, 2);  endif
# 1797 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (allocated(chargePerShell)) then;    deallocate(chargePerShell, stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/pana&
&sas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&,1797,ioerr, 2);  endif
# 1798 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (allocated(UJ)) then;    deallocate(UJ, stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas01/chem/rs15164/DFTB+&
&/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&,1798,ioerr, 2);  endif
# 1799 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (allocated(nUJ)) then;    deallocate(nUJ, stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas01/chem/rs15164/DFT&
&B+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&,1799,ioerr, 2);  endif
# 1800 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (allocated(niUJ)) then;    deallocate(niUJ, stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas01/chem/rs15164/D&
&FTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&,1800,ioerr, 2);  endif
# 1801 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (allocated(iUJ)) then;    deallocate(iUJ, stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas01/chem/rs15164/DFT&
&B+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&,1801,ioerr, 2);  endif
# 1802 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (allocated(kPoint)) then;    deallocate(kPoint, stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas01/chem/rs151&
&64/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&,1802,ioerr, 2);  endif
# 1803 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (allocated(kWeight)) then;    deallocate(kWeight, stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas01/chem/rs1&
&5164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&,1803,ioerr, 2);  endif
# 1804 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (allocated(nEl)) then;    deallocate(nEl, stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas01/chem/rs15164/DFT&
&B+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&,1804,ioerr, 2);  endif
# 1805 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (allocated(specieName)) then;    deallocate(specieName, stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas01/ch&
&em/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&,1805,ioerr, 2);  endif
# 1806 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (allocated(indMovedAtom)) then;    deallocate(indMovedAtom, stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas0&
&1/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&,1806,ioerr, 2);  endif
# 1807 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (allocated(conAtom)) then;    deallocate(conAtom, stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas01/chem/rs1&
&5164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&,1807,ioerr, 2);  endif
# 1808 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (allocated(conVec)) then;    deallocate(conVec, stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas01/chem/rs151&
&64/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&,1808,ioerr, 2);  endif
# 1809 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (allocated(qOutput)) then;    deallocate(qOutput, stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas01/chem/rs1&
&5164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&,1809,ioerr, 2);  endif
# 1810 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (allocated(qInput)) then;    deallocate(qInput, stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas01/chem/rs151&
&64/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&,1810,ioerr, 2);  endif
# 1811 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (allocated(qBlockIn)) then;    deallocate(qBlockIn, stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas01/chem/r&
&s15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&,1811,ioerr, 2);  endif
# 1812 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (allocated(qBlockOut)) then;    deallocate(qBlockOut, stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas01/chem&
&/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&,1812,ioerr, 2);  endif
# 1813 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (allocated(qiBlockIn)) then;    deallocate(qiBlockIn, stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas01/chem&
&/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&,1813,ioerr, 2);  endif
# 1814 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (allocated(qiBlockOut)) then;    deallocate(qiBlockOut, stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas01/ch&
&em/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&,1814,ioerr, 2);  endif
# 1815 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (allocated(q0)) then;    deallocate(q0, stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas01/chem/rs15164/DFTB+&
&/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&,1815,ioerr, 2);  endif
# 1816 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (allocated(qDiffRed)) then;    deallocate(qDiffRed, stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas01/chem/r&
&s15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&,1816,ioerr, 2);  endif
# 1817 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (allocated(qInpRed)) then;    deallocate(qInpRed, stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas01/chem/rs1&
&5164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&,1817,ioerr, 2);  endif
# 1818 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (allocated(qOutRed)) then;    deallocate(qOutRed, stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas01/chem/rs1&
&5164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&,1818,ioerr, 2);  endif
# 1819 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (allocated(iEqOrbitals)) then;    deallocate(iEqOrbitals, stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas01/&
&chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&,1819,ioerr, 2);  endif
# 1820 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (allocated(iEqBlockDFTBU)) then;    deallocate(iEqBlockDFTBU, stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasa&
&s01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&,1820,ioerr, 2);  endif
# 1821 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (allocated(iEqBlockDFTBULS)) then;    deallocate(iEqBlockDFTBULS, stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/pa&
&nasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&,1821,ioerr, 2);  endif
# 1822 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"

    call destroy(neighborList)
    call destroy(pChrgMixer)
    if (tCoordOpt) then
      call destroy(pGeoCoordOpt)
if (allocated(tmpCoords)) then;    deallocate(tmpCoords, stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/panasas01/chem&
&/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&,1827,ioerr, 2);  endif
# 1828 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
    end if
    if (tLatOpt) then
      call destroy(pGeoLatOpt)
    end if
    call destroy(pRanlux)
    if (tMD) then
      call destroy(pMDFrame)
      call destroy(pMDIntegrator)
    end if
    if (tDerivs) then
      call destroy(pDerivDriver)
    end if

    if (tStoreEigvecs) then
      if (tRealHS.and.(.not.t2Component)) then
        do iS = 1, nSpin
          call destruct(storeEigvecsReal(iS))
        end do
      else
        if (t2Component) then
          call destruct(storeEigvecsCplx(1))
        else
          do iS = 1, nSpin
            call destruct(storeEigvecsCplx(iS))
          end do
        end if
      end if
    end if
if (allocated(storeEigvecsReal)) then;    deallocate(storeEigvecsReal, stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/&
&panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&,1856,ioerr, 2);  endif
# 1857 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"
if (allocated(storeEigvecsCplx)) then;    deallocate(storeEigvecsCplx, stat=ioerr);    if (ioerr /= 0) call allocateError("/panfs/&
&panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"&
&,1857,ioerr, 2);  endif
# 1858 "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/prg_dftb/initprogram.F90"

    if (tProjEigenvecs) then
      call destroy(iOrbRegion)
      call destroy(RegionLabels)
    end if

  end subroutine destroyProgramVariables

end module initprogram
