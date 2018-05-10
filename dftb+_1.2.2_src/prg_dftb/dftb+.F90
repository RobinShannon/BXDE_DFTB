!!* The main dftb+ program
program dftb
# include "allocate.h"
  use constants
  use initprogram
  use inputdata_
  use nonscc
  use Eigenvects
  use repulsive
  use eTemp
  use populations
  use densityMatrix
  use forces
  use stress
  use LapackRoutines, only : matinv ! reguired to calculate lattice derivs
  ! from stress tensors
  use SimpleAlgebra, only : determinant33 ! required for cell volumes
  use TaggedOutput
  use SCC
  use ExternalCharges
  use Mixer
  use GeoOpt
  use numDerivs2
  use Spin
  use DFTBplsU
  use FileId
  use FormatOut
  use MDCommon
  use Connectivity
  use Energies
  use Potentials
  use OrbitalEquiv
  use Parser
  use flib_wxml
  use sparse2Dense, only : unpackHS, packHS, packERho, iPackHS
  use blasRoutines, only : hemv
  use HSDUtils
  use Charmanip
  use shift
  use spinorbit
  use angmomentum
  use elecConstraints
!  use EMFields
  use MainIO

  implicit none

  !! Revision control strings
  character(len=*), parameter :: revision = &
      &"$Revision: 4298 $"
  character (len=*), parameter :: headURL = &
      &"$HeadURL: svn+ssh://svn/noodle/tags/_release/_1.2.2/src/prg_dftb/dftb+.F90 $"

  type(inputData)          :: input             ! Contains the parsed input

  integer                  :: nk, iEgy, nSpin2, nK2, iSpin2, iK2
  complex(dp), allocatable :: HSqrCplx(:,:,:,:), SSqrCplx(:,:)
  real(dp),    allocatable :: HSqrReal(:,:,:), SSqrReal(:,:)
  real(dp),    allocatable :: eigen(:,:,:)
  real(dp), allocatable    :: rhoPrim(:,:)
  real(dp), allocatable    :: iRhoPrim(:,:)
  real(dp), allocatable    :: ERhoPrim(:)
  real(dp), allocatable    :: h0(:)

  ! variables for derivatives using the Hellmann-Feynman theorem:
  real(dp), allocatable    :: hprime(:,:) ! for derivatives of H wrt external
  real(dp), allocatable    :: potentialDerivative(:,:) ! for derivatives of V
  real(dp), allocatable    :: dipoleTmp(:,:) ! temporary dipole data
  
  real(dp), allocatable    :: filling(:,:,:)
  real(dp)                 :: Eband(2), TS(2), Ef(2), E0(2), Eold

  type(TEnergies)          :: energy
  type(TPotentials)        :: potential


  real(dp), allocatable    :: derivs(:,:),repulsiveDerivs(:,:),totalDeriv(:,:)
  real(dp), allocatable    :: chrgForces(:,:)
  real(dp) :: elecStress(3,3), repulsiveStress(3,3), virialStress(3,3)
  real(dp) :: dispStress(3,3), cellVolStress(3,3), totalStress(3,3)
  real(dp) :: elecLatDeriv(3,3), repulsiveLatDeriv(3,3), virialLatDeriv(3,3)
  real(dp) :: dispLatDeriv(3,3), totalLatDeriv(3,3)
  real(dp) :: dipoleMoment(3)
  real(dp) :: angularMomentum(3) ! hold total angular momentum vector

  integer                  :: ii, jj, kk

  logical                  :: tConverged

  logical, parameter       :: tDensOn = .false.  ! O(N^2) density mtx creation
  logical, parameter       :: tAppendDetailedOut = .false.

  character(len=*), parameter :: formatEnergy = '(8f12.5)'
  character(len=*), parameter :: formatEigen = '(8f14.8)'
  character(len=*), parameter :: formatHessian = '(4f16.10)'
  character(len=*), parameter :: formatGeoOut = "(I5,F16.8,F16.8,F16.8)"
  character(len=*), parameter :: bxdFormat = "(A10,F16.8,I10)"
  character(len=*), parameter :: bxdFormat3 = "(A8,F16.8,F16.8,F16.8)"
  character(len=*), parameter :: bxdFormat4 = "(A8,F16.8,F16.8,I10)"
  character(len=*), parameter :: bxdFormat2 = "(A8,I10)"
  character(len=*), parameter :: bxdFormat5 = "(A8,F16.10)"
  ! formats for data with 1 or two units, and exponential notation form:
  character(len=*), parameter :: format1U = "(' ',A,':',T32,F18.10,T51,A)"
  character(len=*), parameter :: format2U = &
      &"(' ',A,':',T32,F18.10,T51,A,T54,F16.4,T71,A)"
  character(len=*), parameter :: format1Ue = "(' ',A,':',T37,E13.6,T51,A)"
  character(len=*), parameter :: format2Ue = &
      &"(' ',A,':',T37,E13.6,T51,A,T57,E13.6,T71,A)"
  character(len=*), parameter :: format1U1e = &
      &"(' ',A,':',T32,F18.10,T51,A,T57,E13.6,T71,A)"
  real(dp) :: cellVol, cellPressure

  !! Variables for the geometry optimization
  integer :: iGeoStep                      !* Geometry steps so far
  integer :: consitantWindow = 0
  logical :: end = .false.
  logical :: RealCOM = .false.
  integer :: iLatGeoStep                   !* Lattice geometry steps so far
  logical :: tGeomEnd                      !* Do we have the final geometry?
  logical :: tCoordEnd                     !* Has this completed?
  logical :: tCoordStep                    !* do we take an optimization step
  !* on the lattice or the internal coordinates if optimizing both in a
  !* periodic geometry
  real(dp) :: invLatVec(3,3)
  real(dp) :: derivCellVol(3,3)            !* derivative of cell volume wrt to
  ! lattice vectors
  real(dp), allocatable, target :: coord0Fold(:,:) !* Folded coords (3, nAtom)
  real(dp), pointer :: pCoord0Out(:,:)  ! Coordinates to print out
  real(dp), allocatable :: new3Coord(:,:)     !* New coordinates returned by
  !* the MD routines
  real(dp) :: tmpLatVecs(9), newLatVecs(9) !* lattice vectors returned by
  ! the optimizer
  real(dp) :: tmpLat3Vecs(3,3)
  real(dp), allocatable :: velocities(:,:) !* MD velocities
  real(dp), allocatable :: movedVelo(:,:)  !* MD velocities for moved atoms
  real(dp), allocatable :: movedAccel(:,:) !* MD acceleration for moved atoms
  real(dp), allocatable :: movedMass(:,:)  !* Mass of the moved atoms
  real(dp) :: KE                           !* MD Kinetic energy
  real(dp) :: kT                           !* MD instantaneous thermal energy

  !* COM parameters
  real(dp) :: distance
  real(dp) :: oldDistance = 15000_dp
  real(dp) :: COMup
  real(dp), allocatable :: constraint(:,:)
  real(dp), allocatable :: COM1(:)
  real(dp), allocatable :: COM2(:)
  logical :: COMreverse = .false.

  real(dp) :: highE = -500000.0_dp
  integer :: adaptiveCount = 0
  

  real(dp) :: Efield(3), absEfield !* external electric field
  
  real(dp) :: diffGeo                      !* Difference between last calculated
  !* and new geometry.

  !!* Loop variables
  integer :: iSCCIter, iSpin, iAtom, iNeigh
  integer :: fdTagged  !!* File descriptor for the tagged writer
  integer :: fdUser    !!* File descriptor for the human readable output
  integer :: fdBand    !!* File descriptor for the band structure output
  integer :: fdEigvec  !!* File descriptor for the eigenvector output
  integer :: fdResultsTag !!* File descriptor for detailed.tag
  integer :: fdMD      !!* File descriptor for extra MD output
  integer :: fdDOS
  integer :: fdDOS2

  !!* Name of the human readable file
  character(*), parameter :: taggedOut = "autotest.tag"
  character(*), parameter :: userOut = "detailed.out"
  character(*), parameter :: bandOut = "band.out"
  character(*), parameter :: mdOut = "md.out"
  character(*), parameter :: resultsTag = "results.tag"

  !!* New files
  character(*), parameter :: bxdOut = "bxd.out"
  character(*), parameter :: conformOut = "conform.out"
  character(*), parameter :: transOut = "transition.xyz"
  character(*), parameter :: newMinOut = "newMin.xyz"
  character(*), parameter :: DOSout = "DOS.out"
  character(*), parameter :: DOSout2 = "DOS2.out"
  integer :: fdBXD      !!* File descriptor for BXD output
  integer :: fdConform      !!* File descriptor for Confortmational output
  integer :: fdTrans      !!* File descriptor for transition state output
  integer :: fdNewMin      !!* File descriptor for new minima output

  real(dp), allocatable :: totalOldDeriv(:,:)  !* Energy derivatives at previous step
  real(dp) :: bxdLow
  real(dp) :: bxdUp
  integer :: bxdCount  
  integer :: rxnCountDown 
  integer :: boxCount 
  integer, allocatable :: revIdx(:) 
  logical :: bxdBounded = .false.
  logical :: boxed = .false. 
  logical :: transition = .false.
  logical :: dissociation = .false.
  logical :: reversal = .false.
  logical :: bxdStarted = .false.

  !!* DOS related parameters
  logical :: DOSreverse = .false.
  integer :: DOSstuck = 0
  integer :: DOSiterations = 0
  logical :: hitLimit = .true.
  integer :: downCount = 0
  integer :: upCount = 0
  integer :: rateDown = 0
  integer :: rateUp = 0
  integer :: ratePrev = 0
  integer :: rateUpPrev = 0
  integer :: currentBox = 0
  integer, allocatable :: DOScounts(:,:)

  character(*), parameter :: hessianOut = "hessian.out"
  integer :: fdHessian !!* File descriptor

  real(dp) :: sccErrorQ           !* Charge error in the last iterations
  real(dp) :: rTmp !, r3Tmp(3)
  real(dp), allocatable :: tmpDerivs(:)
  real(dp), allocatable :: tmpMatrix(:,:)
  real(dp), allocatable :: orbitalL(:,:,:), orbitalLPart(:,:,:)
  real(dp), allocatable    :: rVecTemp(:)
  character(lc) :: lcTmp

  real(dp) :: EbandTmp, EfTmp, TSTmp, E0Tmp

  character(lc) :: tmpStr  !!* temporary character variable

  real(dp), pointer :: pDynMatrix(:,:)


  logical :: tWriteRestart = .false. !* flag to write out geometries (and
  !* charge data if scc) when moving atoms about - in the case of conjugate
  !* gradient/steepest descent the geometries are written anyway
  integer :: minSCCIter                   !* Minimal number of SCC iterations

  type(xmlf_t) :: xf
  real(dp), allocatable :: bufferRealR2(:,:)
  logical :: tStopSCC, tStopDriver   ! if scf/driver should be stopped
  integer :: ang, iSh1, iSp1
  integer :: stuckCount = 0
  integer :: idos 

  real(dp), allocatable :: shift3rd(:)

  call printDFTBHeader(revision, headURL, parserVersion)
  write (*,'(/A/A/A/)') repeat("*", 80), "** Parsing and initializing", &
      &repeat("*", 80)

  !! Parse input and set the variables in the local scope according the input.
  !! These variables are defined in the initprogram module.

  call parseHSDInput(input)
  write (*,"(/A)") "Starting initialization..."
  write (*,"(A80)") repeat("-", 80)
  call initProgramVariables(input)
  write (*,"(A80)") repeat("-", 80)
  call destroy(input)
  write (*,*)

  print *, "adaptive"
  print *, adaptive
  elecStress = 0.0_dp
  repulsiveStress = 0.0_dp
  virialStress = 0.0_dp
  cellVolStress = 0.0_dp
  totalStress = 0.0_dp
  elecLatDeriv = 0.0_dp
  repulsiveLatDeriv = 0.0_dp
  virialLatDeriv = 0.0_dp
  totalLatDeriv = 0.0_dp
  derivCellVol = 0.0_dp
  if (tWriteTagged) then
    !! Initialize tagged writer and human readable output
    call initTaggedWriter()
    fdTagged = getFileId()
    open(fdTagged, file=taggedOut, position="rewind", status="replace")
    close(fdTagged)
  end if
  if (tWriteBandDat) then
    fdBand = getFileId()
  end if
  fdEigvec = getFileId()

  if (tDerivs) then
    fdHessian = getFileId()
    open(fdHessian, file=hessianOut, position="rewind", status="replace")
  end if

  ! initially empty file
  if (tWriteDetailedOut) then
    fdUser = getFileId()
    open(fdUser, file=userOut, position="rewind", status="replace")
  end if

  ! initially open to file to be empty
  if (tMD) then
    fdMD = getFileId()
    open(fdMD, file=mdOut, position="rewind", status="replace")
  end if

  ALLOCATE_(rhoPrim, (0, nSpin))
  ALLOCATE_(h0, (0))
  ALLOCATE_(iRhoPrim, (0, nSpin))

  if (tForces) then
    ALLOCATE_(ERhoPrim, (0))
  end if

  if (tForces) then
    ALLOCATE_(derivs,(3,nAtom))
    ALLOCATE_(repulsiveDerivs,(3,nAtom))
    ALLOCATE_(totalDeriv, (3,nAtom))
    if (tExtChrg) then
      ALLOCATE_(chrgForces, (3, nExtChrg))
    end if
  end if

  call create(energy,nAtom)

  call create(potential,orb,nAtom,nSpin)
  ALLOCATE_(shift3rd, (nAtom))

  if (nSpin <= 2) then
    ALLOCATE_(eigen, (nOrb, nKPoint,nSpin))
    ALLOCATE_(filling,(nOrb,nKpoint,nSpin))
  else
    ALLOCATE_(eigen, (2*nOrb, nKPoint,1))
    ALLOCATE_(filling,(2*nOrb,nKpoint,1))
  end if

  ALLOCATE_(coord0Fold, (3, nAtom))
  if (tShowFoldedCoord) then
    pCoord0Out => coord0Fold
  else
    pCoord0Out => coord0
  end if


  if (tMD.or.tDerivs) then
    ALLOCATE_(new3Coord, (3, nMovedAtom))
  end if

  if (tCoordOpt) then
    ALLOCATE_(tmpDerivs,(size(tmpCoords)))
  else
    ALLOCATE_(tmpDerivs,(0))
  end if

  if ((tMulliken .and. tSpinOrbit) .or. tImHam) then
    ALLOCATE_(orbitalL,(3,orb%mShell,nAtom))
    orbitalL = 0.0_dp
  else
    ALLOCATE_(orbitalL,(0,0,0))      
  end if

  if ((tMulliken .and. tSpinOrbit) .and. .not.  tDualSpinOrbit) then
    ALLOCATE_(orbitalLPart,(3,orb%mShell,nAtom))
    orbitalLPart = 0.0_dp
  else
    ALLOCATE_(orbitalLPart,(0,0,0))
  end if
  eigen(:,:,:) = 0.0_dp
  filling(:,:,:) = 0.0_dp
  TS(:) = 0.0_dp
  Eband(:) = 0.0_dp
  Ef(:) = 0.0_dp
  E0(:) = 0.0_dp

  if (tStoreEigvecs) then
    nSpin2 = 1
    nK2 = 1
  else
    nSpin2 = nSpin
    nK2 = nKPoint
  end if

  ! If only H/S should be printed, no allocation for square HS is needed
  if (.not. (tWriteRealHS .or. tWriteHS)) then
    if (t2Component) then
      ALLOCATE_(HSqrCplx, (2*nOrb, 2*nOrb, nK2, 1))
      ALLOCATE_(SSqrCplx, (2*nOrb, 2*nOrb))
    elseif (tRealHS) then
      ALLOCATE_(HSqrReal, (nOrb, nOrb, nSpin2))
      ALLOCATE_(SSqrReal, (nOrb, nOrb))
    else
      ALLOCATE_(HSqrCplx, (nOrb, nOrb, nK2, nSpin2))
      ALLOCATE_(SSqrCplx, (nOrb, nOrb))
    end if
  end if

  if (tMD) then
    ALLOCATE_(velocities,(3,nAtom))
    ALLOCATE_(movedVelo, (3, nMovedAtom))
    ALLOCATE_(movedAccel, (3, nMovedAtom))
    ALLOCATE_(movedMass, (3, nMovedAtom))
    ALLOCATE_(totalOldDeriv, (3, nMovedAtom))
    ALLOCATE_(revIdx, (nMovedAtom))
    ALLOCATE_(constraint, (3, nMovedAtom))
    ALLOCATE_(COM1, (3))
    ALLOCATE_(COM2, (3))
    ALLOCATE_(DOScounts, (bxdBoxes, 2))
    DOScounts = 0
    bxdCount = 0 
    rxnCountDown = 0
    boxCount = 0
    movedMass(:,:) = spread(mass(specie0(indMovedAtom)),1,3)
    velocities(:,:) = 0.0_dp
    fdBXD = getFileId()
    open(fdBXD, file=bxdOut, position="rewind", status="replace")
    fdConform = getFileId()
    open(fdConform, file=conformOut, position="rewind", status="replace")
    fdTrans = getFileId()
    open(fdTrans, file=transOut, position="rewind", status="replace")
    fdNewMin = getFileId()
    open(fdNewMin, file=newMinOut, position="rewind", status="replace")
    if (DOS .eqv. .true.) then
        fdDOS =getFileId()
        open(fdDOS, file=DOSout, position="rewind" , status="replace")
        fdDOS2 =getFileId()
        open(fdDOS2, file=DOSout2, position="rewind" , status="replace")
    end if
    call setRefBonds(pConnect,specie0,nMovedAtom,COMindex)
    call initialise(pConnect,coord0)
  end if

  !! Document initial variables in tagged form
  if (tWriteTagged) then
    open(fdTagged, file=taggedOut, position="append")
    call writeTagged(fdTagged, tag_SCC, tSCC)
    call writeTagged(fdTagged, tag_spin, nSpin)
    call writeTagged(fdTagged, tag_DFTBU, tDFTBU)
    if (tDFTBU) then
      call writeTagged(fdTagged, tag_nDFTBU, nDFTBUfunc)
    end if
    call writeTagged(fdTagged, tag_LS, tSpinOrbit)
    if (tSpinOrbit) then
      call writeTagged(fdTagged, tag_LSdual, tDualSpinOrbit)
    end if    
    call writeTagged(fdTagged, tag_dispersn, tDispersion)
    !! Only temporary hack, as long as autotest is not recalculated
    call writeTagged(fdTagged, tag_mAngSpecie, orb%nShell(:)-1)
    call writeTagged(fdTagged, tag_kPoint, kPoint)
    call writeTagged(fdTagged, tag_kWeight, kWeight)
    call writeTagged(fdTagged, tag_atomEigVal, atomEigVal)
    call writeTagged(fdTagged, tag_hubbU, hubbU)
    call writeTagged(fdTagged, tag_distribFn, iDistribFn)
    !! change this behavior to just write whole array as tag tag_nEl at the end
    !! of code changes and then reset autotests
    if (size(nEl) == 1) then
      call writeTagged(fdTagged, tag_nElUp, 0.5_dp*nEl(1))
      call writeTagged(fdTagged, tag_nElDown, 0.5_dp*nEl(1))
    else
      call writeTagged(fdTagged, tag_nElUp, nEl(1))
      call writeTagged(fdTagged, tag_nElDown, nEl(2))
    end if
    close(fdTagged)
  end if

  E0(:) = 0.0_dp


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Geometry loop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  tGeomEnd = nGeoSteps == 0

  tCoordStep = .false.
  if (tCoordOpt) then
    tCoordStep = .true.
    tCoordEnd = .false.
  end if

  iGeoStep = 0
  iLatGeoStep = 0
  tStopDriver = .false.

  lpGeomOpt: do while ((iGeoStep <= nGeoSteps) .and. (end .eqv. .false.))

    if (restartFreq > 0 .and. (tGeoOpt .or. tMD)) then
      tWriteRestart = (iGeoStep == nGeoSteps .or. &
          & (mod(iGeoStep, restartFreq) == 0 ))
    else
      tWriteRestart = .false.
    end if        
    
    if (tMD.and.tWriteRestart) then
      write(fdMD,*)"MD step:",iGeoStep
      call state(pMDIntegrator,fdMD)
    end if

    !! Write out geometry information
    write (*,'(/A)') repeat("*", 80)
    if (tCoordOpt .and. tLatOpt) then
      write (*, "('** Geometry step: ',I0,', Lattice step: ',I0)") &
          & iGeoStep,iLatGeoStep
    else
      write (*, "('** Geometry step: ',I0)") iGeoStep
    end if
    write (*,'(A/)') repeat("*", 80)
    if (tGeoOpt .or. tMD) then
      write (lcTmp, "(A,A)") trim(geoOutFile), ".gen"
      call clearFile(trim(lcTmp))
      if (tPeriodic) then
        call writeGenFormat(trim(lcTmp), pCoord0Out, specie0, specieName, &
            &latVec, tFracCoord)
      else
        call writeGenFormat(trim(lcTmp), coord0, specie0, specieName)
      end if
      write (lcTmp, "(A,A)") trim(geoOutFile), ".xyz"
      if (iGeoStep == 0) then
        call clearFile(trim(lcTmp))
      end if
    end if

    if (tPeriodic) then
      invLatVec = transpose(latVec)
      call matinv(invLatVec)
      cellVol = abs(determinant33(latVec))
      if (tStress.and.pressure/=0.0_dp) then
        call derivDeterminant33(derivCellVol,latVec)
        derivCellVol(:,:) = pressure * derivCellVol(:,:)
        cellVolStress(:,:) = -matmul(derivCellVol,transpose(latVec))/cellVol
      end if
    end if

    !! Save old coordinates and fold coords to unit cell
    coord0Fold(:,:) = coord0
    if (tPeriodic) then
      call foldCoordToUnitCell(coord0Fold, latVec, recVec2p)
    end if

    !! Initialize neighborlists
    call updateNeighborListAndSpecies(coord, specie, img2CentCell, iCellVec, &
        &neighborList, nAllAtom, coord0Fold, specie0, mCutoff, rCellVec)
    nAllOrb = sum(orb%nOrbSpecie(specie(1:nAllAtom)))

    !! Calculate neighborlist for SK and repulsive calculation
    call getNrOfNeighborsForAll(nNeighbor, neighborList, skRepCutoff)

    !! Reallocate Hamiltonian and overlap based on the new neighbor list
    call reallocateHS(ham, over, iPair, neighborList%iNeighbor, nNeighbor, &
        &orb, img2CentCell)

    !! Reallocate density matrixes if necessary
    if (size(ham, dim=1) > size(rhoPrim, dim=1)) then
      DEALLOCATE_(H0)
      ALLOCATE_(H0,(size(ham,dim=1)))
      DEALLOCATE_(rhoPrim)
      ALLOCATE_(rhoPrim,(size(ham,dim=1),nSpin))
      if (tImHam) then
        DEALLOCATE_(iRhoPrim)
        ALLOCATE_(iRhoPrim,(size(ham,dim=1),nSpin))
        DEALLOCATE_PARR(iHam)
        INITALLOCATE_PARR(iHam,(size(ham,dim=1),4))
      end if
      if (tForces) then
        DEALLOCATE_(ERhoPrim)
        ALLOCATE_(ERhoPrim,(size(ham,dim=1)))
      end if
    end if

    !! (Re)Initialize mixer
    if (tSCC) then
      call reset(pChrgMixer, nMixElements)
    end if

    if (tWriteTagged) then
      !! Write out initialization information
      if (iGeoStep == 0) then
        open(fdTagged, file=taggedOut, position="append")
        call writeTagged(fdTagged, tag_initCoord, coord(:,:nAtom))
        call writeTagged(fdTagged, tag_specie, specie(:nAtom))
        call writeTagged(fdTagged, tag_nNeighbor, neighborList%nNeighbor)
        call writeTagged(fdTagged, tag_iNeighbor, neighborList%iNeighbor)
        close(fdTagged)
      end if
    end if

    !! Notify various modules about coordinate changes
    if (tSCC) then
      call updateCoords_SCC(coord, specie, neighborList, img2CentCell)
    end if
    if (tDispersion) then
      call updateCoords(myDispersion, neighborList, img2CentCell, coord, &
          &specie0)
    end if
    if (t3rdFull) then
      call updateCoords(thirdOrd, neighborList, specie)
    end if

    !! Build non-scc Hamiltonian
    call buildH0S(H0,over, skHamCont, skOverCont, atomEigVal, coord, &
        &nNeighbor, neighborList%iNeighbor, specie, iPair, orb)

    !! Adapt electron temperature to MD, if necessary
    if (tSetFillingTemp) then
      call getTemperature(pTempProfile, tempElec)
    end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! SCC-loop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (tSCC .and. (.not. tAppendDetailedOut)) then
      if (.not. tDFTBU) then
        write(*,"(' ',A5,A18,A18,A18)") "iSCC", " Total electronic ", &
            & "  Diff electronic ", "     SCC error    "
      else
        write(*,"(' ',A5,A18,A18,A18)") "iSCC", " Total electronic ", &
            & "  Diff electronic ", "     SCC error    "
      end if
    end if

    tConverged = .false.
    if (tSCC) then
      if (tDFTBU) then
        minSCCIter = 2
      else
        if (nSpin == 1) then
          minSCCIter = 1
        else
          minSCCIter = 2
        end if
      end if
    else
      minSCCIter = 1
    end if

    energy%ETotal = 0.0_dp
    energy%atomTotal(:) = 0.0_dp

    !! Calculate repulsive energy
    call getERep(energy%atomRep, coord, nNeighbor, neighborList%iNeighbor, &
        &specie, pRepCont, img2CentCell, tUseBuggyRepSum)
    energy%Erep = sum(energy%atomRep)

    if (tDispersion) then
      call getEnergies(myDispersion, energy%atomDisp)
      energy%eDisp = sum(energy%atomDisp)
    else
      energy%atomDisp(:) = 0.0_dp
    end if
    
    potential%extAtom = 0.0_dp
    potential%extShell = 0.0_dp
    potential%extBlock = 0.0_dp
    
    if (tSCC) then
      if (tEField) then
        Efield(:) = EFieldStrength * EfieldVector(:)         
        if (tTDEfield) then
          Efield(:) = Efield(:) &
              & * sin(EfieldOmega*deltaT*real(iGeoStep+EfieldPhase,dp))
        end if
        absEfield = sqrt(sum(Efield**2))
        if (tPeriodic) then
          do iAtom = 1, nAtom
            do iNeigh = 1, nNeighbor(iAtom)
              ii = neighborList%iNeighbor(iNeigh,iAtom)
              if (iCellVec(ii) /= 0) then ! overlap between atom in central
                !  cell and non-central cell
                if (abs(dot_product(cellVec(:,iCellVec(ii)),EfieldVector))&
                    & /= 0.0_dp) then ! component of electric field projects
                  ! onto vector between cells
                  write(tmpStr,"('Interaction between atoms ',I0,' and ', I0,&
                      &' crosses the saw-tooth discontinuity in the &
                      &electric field.')")iAtom,img2centcell(ii)
                  call error(tmpStr)
                end if
              end if
            end do
          end do
          do iAtom = 1, nAtom
            potential%extAtom(iAtom,1)=dot_product(coord0Fold(:,iAtom),Efield)
          end do
        else
          do iAtom = 1, nAtom
            potential%extAtom(iAtom,1)=dot_product(coord(:,iAtom),Efield)
          end do
        end if
      else
        Efield = 0.0_dp
      end if
            
      call total_shift(potential%extShell,potential%extAtom, &
          & potential%extShell,orb,specie)        
      call total_shift(potential%extBlock,potential%extShell, &
          & potential%extBlock,orb,specie)
    end if
    
    iSCCIter = 1
    tStopSCC = .false.
    lpSCC: do while (iSCCiter <= nSCCIter)

      rhoPrim(:,:) = 0.0_dp
      if (tImHam) then
        iRhoPrim(:,:) = 0.0_dp
      end if

      ham(:,:) = 0.0_dp
      do ii = 1, size(H0)
        ham(ii,1) = h0(ii)
      end do

      !! Build various contribution to the Hamiltonian
      
      potential%iorbitalBlock = 0.0_dp
      if (tDualSpinOrbit) then
        call shiftLS(potential%iorbitalBlock,xi,orb,specie)
      end if
            
      if (.not. tSCC) then
        tConverged = .true.
        potential%intBlock = 0.0_dp
      else
        chargePerShell(:,:,:) = 0.0_dp ! hack for the moment to get charge
        ! and magnetization
        do iAtom = 1, nAtom
          iSp1 = specie(iAtom)
          do iSh1 = 1, orb%nShell(iSp1)
            chargePerShell(iSh1,iAtom,1:nSpin) = &
                & chargePerShell(iSh1,iAtom,1:nSpin) + &
                & sum(qInput(orb%posShell(iSh1,iSp1): &
                & orb%posShell(iSh1+1,iSp1)-1,iAtom,1:nSpin),dim=1)
          end do
        end do

        potential%intAtom = 0.0_dp
        potential%intShell = 0.0_dp
        potential%intBlock = 0.0_dp
        call updateCharges_SCC(qInput, q0, orb, specie, &
            &neighborList%iNeighbor, img2CentCell)
        call getShiftPerAtom(potential%intAtom)
        call getShiftPerL(potential%intShell)
        
        if (t3rdFull) then
          call updateCharges(thirdOrd, specie0, neighborList, qInput, q0, &
              &img2CentCell)
          call getShiftAtom(thirdOrd, shift3rd)
          potential%intAtom(:,1) = potential%intAtom(:,1) + shift3rd
        end if
        
        call total_shift(potential%intShell,potential%intAtom, &
            & potential%intShell,orb,specie)
        
        !! Build spin contribution (if necessary)
        if (tSpin) then
          call addSpinShift(potential%intShell,chargePerShell,specie,orb,W)
        end if

        call total_shift(potential%intBlock,potential%intShell, &
            & potential%intBlock,orb,specie)

        if (tDFTBU) then !! Apply LDA+U correction (if necessary)
          potential%orbitalBlock = 0.0_dp
          if (tImHam) then
            call shift_DFTBU(potential%orbitalBlock,potential%iorbitalBlock, &
                & qBlockIn, qiBlockIn, specie,orb, nDFTBUfunc, &
                & UJ, nUJ, niUJ, iUJ)
          else
            call shift_DFTBU(potential%orbitalBlock,qBlockIn,specie,orb, &
                & nDFTBUfunc, UJ, nUJ, niUJ, iUJ)            
          end if
          potential%intBlock = potential%intBlock + potential%orbitalBlock
        end if

      end if
      
      potential%intBlock = potential%intBlock + potential%extBlock
      
      call add_shift(ham,over,nNeighbor, neighborList%iNeighbor, &
          & specie,orb,iPair,nAtom,img2CentCell,potential%intBlock)

      if (tImHam) then
        iHam = 0.0_dp
        call add_shift(iHam,over,nNeighbor, neighborList%iNeighbor, &
            & specie,orb,iPair,nAtom,img2CentCell,potential%iorbitalBlock)
        iHam(:,:) = 2.0_dp*iHam(:,:)              
      end if

      ! hack due to not using Pauli-type structure for diagonalisation
      ! etc.
      if (nSpin>1) then
        ham(:,:) = 2.0_dp*ham(:,:)
      end if

      if (nSpin /= 4) then

        if (nSpin == 2) then
          call qm2ud(ham)
        end if

        ! Write out matrices if necessary and quit.
        call writeHS(tWriteHS, tWriteRealHS, ham, over, neighborList%iNeighbor, &
            &nNeighbor, iAtomStart, iPair, img2CentCell, kPoint, iCellVec, &
            &cellVec)
        if (tWriteRealHS .or. tWriteHS) then
          write (*, "(A)") "Hamilton/Overlap written, exiting program."
          stop
        end if

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Spin loop
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        lpSpin: do iSpin = 1, nSpin


          if (tStoreEigvecs) then
            iSpin2 = 1
            if (tRealHS) then
              call reset(storeEigvecsReal(iSpin), (/ nOrb, nOrb /))
            else
              call reset(storeEigvecsCplx(iSpin), (/ nOrb, nOrb /))
            end if
          else
            iSpin2 = iSpin
          end if


          !! Solve eigenproblem for real or complex Hamiltionian
          if (tRealHS) then
            call diagonalize(HSqrReal(:,:,iSpin2), SSqrReal, eigen(:,1,iSpin), &
                &ham(:,iSpin), over, neighborList%iNeighbor, nNeighbor, &
                &iAtomStart, iPair, img2CentCell, solver, 'V')

            if (tStoreEigvecs) then
              call push(storeEigvecsReal(iSpin), HSqrReal(:,:,iSpin2))
            end if

            !! Occupy eigenlevels including possible spin polarisation
            call Efilling(Eband(iSpin),Ef(iSpin),TS(iSpin),&
                & E0(iSpin),&
                & filling(:,:,iSpin:iSpin),&
                & eigen(:,:,iSpin:iSpin),&
                & 0.5_dp*real(nSpin,dp)*nEl(iSpin), tempElec,kWeight,iDistribFn)

            !! prefactor 2.0_dp for spin averaged, 1.0_dp for spin polarised
            Eband(iSpin) = real(3-nSpin,dp)*Eband(iSpin)
            E0(iSpin) = real(3-nSpin,dp)*E0(iSpin)
            TS(iSpin) = real(3-nSpin,dp)*TS(iSpin)
            filling(:,:,iSpin:iSpin) = real(3-nSpin,dp)* &
                & filling(:,:,iSpin:iSpin)

            if (tDensON) then
              call makeDensityMatrix(SSqrReal, HSqrReal(:,:,iSpin2), &
                  &filling(:,1,iSpin), neighborlist%iNeighbor, nNeighbor, orb, &
                  &iAtomStart, img2CentCell)
            else
              call makeDensityMatrix(SSqrReal, HSqrReal(:,:,iSpin2), &
                  &filling(:,1,iSpin))
            end if
            call packHS(rhoPrim(:,iSpin), SSqrReal, neighborlist%iNeighbor, &
                &nNeighbor, orb%mOrb, iAtomStart, iPair, img2CentCell)
          else
            !! k-point calculation
            !! Do unpacking/diagonalizing/density mtx creating in each K-point

            nkLoop: do nk = 1, nKPoint
              if (tStoreEigvecs) then
                iK2 = 1
              else
                iK2 = nK
              end if
              call diagonalize(HSqrCplx(:,:,iK2,iSpin2), SSqrCplx(:,:), &
                  &eigen(:,nk,iSpin), ham(:,iSpin), over,kPoint(:,nk), &
                  &neighborList%iNeighbor, nNeighbor, iCellVec, cellVec, &
                  &iAtomStart, iPair, img2CentCell, solver, "V")
              if (tStoreEigvecs) then
                call push(storeEigvecsCplx(iSpin), HSqrCplx(:,:,iK2,iSpin2))
              end if

            end do nkLoop

            !! Fill up levels with electrons
            if (.not. tFillKSep) then
              call Efilling(Eband(iSpin),Ef(iSpin),TS(iSpin),&
                  & E0(iSpin), &
                  & filling(:,:,iSpin:iSpin),&
                  & eigen(:,:,iSpin:iSpin),&
                  & 0.5_dp*real(nSpin,dp)*nEl(iSpin),tempElec,kWeight, &
                  & iDistribFn)
            else
              !! Make aufbau principle in each K-point separately
              do nk = 1, nKPoint
                call Efilling(EbandTmp,EfTmp,TSTmp,&
                    & E0Tmp, &
                    & filling(:,nk:nk,iSpin:iSpin),&
                    & eigen(:,nk:nk,iSpin:iSpin),&
                    & 0.5_dp*real(nSpin,dp)*nEl(iSpin),tempElec,kWeight(nk:nk),&
                    & iDistribFn)
                Eband(iSpin) = eBand(iSpin) + eBandTmp
                Ef(iSpin) = Ef(iSpin) + EfTmp
                TS(iSpin) = TS(iSpin) + TSTmp
                E0(iSpin) = E0(iSpin) + E0Tmp
              end do
              !! Fermi energy meaningless in that case -> put average Ef
              Ef(iSpin) = Ef(iSpin) / real(nKPoint, dp)
            end if

            ! prefactor 2.0_dp for spin averaged, 1.0_dp for spin polarised
            Eband(iSpin) = real(3-nSpin,dp)*Eband(iSpin)
            E0(iSpin) = real(3-nSpin,dp)*E0(iSpin)
            TS(iSpin) = real(3-nSpin,dp)*TS(iSpin)
            filling(:,:,iSpin:iSpin) = real(3-nSpin,dp)* &
                & filling(:,:,iSpin:iSpin)

            nkLoop2: do nk = 1, nKPoint

              if (tStoreEigvecs) then
                iK2 = 1
                call get(storeEigvecsCplx(iSpin), HSqrCplx(:,:,iK2,iSpin2))
              else
                iK2 = nK
              end if

              if (tDensON) then
                call makeDensityMatrix(SSqrCplx, HSqrCplx(:,:,iK2,iSpin2), &
                    &filling(:,nK,iSpin), neighborlist%iNeighbor, nNeighbor, &
                    & orb, iAtomStart, img2CentCell)
              else
                call makeDensityMatrix(SSqrCplx, HSqrCplx(:,:,iK2,iSpin2), &
                    &filling(:,nK,iSpin))
              end if
              call packHS(rhoPrim(:,iSpin), SSqrCplx, kPoint(:,nK), &
                  &kWeight(nk), neighborList%iNeighbor, nNeighbor, orb%mOrb, &
                  &iCellVec, cellVec, iAtomStart, iPair, img2CentCell)

            end do nkLoop2
          end if

          if (tWriteBandDat) then
            if (iSpin == 1) then
              open(unit=fdBand, file=bandOut, action="write")
            end if
            do nk=1,nKPoint
              if (nSpin <= 2) then
                write(fdBand,*)'KPT ',nk,' SPIN ', iSpin, &
                    &' KWEIGHT ', kweight(nK)
              else
                write(fdBand,*)'KPT ', nk, ' KWEIGHT ', kweight(nK)
              end if
              do iEgy=1,nOrb
                write(fdBand, formatEnergy) Hartree__eV*eigen(iEgy,nk,iSpin),&
                    & filling(iEgy,nk,iSpin)
              end do
              write(fdBand,*)
            end do
            if (iSpin == nSpin) then
              close(fdBand)
            end if
          end if


        end do lpSpin

        call ud2qm(rhoPrim)

      else

        if (tStoreEigvecs) then
          call reset(storeEigvecsCplx(1), (/ 2*nOrb, 2*nOrb /))
        end if
        if (tRealHS) then
          if (tImHam) then
            call diagonalize(HSqrCplx(:,:,1,1), SSqrCplx, eigen(:,1,1), &
                &ham(:,:), over, neighborList%iNeighbor, nNeighbor, &
                &iAtomStart, iPair, img2CentCell, solver, 'V', iHam=iHam)
          else if (tSpinOrbit) then            
            call diagonalize(HSqrCplx(:,:,1,1), SSqrCplx, eigen(:,1,1), &
                &ham(:,:), over, neighborList%iNeighbor, nNeighbor, &
                &iAtomStart, iPair, img2CentCell, solver, 'V', xi,orb,specie)
          else
            call diagonalize(HSqrCplx(:,:,1,1), SSqrCplx, eigen(:,1,1), &
                &ham(:,:), over, neighborList%iNeighbor, nNeighbor, &
                &iAtomStart, iPair, img2CentCell, solver, 'V')
          end if
          if (tStoreEigvecs) then
            call push(storeEigvecsCplx(1), HSqrCplx(:,:,1,1))
          end if
        else
          nkLoop3: do nk = 1, nKPoint
            if (tStoreEigvecs) then
              iK2 = 1
            else
              iK2 = nK
            end if
            if (tImHam) then
              call diagonalize(HSqrCplx(:,:,iK2,1), SSqrCplx, eigen(:,nk,1), &
                  & ham(:,:), over, kPoint(:,nk),neighborList%iNeighbor, &
                  & nNeighbor, iCellVec, cellVec,iAtomStart, iPair, &
                  & img2CentCell, solver, 'V', iHam=iHam)
            else if (tSpinOrbit) then             
              call diagonalize(HSqrCplx(:,:,iK2,1), SSqrCplx, eigen(:,nk,1), &
                  & ham(:,:), over, kPoint(:,nk),neighborList%iNeighbor, &
                  & nNeighbor, iCellVec, cellVec,iAtomStart, iPair, &
                  & img2CentCell, solver, 'V', xi,orb,specie)
            else
              call diagonalize(HSqrCplx(:,:,iK2,1), SSqrCplx, eigen(:,nk,1), &
                  & ham(:,:), over, kPoint(:,nk),neighborList%iNeighbor, &
                  & nNeighbor, iCellVec, cellVec,iAtomStart, iPair, &
                  & img2CentCell, solver, 'V')
            end if
            if (tStoreEigvecs) then
              call push(storeEigvecsCplx(1), HSqrCplx(:,:,iK2,1))
            end if
          end do nkLoop3
        end if
        
        !! Fill up levels with electrons
        if (.not. tFillKSep) then
          call Efilling(Eband(1),Ef(1),TS(1),E0(1), filling(:,:,1:1),&
              & eigen(:,:,1:1),nEl(1),tempElec,kWeight, iDistribFn)
        else
          !! Make aufbau principle in each K-point separately
          do nk = 1, nKPoint
            call Efilling(EbandTmp,EfTmp,TSTmp,E0Tmp, filling(:,nk:nk,1:1),&
                & eigen(:,nk:nk,1:1),nEl(1),tempElec,kWeight(nk:nk), &
                & iDistribFn)
            Eband(1) = eBand(1) + eBandTmp
            Ef(1) = Ef(1) + EfTmp
            TS(1) = TS(1) + TSTmp
            E0(1) = E0(1) + E0Tmp
          end do
          !! Fermi energy meaningless in that case -> put average Ef
          Ef(1) = Ef(1) / real(nKPoint, dp)
        end if

        filling(:,1:nKPoint,1) = 2.0_dp * filling(:,1:nKPoint,1)
        SSqrCplx = 0.0_dp
        if (tSpinOrbit) then
          energy%atomLS = 0.0_dp
          if (.not.tDualSpinOrbit) then
            ALLOCATE_(rVecTemp,(nAtom))
          end if
        end if
        if (tImHam .and. tMulliken) then
          orbitalL = 0.0_dp
        end if
        nkLoop4: do nk = 1, nKPoint
          if (tStoreEigvecs) then
            iK2 = 1
            call get(storeEigvecsCplx(1), HSqrCplx(:,:,iK2,1))
          else
            iK2 = nK
          end if
          if (tDensON) then
            call error("Currently missing.")
          else
            call makeDensityMatrix(SSqrCplx, HSqrCplx(:,:,iK2,1), &
                &filling(:,nK,1))
          end if
          
          if (tSpinOrbit .and. .not. tDualSpinOrbit) then
            rVecTemp = 0.0_dp
            call getEnergySpinOrbit(rVecTemp,SSqrCplx,iAtomStart, &
                & xi, orb, specie)
            energy%atomLS = energy%atomLS + kWeight(nk)*rVecTemp
            if (tMulliken) then
              orbitalLPart = 0.0_dp
              call getL(orbitalLPart,SSqrCplx,iAtomStart, orb, specie)
              orbitalL = orbitalL + kWeight(nk) * orbitalLPart
            end if
          end if

          if (tRealHS) then
            call packHS(rhoPrim(:,:), SSqrCplx, neighborlist%iNeighbor, &
                &nNeighbor, orb%mOrb, iAtomStart, iPair, img2CentCell)
            if (tImHam) then              
              call iPackHS(iRhoPrim(:,:), SSqrCplx, neighborlist%iNeighbor, &
                  &nNeighbor, orb%mOrb, iAtomStart, iPair, img2CentCell)
            end if
          else
            call packHS(rhoPrim(:,:), SSqrCplx, kPoint(:,nK), &
                &kWeight(nk), neighborList%iNeighbor, nNeighbor, orb%mOrb, &
                &iCellVec, cellVec, iAtomStart, iPair, img2CentCell)
            if (tImHam) then              
              call iPackHS(iRhoPrim(:,:), SSqrCplx, kPoint(:,nK), &
                  & kWeight(nk), neighborlist%iNeighbor, &
                  & nNeighbor, orb%mOrb, iCellVec, cellVec, iAtomStart, &
                  & iPair, img2CentCell)
            end if
          end if

        end do nkLoop4

        if (tSpinOrbit .and. .not. tDualSpinOrbit) then
          DEALLOCATE_(rVecTemp)
          energy%ELS = sum(energy%atomLS(:))
        end if
        filling(:,1:nKPoint,1) = 0.5_dp * filling(:,1:nKPoint,1)


        if (tWriteBandDat) then
          open(unit=fdBand, file=bandOut, action="write")
          do nk=1,nKPoint
            if (nSpin <= 2) then
              write(fdBand,*)'KPT ',nk,' SPIN ', iSpin, &
                  &' KWEIGHT ', kweight(nk)
            else
              write(fdBand,*)'KPT ',nk, ' KWEIGHT ', kweight(nk)
            end if
            do iEgy=1,2*nOrb
              write(fdBand, formatEnergy) Hartree__eV*eigen(iEgy,nk,1),&
                  & filling(iEgy,nk,1)
            end do
            write(fdBand,*)
          end do
          close(fdBand)
        end if

      end if

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Mulliken analysis
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      if (tMulliken) then
        qOutput(:,:,:) = 0.0_dp
        do iSpin = 1, nSpin
          call mulliken(qOutput(:,:,iSpin), over, rhoPrim(:,iSpin), &
              &orb, neighborList%iNeighbor, nNeighbor, img2CentCell, iPair)
        end do
      end if

      if (tImHam) then
        qiBlockOut(:,:,:,:) = 0.0_dp
        energy%atomLS = 0.0_dp
        do iSpin = 1, nSpin
          call skewMulliken(qiBlockOut(:,:,:,iSpin), over, iRhoPrim(:,iSpin), &
              &orb, neighborList%iNeighbor, nNeighbor, img2CentCell, iPair)
        end do
        call getL(orbitalL,qiBlockOut,orb,specie)
        if (tDualSpinOrbit) then
          call getEnergySpinOrbit(energy%atomLS,qiBlockOut,xi,orb,specie)
          energy%ELS = sum(energy%atomLS(:))
        end if
        qBlockOut(:,:,:,:) = 0.0_dp
      end if

      if (tDFTBU) then
        qBlockOut(:,:,:,:) = 0.0_dp
        do iSpin = 1, nSpin
          call mulliken(qBlockOut(:,:,:,iSpin), over, rhoPrim(:,iSpin), &
              &orb, neighborList%iNeighbor, nNeighbor, img2CentCell, iPair)
        end do
      end if

      if (tSCC) then

        chargePerShell(:,:,:) = 0.0_dp ! hack for the moment to get charge
        ! and magnetization
        do iAtom = 1, nAtom
          iSp1 = specie(iAtom)
          do iSh1 = 1, orb%nShell(iSp1)
            chargePerShell(iSh1,iAtom,1:nSpin) = &
                & chargePerShell(iSh1,iAtom,1:nSpin) + &
                & sum(qOutput(orb%posShell(iSh1,iSp1): &
                & orb%posShell(iSh1+1,iSp1)-1,iAtom,1:nSpin),dim=1)
          end do
        end do

        ! recalculate the SCC shifts for the output charge.

        !! SCC contribution is calculated with the output charges.
        call updateCharges_SCC(qOutput, q0, orb, specie, &
            &neighborList%iNeighbor, img2CentCell)

        potential%intAtom = 0.0_dp
        potential%intShell = 0.0_dp
        potential%intBlock = 0.0_dp

        call getShiftPerAtom(potential%intAtom)
        call getShiftPerL(potential%intShell)
        if (t3rdFull) then
          call updateCharges(thirdOrd, specie0, neighborList, qOutput, q0, &
              &img2CentCell)
          call getShiftAtom(thirdOrd, shift3rd)
          potential%intAtom(:,1) = potential%intAtom(:,1) + shift3rd
        end if

        call total_shift(potential%intShell,potential%intAtom, &
            & potential%intShell,orb,specie)

        !! Build spin contribution (if necessary)
        if (tSpin) then
          call addSpinShift(potential%intShell,chargePerShell,specie,orb,W)
        end if

        call total_shift(potential%intBlock,potential%intShell, &
            & potential%intBlock,orb,specie)
      end if
      
      !! Calculate energies
      
      ! non-SCC part
      energy%EnonSCC = 0.0_dp
      energy%atomNonSCC(:) = 0.0_dp
      
      call mulliken(energy%atomNonSCC(:), rhoPrim(:,1), H0,orb,&
          &neighborList%iNeighbor, nNeighbor, img2CentCell, iPair)
      energy%EnonSCC =  sum(energy%atomNonSCC)
      
      if (tEfield) then ! energy in external field
        energy%atomExt = -sum( q0(:, :, 1) - qOutput(:, :, 1),dim=1) &
            & * potential%extAtom(:,1)
        energy%Eext =  sum(energy%atomExt)
      else
        energy%Eext = 0.0_dp
        energy%atomExt = 0.0_dp
      end if
      
      if (tSCC) then
        call getEnergyPerAtom_SCC(energy%atomSCC, specie, &
            &neighborList%iNeighbor, img2CentCell)
        energy%eSCC = sum(energy%atomSCC)
        if (t3rdFull) then
          call getEnergyPerAtom(thirdOrd, energy%atom3rd)
          energy%e3rd = sum(energy%atom3rd)
        end if

        if (nSpin > 1) then
          energy%atomSpin = 0.5_dp * sum(sum(potential%intShell(:,:,2:nSpin)&
              & * chargePerShell(:,:,2:nSpin), dim=1),dim=2)
          energy%Espin = sum(energy%atomSpin)
        else
          energy%atomSpin = 0.0_dp
          energy%eSpin = 0.0_dp
        end if
      end if

      potential%iorbitalBlock = 0.0_dp
      if (tDualSpinOrbit) then
        call shiftLS(potential%iorbitalBlock,xi,orb,specie)
      end if
      
      if (tDFTBU) then
        energy%atomDftbu(:) = 0.0_dp
        if (.not. tImHam) then
          call E_DFTBU(energy%atomDftbu,qBlockOut,specie,orb, &
              & nDFTBUfunc, UJ, nUJ, niUJ, iUJ)
        else
          call E_DFTBU(energy%atomDftbu,qBlockOut,specie,orb, &
              & nDFTBUfunc, UJ, nUJ, niUJ, iUJ, qiBlockOut)
        end if
        energy%Edftbu = sum(energy%atomDftbu(:))
        potential%orbitalBlock = 0.0_dp
        
        if (tImHam) then  
          call shift_DFTBU(potential%orbitalBlock,potential%iorbitalBlock, &
              & qBlockOut,qiBlockOut, specie,orb, nDFTBUfunc, UJ, nUJ, niUJ, &
              & iUJ)
        else
          call shift_DFTBU(potential%orbitalBlock,qBlockOut,specie,orb, &
              & nDFTBUfunc, UJ, nUJ, niUJ, iUJ)
        end if
      else
        energy%Edftbu = 0.0_dp
      end if
      
      energy%Eelec = energy%EnonSCC + energy%ESCC + energy%Espin &
          & + energy%ELS + energy%Edftbu + energy%Eext + energy%e3rd 

      energy%atomElec(:) = energy%atomNonSCC(:) &
          & + energy%atomSCC(:) + energy%atomSpin(:) + energy%atomDftbu(:) &
          & + energy%atomLS(:) + energy%atomExt(:) + energy%atom3rd(:)

      energy%atomTotal(:) = energy%atomElec(:) + energy%atomRep(:) + &
          & energy%atomDisp(:)

      energy%Etotal = energy%Eelec + energy%Erep + energy%eDisp

      !! Stop SCC if appropriate stop file is present (We need this query here
      !! since the following block contains a check iSCCIter /= nSCCIter)
      inquire(file=fStopSCC, exist=tStopSCC)
      if (tStopSCC) then
        write (*,*) "Stop file '" // fStopSCC // "' found."
        nSCCIter = iSCCIter
        write (*,*) "Setting max number of scc cycles to current cycle."
      end if

      !! Mix charges
      if (tSCC) then
        qOutRed = 0.0_dp
        if (nSpin == 2) then
          call qm2ud(qOutput)
          if (tDFTBU) then
            call qm2ud(qBlockOut)
          end if
        end if
        call OrbitalEquiv_reduce(qOutput, iEqOrbitals, orb, &
            & qOutRed(1:nIneqOrb))
        if (tDFTBU) then          
          call AppendBlock_reduce( qBlockOut, iEqBlockDFTBU, orb, &
              & qOutRed )
          if (tImHam) then
            call AppendBlock_reduce( qiBlockOut, iEqBlockDFTBULS, orb, &
                & qOutRed, skew=.true. )
          end if
        end if
        if (nSpin == 2) then
          call ud2qm(qOutput)
          if (tDFTBU) then
            call ud2qm(qBlockOut)
          end if
        end if

        qDiffRed(:) = qOutRed(:) - qInpRed(:)        
        sccErrorQ = maxval(abs(qDiffRed))

        tConverged = (sccErrorQ < sccTol) .and. &
            & (iSCCiter >= minSCCIter .or. tReadChrg)
        if ((.not. tConverged) .and. iSCCiter /= nSCCiter) then
          !! Avoid mixing of spin unpolarised density for spin polarised
          !! cases, this is only a problem in iteration 1, as there is
          !! only the (spin unpolarised!) atomic input density at that
          !! point. (Unless charges had been initialized externally)
          if (iSCCIter == 1 .and. (nSpin > 1.or.tDFTBU) .and. .not.tReadChrg)&
              & then
            qInput(:,:,:) = qOutput(:,:,:)
            qInpRed(:) = qOutRed(:)
            if (tDFTBU) then
              qBlockIn(:,:,:,:) = qBlockOut(:,:,:,:)
              if (tSpinOrbit) then
                qiBlockIn(:,:,:,:) = qiBlockOut(:,:,:,:)
              end if
            end if
          else
            call mix(pChrgMixer, qInpRed, qDiffRed)
            call OrbitalEquiv_expand(qInpRed(1:nIneqOrb), iEqOrbitals, &
                & orb, qInput)
            if (tDFTBU) then
              qBlockIn = 0.0_dp
              call Block_expand( qInpRed ,iEqBlockDFTBU, orb, &
                  & qBlockIn, specie0, nUJ, niUJ, iUJ, orbEquiv=iEqOrbitals )
              if (tSpinOrbit) then
                call Block_expand( qInpRed ,iEqBlockDFTBULS, orb, &
                    & qiBlockIn, specie0, nUJ, niUJ, iUJ, skew=.true. )
              end if
            end if
            if (nSpin == 2) then
              call ud2qm(qInput)
              if (tDFTBU) then
                call ud2qm(qBlockIn)
              end if
            end if
          end if
        end if
        
      end if

      !! Clear detailed.out if necessary
      if (tWriteDetailedOut .and. .not. tAppendDetailedOut) then
        close(fdUser)
        open(fdUser, file=userOut, position="rewind", status="replace")
        select case(iDistribFn)
        case(0)
          write(fdUser,*)'Fermi distribution function'
        case(1)
          write(fdUser,*)'Gaussian distribution function'
        case default
          write(fdUser,*)'Methfessel-Paxton distribution function order'&
              &,iDistribFn
        end select
        write(fdUser,*) ""
      end if

      if (tSCC) then
        if (iSCCiter > 1) then
          rTmp = (energy%Eelec-Eold)
        else
          rTmp = 0.0_dp
        end if
        Eold = energy%Eelec
        
        if (tWriteDetailedOut) then
          if (nGeoSteps > 0) then
            if (tCoordOpt .and. tLatOpt) then
              write (fdUser, "('Geometry optimization step: ',I0, &
                  & ', Lattice step: ',I0)") &
                  & iGeoStep,iLatGeoStep
            else
              write(fdUser,"('Geometry optimization step: ',I0)") iGeoStep
            end if
          else
            write(fdUser,*) "Calculation with static geometry"
          end if
          write(fdUser,*)''
          write (fdUser,'(/A)') repeat("*", 80)
          if (tDFTBU) then
            write(fdUser,"(' ',A5,A18,A18,A18)") "iSCC", " Total electronic ", &
                & "  Diff electronic ", "     SCC error    "
            write(fdUser,"(I5,E18.8,E18.8,E18.8,E18.8)") iSCCIter, &
                & energy%Eelec, rTmp, sccErrorQ
          else
            write(fdUser,"(' ',A5,A18,A18,A18)") "iSCC", " Total electronic ", &
                & "  Diff electronic ", "     SCC error    "
            write(fdUser,"(I5,E18.8,E18.8,E18.8)") iSCCIter, energy%Eelec, &
                & rTmp, sccErrorQ
          end if
          write (fdUser,'(A)') repeat("*", 80)
          write(fdUser,*)""
        end if
        
        if (tDFTBU) then
          write(*,"(I5,E18.8,E18.8,E18.8)") iSCCIter, energy%Eelec, rTmp, &
              &sccErrorQ
        else
          write(*,"(I5,E18.8,E18.8,E18.8)") iSCCIter, energy%Eelec, rTmp, &
              & sccErrorQ
        end if
        
      end if

      !! Not writing any restarting info if not converged and minimal number of
      !! SCC iterations not done.
      if (restartFreq > 0 .and. .not.(tMD .or. tGeoOpt .or. tDerivs)) then
        if (tConverged .or. ((iSCCIter >= minSCCIter &
            & .or. tReadChrg .or. iGeoStep > 0) &
            &.and. (iSCCIter == nSCCIter .or. mod(iSCCIter, restartFreq) ==&
            & 0))) then
          if (tMulliken.and.tSCC) then
            if (tDFTBU) then
              if (tSpinOrbit) then
                call writeQToFile(qInput, fChargeIn, orb, qBlockIn, qiBlockIn)
              else
                call writeQToFile(qInput, fChargeIn, orb, qBlockIn)
              end if
            else
              call writeQToFile(qInput, fChargeIn, orb)
            end if
            print "('>> Charges saved for restart in ',A)", fChargeIn
          end if
        end if
      end if

      if (tWriteDetailedOut) then
        if (nMovedAtom > 0 .and. .not. tDerivs) then
          write (fdUser,*) "Coordinates of moved atoms (au):"
          do ii = 1, nMovedAtom
            write(fdUser,formatGeoOut) indMovedAtom(ii), &
                &pCoord0Out(:, indMovedAtom(ii))
          end do
          write (fdUser,*) ""
        end if
        
        !! Write out atomic charges
        if (tMulliken) then
          write (fdUser, "(/,A)") " Net atomic charges (e)"
          write (fdUser, "(1X,A5,1X,A16)")" Atom", " Net charge"
          do ii = 1, nAtom
            write (fdUser, "(1X,I5,1X,F16.8)") ii, &
                &sum(q0(:, ii, 1) - qOutput(:, ii, 1))
          end do
          write(fdUser,*) ""
        end if
        lpSpinPrint: do iSpin = 1, mod(nSpin,3)
          if (nSpin == 2) then
            write(fdUser,*) 'COMPONENT = ',trim(spinName(iSpin))
          else
            write(fdUser,*) 'COMPONENT = ',trim(quaternionName(iSpin))
          end if
          write(fdUser,*) ' '
          write(fdUser,*)'Eigenvalues /H'
          do iEgy=1,nOrb*ceiling(0.5*real(nSpin))
            write(fdUser, formatEigen) (eigen(iEgy,ii,iSpin),ii=1,nKPoint)
          end do
          write(fdUser,*) ""
          write(fdUser,*)'Eigenvalues /eV'
          do iEgy=1,nOrb*ceiling(0.5*real(nSpin))
            write(fdUser, formatEigen) &
                &(Hartree__eV*eigen(iEgy,ii,iSpin),ii=1,nKPoint)
          end do
          write (fdUser,*) ''
          write(fdUser,*)'Fillings'
          do iEgy=1,nOrb*ceiling(0.5*real(nSpin))
            write(fdUser, formatEnergy) (filling(iEgy,ii,iSpin),ii=1,nKPoint)
          end do
          write (fdUser,*) ""
        end do lpSpinPrint
        if (nSpin == 4) then
          if (tMulliken) then
            do jj = 1, 4
              write (fdUser,"(' Nr. of electrons (',A,'):',F16.8)") &
                  & quaternionName(jj),sum(qOutput(:, :,jj))
              write (fdUser,*) ""
              write (fdUser, "(' Atom populations (',A,')')") quaternionName(jj)
              write (fdUser, "(1X,A5,1X,A16)")" Atom", " Population"

              do ii = 1, nAtom
                write (fdUser, "(1X,I5,1X,F16.8)") ii, sum(qOutput(:, ii, jj))
              end do
              write (fdUser,*) ''
              write (fdUser, "(' l-shell populations (',A,')')") quaternionName(jj)
              write (fdUser, "(1X,A5,1X,A3,1X,A3,1X,A16)")" Atom", "Sh.", &
                  &"  l", " Population"
              do ii = 1, nAtom
                iSp1 = specie(ii)
                do iSh1 = 1, orb%nShell(iSp1)
                  write (fdUser, "(1X,I5,1X,I3,1X,I3,1X,F16.8)") ii, iSh1, &
                      &orb%angShell(iSh1,iSp1), &
                      &sum(qOutput(orb%posShell(iSh1,iSp1)&
                      &:orb%posShell(iSh1+1,iSp1)-1, ii, jj))
                end do
              end do
              write (fdUser,*) ''
              write (fdUser, "(' Orbital populations (',A,')')") &
                  & quaternionName(jj)
              write (fdUser, "(1X,A5,1X,A3,1X,A3,1X,A3,1X,A16)") " Atom", &
                  & "Sh.","  l","  m", " Population"
              do ii = 1, nAtom
                iSp1 = specie(ii)
                do iSh1 = 1, orb%nShell(iSp1)
                  ang = orb%angShell(iSh1, iSp1)
                  do kk = 0, 2 * ang
                    write (fdUser, "(' ',I5,1X,I3,1X,I3,1X,I3,1X,F16.8)") &
                        &ii, iSh1, ang, kk - ang, &
                        &qOutput(orb%posShell(iSh1,iSp1)+kk, ii, jj)
                  end do
                end do
              end do
              write (fdUser,*) ''
            end do
          end if
          if (tDFTBU) then
            do jj = 1, 4
              write (fdUser, "(' Block populations (',A,')')") &
                  & quaternionName(jj)
              do ii = 1, nAtom
                iSp1 = specie(ii)
                write(fdUser,*)'Atom',ii
                do kk = 1, orb%nOrbSpecie(iSp1)
                  write(fdUser,"(16F8.4)") &
                      & qBlockOut(1:orb%nOrbSpecie(iSp1),kk,ii,jj)
                end do
                write (fdUser,*) ''
              end do
            end do
          end if
          if (tImHam .and. tMulliken) then
            write (fdUser,*) ''
            write (fdUser,*) ' Electron angular momentum (mu_B/hbar)'
            write (fdUser, "(2X,A5,T10,A3,T14,A1,T20,A1,T35,A9)")"Atom", "Sh.",&
                &"l", "S", "Momentum"
            do ii = 1, nAtom
              iSp1 = specie(ii)
              do iSh1 = 1, orb%nShell(iSp1)
                write (fdUser, "(1X,I5,1X,I3,1X,I3,1X,F14.8,' :',3F14.8)") &
                    & ii, iSh1, orb%angShell(iSh1,iSp1), &
                    & 0.5_dp*sqrt(sum(sum(qOutput(orb%posShell(iSh1,iSp1)&
                    & :orb%posShell(iSh1+1,iSp1)-1, ii, 2:4),dim=1)**2)), &
                    & -gfac*0.25_dp*sum(qOutput(orb%posShell(iSh1,iSp1)&
                    & :orb%posShell(iSh1+1,iSp1)-1, ii, 2:4),dim=1)
              end do
            end do
            write (fdUser,*) ''
            write (fdUser,*) ' Orbital angular momentum (mu_B/hbar)'
            write (fdUser, "(2X,A5,T10,A3,T14,A1,T20,A1,T35,A9)")"Atom", "Sh.",&
                &"l", "L", "Momentum"
            do ii = 1, nAtom
              iSp1 = specie(ii)
              do iSh1 = 1, orb%nShell(iSp1)                
                write (fdUser, "(1X,I5,1X,I3,1X,I3,1X,F14.8,' :',3F14.8)") &
                    & ii, iSh1, orb%angShell(iSh1,iSp1), &
                    & sqrt(sum(orbitalL(1:3,iSh1,ii)**2)),-orbitalL(1:3,iSh1,ii)
              end do
            end do
            write (fdUser,*) ''
            write (fdUser,*) ' Total angular momentum (mu_B/hbar)'
            write (fdUser, "(2X,A5,T10,A3,T14,A1,T20,A1,T35,A9)")"Atom", "Sh.",&
                &"l", "J", "Momentum"
            angularMomentum = 0.0_dp
            do ii = 1, nAtom
              iSp1 = specie(ii)
              do iSh1 = 1, orb%nShell(iSp1)
                write (fdUser, "(1X,I5,1X,I3,1X,I3,1X,F14.8,' :',3F14.8)") &
                    & ii, iSh1, orb%angShell(iSh1,iSp1), sqrt(sum((&
                    & orbitalL(1:3,iSh1,ii) &
                    & +sum(0.5_dp*qOutput(orb%posShell(iSh1,iSp1)&
                    & :orb%posShell(iSh1+1,iSp1)-1, ii, 2:4),dim=1))**2)), &
                    & -orbitalL(1:3,iSh1,ii) &
                    & -gfac*0.25_dp*sum(qOutput(orb%posShell(iSh1,iSp1)&
                    & :orb%posShell(iSh1+1,iSp1)-1, ii, 2:4),dim=1)
                angularMomentum(1:3) = angularMomentum(1:3) &
                    &  -orbitalL(1:3,iSh1,ii) &
                    & -gfac*0.25_dp*sum(qOutput(orb%posShell(iSh1,iSp1)&
                    & :orb%posShell(iSh1+1,iSp1)-1, ii, 2:4),dim=1)
              end do
            end do
            write (fdUser,*) ''
          end if
        else
          if (nSpin == 2) then
            call qm2ud(qOutput)
            if (tDFTBU) then
              call qm2ud(qBlockOut)
            end if
          end if
          lpSpinPrint2: do iSpin = 1, nSpin
            if (tMulliken) then
              write (fdUser,"(' Nr. of electrons (',A,'):',F16.8)") &
                  & trim(spinName(iSpin)),sum(qOutput(:, :,iSpin))
              write (fdUser, "(' Atom populations (',A,')')")&
                  & trim(spinName(iSpin))
              write (fdUser, "(1X,A5,1X,A16)")" Atom", " Population"

              do ii = 1, nAtom
                write (fdUser, "(1X,I5,1X,F16.8)")ii, sum(qOutput(:, ii, iSpin))
              end do
              write (fdUser,*) ''
              write (fdUser, "(' l-shell populations (',A,')')") &
                  & trim(spinName(iSpin))
              write (fdUser, "(1X,A5,1X,A3,1X,A3,1X,A16)")" Atom", "Sh.", &
                  &"  l", " Population"
              do ii = 1, nAtom
                iSp1 = specie(ii)
                do iSh1 = 1, orb%nShell(iSp1)
                  write (fdUser, "(1X,I5,1X,I3,1X,I3,1X,F16.8)") ii, iSh1, &
                      &orb%angShell(iSh1,iSp1), &
                      &sum(qOutput(orb%posShell(iSh1,iSp1)&
                      &:orb%posShell(iSh1+1,iSp1)-1, ii, iSpin))
                end do
              end do
              write (fdUser,*) ''
              write (fdUser, "(' Orbital populations (',A,')')") &
                  & trim(spinName(iSpin))
              write (fdUser, "(1X,A5,1X,A3,1X,A3,1X,A3,1X,A16)")" Atom", "Sh.",&
                  &"  l","  m", " Population"
              do ii = 1, nAtom
                iSp1 = specie(ii)
                do iSh1 = 1, orb%nShell(iSp1)
                  ang = orb%angShell(iSh1, iSp1)
                  do kk = 0, 2 * ang
                    write (fdUser, "(' ',I5,1X,I3,1X,I3,1X,I3,1X,F16.8)") &
                        &ii, iSh1, ang, kk - ang, &
                        &qOutput(orb%posShell(iSh1,iSp1)+kk, ii, iSpin)
                  end do
                end do
              end do
              write (fdUser,*) ''
            end if
            if (tDFTBU) then
              write (fdUser, "(' Block populations (',A,')')") &
                  & trim(spinName(iSpin))
              do ii = 1, nAtom
                iSp1 = specie(ii)
                write(fdUser,*)'Atom',ii
                do kk = 1, orb%nOrbSpecie(iSp1)
                  write(fdUser,"(16F8.4)") &
                      & qBlockOut(1:orb%nOrbSpecie(iSp1),kk,ii,iSpin)
                end do
              end do
              write (fdUser,*) ''
            end if
          end do lpSpinPrint2
          if (nSpin == 2) then
            call ud2qm(qOutput)
            if (tDFTBU) then
              call ud2qm(qBlockOut)
            end if
          end if
        end if

        if (nSpin == 2) then
          call qm2ud(qOutput)
          call qm2ud(qInput)
        end if
        lpSpinPrint3: do iSpin = 1, mod(nSpin,3)
          if (nSpin == 2) then
            write(fdUser,*)'Spin ',trim(spinName(iSpin))
          end if
          write(fdUser,format2U) 'Fermi energy', Ef(iSpin),"H", &
              & Hartree__eV*Ef(iSpin),'eV'
          write(fdUser,format2U) 'Band energy', Eband(iSpin),"H", &
              &Hartree__eV*Eband(iSpin),'eV'
          write(fdUser,format2U)'TS', TS(iSpin),"H", Hartree__eV*TS(iSpin),'eV'
          write(fdUser,format2U) 'Band free energy (E-TS)',&
              &Eband(iSpin)-TS(iSpin),"H",&
              &Hartree__eV*(Eband(iSpin)-TS(iSpin)),'eV'
          write(fdUser,format2U)'Extrapolated E(0K)',E0(iSpin),"H",&
              &Hartree__eV*(E0(iSpin)),'eV'
          if (tMulliken) then
            if (nSpin == 2) then
              write(fdUser, &
                  & "(' Input/Output electrons (',A,'):',F16.8,F16.8)") &
                  &trim(spinName(iSpin)), sum(qInput(:, :, iSpin)), &
                  & sum(qOutput(:, :,iSpin))
            else
              write(fdUser, &
                  & "(' Input/Output electrons (',A,'):',F16.8,F16.8)") &
                  & quaternionName(iSpin), sum(qInput(:, :, iSpin)), &
                  & sum(qOutput(:, :,iSpin))
            end if
          end if
          write (fdUser,*) ''

        end do lpSpinPrint3
        if (nSpin == 2) then
          call ud2qm(qOutput)
          call ud2qm(qInput)
        end if

        write(fdUser,format2U) 'Energy H0', energy%EnonSCC,'H',  &
            & energy%EnonSCC*Hartree__eV,'eV'

        if (tSCC) then
          write (fdUser,format2U) 'Energy SCC', energy%ESCC,'H', &
              & energy%ESCC*Hartree__eV,'eV'
          if (tSpin) then
            write (fdUser,format2U) 'Energy SPIN', energy%Espin,'H', &
                & energy%Espin*Hartree__eV,'eV'
          end if
          if (tDFTBU) then
            write (fdUser,format2U) 'Energy DFTB+U', energy%Edftbu,'H', &
                &energy%Edftbu*Hartree__eV,'eV'
          end if
        end if
        if (tSpinOrbit) then
          write(fdUser,format2U) 'Energy L.S', energy%ELS,'H', &
              & energy%ELS*Hartree__eV,'eV'
        end if
        if (tEfield) then
          write(fdUser,format2U) 'Energy ext. field', energy%Eext,'H', &
              & energy%Eext*Hartree__eV,'eV'
        end if
        
        write(fdUser,format2U) 'Total Electronic energy', energy%Eelec,'H', &
            & energy%Eelec*Hartree__eV,'eV'

        !! Write out repulsive related data

        write(fdUser,format2U) 'Repulsive energy', energy%Erep,'H', &
            & energy%Erep*Hartree__eV,'eV'
        if (tDispersion) then
          write (fdUser,format2U) 'Dispersion energy', energy%eDisp,'H', &
              &energy%eDisp * Hartree__eV,'eV'
        end if
        write(fdUser,format2U) 'Total energy', energy%Etotal,'H', &
            &(energy%Etotal) * Hartree__eV,'eV'
        write(fdUser,format2U) 'Total Mermin free energy', energy%Etotal -&
            & sum(TS(:)),'H', (energy%Etotal - sum(TS(:))) * Hartree__eV,'eV'
        if (tPeriodic.and.pressure/=0.0_dp) then
          write(fdUser,format2U) 'Gibbs free energy', energy%Etotal -&
              & sum(TS(:))+cellVol*pressure,'H', Hartree__eV * &
              & (energy%Etotal - sum(TS(:))+cellVol*pressure),'eV'
        end if
        write (fdUser,*) ''

        if (tAtomicEnergy) then
          write (fdUser, "(' Atom resolved electronic energies ')")
          do ii = 1, nAtom
            write (fdUser, "(1X,I5,F16.8,' H',F16.6,' eV')") ii, &
                & energy%atomElec(ii), Hartree__eV*energy%atomElec(ii)
          end do
          write (fdUser,*) ''
        end if

        if (tAtomicEnergy) then
          write (fdUser, "(' Atom resolved repulsive energies ')")
          do ii = 1, nAtom
            write (fdUser, "(1X,I5,F16.8,' H',F16.6,' eV')") ii, &
                & energy%atomRep(ii), Hartree__eV*energy%atomRep(ii)
          end do
          write (fdUser,*) ''
          write (fdUser, "(' Atom resolved total energies ')")
          do ii = 1, nAtom
            write (fdUser, "(1X,I5,F16.8,' H',F16.6,' eV')") ii, &
                & energy%atomTotal(ii), Hartree__eV*energy%atomTotal(ii)
          end do
          write (fdUser,*) ''
        end if

      end if

      if (tConverged) then
        exit lpSCC
      end if

      iSCCIter = iSCCIter + 1

    end do lpSCC

    if (tEigenvecs) then
      if (nSpin == 4) then
        call warning("Eigenvector writing not supported for non-colinear spin.")
      else
        if (tRealHS) then
          if (tStoreEigvecs) then
            call writeRealEigvecs(fdEigvec, runId, nAtom, nSpin, neighborList, &
                &nNeighbor, iAtomStart, iPair, img2CentCell, orb, specie, &
                &specieName, over, HSqrReal, SSqrReal, storeEigvecsReal)
          else
            call writeRealEigvecs(fdEigvec, runId, nAtom, nSpin, neighborList, &
                &nNeighbor, iAtomStart, iPair, img2CentCell, orb, specie, &
                &specieName, over, HSqrReal, SSqrReal)
          end if          
        else
          if (tStoreEigvecs) then
            call writeCplxEigvecs(fdEigvec, runId, nAtom, nSpin, neighborList, &
                &nNeighbor, cellVec, iCellVec, iAtomStart, iPair, img2CentCell,&
                &orb, specie, specieName, over, kpoint, HSqrCplx, SSqrCplx, &
                &storeEigvecsCplx)
          else
            call writeCplxEigvecs(fdEigvec, runId, nAtom, nSpin, neighborList,&
                &nNeighbor, cellVec, iCellVec, iAtomStart, iPair, img2CentCell,&
                &orb, specie, specieName, over, kpoint, HSqrCplx, SSqrCplx)
          end if
        end if
      end if
    end if

    if (tProjEigenvecs) then
      if (nSpin == 4) then
        call warning("Eigenvector projection not supported for non-colinear &
            &spin.")
      else
        if (tRealHS) then
          if (tStoreEigvecs) then
            call writeProjRealEigvecs(regionLabels, eigen, nSpin, neighborList, &
                & nNeighbor, iAtomStart, iPair, img2CentCell, orb, over, &
                & HSqrReal, SSqrReal, iOrbRegion, storeEigvecsReal)
          else
            call writeProjRealEigvecs(regionLabels, eigen, nSpin, neighborList, &
                & nNeighbor, iAtomStart, iPair, img2CentCell, orb, over, &
                & HSqrReal, SSqrReal, iOrbRegion)
          end if
        else
          if (tStoreEigvecs) then
            call writeProjCplxEigvecs(regionLabels, eigen, nSpin, neighborList, &
                & nNeighbor, cellVec, iCellVec, iAtomStart, iPair, &
                & img2CentCell, orb, over, kpoint, kWeight, HSqrCplx, &
                & SSqrCplx, iOrbRegion, storeEigvecsCplx)
          else
            call writeProjCplxEigvecs(regionLabels, eigen, nSpin, neighborList, &
                & nNeighbor, cellVec, iCellVec, iAtomStart, iPair, &
                & img2CentCell, orb, over, kpoint, kWeight, HSqrCplx, &
                & SSqrCplx, iOrbRegion)
          end if
        end if
      end if
    end if
    
    if (tGeoOpt) then
      write (lcTmp, "(A,A)") trim(geoOutFile), ".xyz"
      if (.not. tAppendGeo) then
        call clearFile(trim(lcTmp))
      end if
      if (tLatOpt) then
        write (tmpStr, "('** Geometry step: ',I0,', Lattice step: ',I0)") &
          & iGeoStep,iLatGeoStep
      else
        write(tmpStr,"('Geometry Step: ',i0)")iGeoStep
      end if

      if (tMulliken) then
        if (nSpin == 4) then
          ALLOCATE_(tmpMatrix,(3,nAtom))
          do jj = 1, nAtom
            do ii = 1, 3
              tmpMatrix(ii,jj) = sum(qOutput(:,jj,ii+1))
            end do
          end do
          ! convert by the inverse of the scaling used in writeXYZFormat :
          tmpMatrix(:,:) = tmpMatrix(:,:) * au__fs / (1000_dp * Bohr__AA)
          call writeXYZFormat(trim(lcTmp), pCoord0Out, specie0, specieName, &
              &charges=sum(qOutput(:,:,1),dim=1), velocities = tmpMatrix, &
              & comment=trim(tmpStr))
          DEALLOCATE_(tmpMatrix)
        else
          call writeXYZFormat(trim(lcTmp), pCoord0Out, specie0, specieName, &
              &comment=trim(tmpStr))
        end if
      else
        call writeXYZFormat(trim(lcTmp), pCoord0Out, specie0, specieName, &
            &comment=trim(tmpStr))
      end if
    end if

    write (*,*)
    write (*, format1U) "Total Energy", energy%Etotal,"H"
    write (*, format1U) "Total Mermin free energy", &
        &energy%Etotal - sum(TS),"H"

    if (tDipole) then
      dipoleMoment(:) = 0.0_dp
      do iAtom = 1, nAtom
        dipoleMoment(:) = dipoleMoment(:) &
            & + sum(q0(:, iAtom, 1) - qOutput(:, iAtom, 1)) * coord(:,iAtom)
      end do
#if DEBUG >= 1
      ! extra test for the potential in the code, does the dipole from
      ! charge positions match the derivative of energy wrt an external E field?
      ALLOCATE_(hprime,(size(h0),1))
      ALLOCATE_(dipoleTmp,(size(qOutput,dim=1),nAtom))
      ALLOCATE_(potentialDerivative,(nAtom,1))
      write(*,"(A)",advance='no')'Hellmann Feynman dipole:'
      do ii = 1, 3 ! loop over directions
        potentialDerivative = 0.0_dp
        potentialDerivative(:,1) = -coord(ii,:) ! Potential from dH/dE
        hprime = 0.0_dp
        dipoleTmp = 0.0_dp
        call add_shift(hprime,over,nNeighbor, neighborList%iNeighbor, &
            & specie,orb,iPair,nAtom,img2CentCell,potentialDerivative)
        ! evaluate <psi| dH/dE | psi>
        call mulliken(dipoleTmp, hprime(:,1), rhoPrim(:,1), & 
            &orb, neighborList%iNeighbor, nNeighbor, img2CentCell, iPair)
        ! add nuclei term for derivative wrt E
        do iAtom = 1, nAtom
          dipoleTmp(1,iAtom) = dipoleTmp(1,iAtom) &
              & + sum(q0(:,iAtom,1))*coord0(ii,iAtom)
        end do
        write(*,"(f12.8)",advance='no')sum(dipoleTmp)
      end do
      write(*,*)" au"
      DEALLOCATE_(potentialDerivative)
      DEALLOCATE_(hprime)
      DEALLOCATE_(dipoleTmp)
#endif
    else
      dipoleMoment(:) = 0.0_dp
    end if

    !! Calculate energy weighted density matrix
    if (tForces) then
      
      !! Calculate the identity part of the energy weighted density matrix
      ERhoPrim(:) = 0.0_dp

      if (nSpin == 4) then
        do nK = 1, nKPoint
          !! Calculate eigenvectors, if necessary. Eigenvectors for the last
          !! last k-points are still there, so use those directly
          if (tStoreEigvecs) then
            iK2 = 1
            call get(storeEigvecsCplx(1), HSqrCplx(:,:,iK2, 1))
          else
            iK2 = nK
          end if

          if (tDensON) then
            call error("Currently missing.")
          else
            call makeDensityMatrix(SSqrCplx, HSqrCplx(:,:,iK2,1), &
                &filling(:,nk,1), eigen(:,nk,1))
          end if
          if (tRealHS) then
            call packERho(ERhoPrim(:), SSqrCplx, neighborList%iNeighbor, &
                &nNeighbor, orb%mOrb, iAtomStart, iPair, img2CentCell)
          else
            call packERho(ERhoPrim(:), SSqrCplx, kPoint(:,nk), &
                &kWeight(nk), neighborList%iNeighbor, nNeighbor, orb%mOrb, &
                &iCellVec, cellVec, iAtomStart, iPair, img2CentCell)
          end if
        end do
      else
        if (tRealHS) then
          do iSpin = 1, nSpin

            if (tStoreEigvecs) then
              iSpin2 = 1
              call get(storeEigvecsReal(mod(iSpin,3)), &
                  & HSqrReal(:,:,mod(iSpin2,3)))
            else
              iSpin2 = iSpin
            end if

            ! Build energy weighted density matrix
            if (tDensON) then
              call makeDensityMatrix(SSqrReal, HSqrReal(:,:,iSpin2), &
                  &filling(:,1,iSpin), eigen(:,1,iSpin), &
                  & neighborlist%iNeighbor, nNeighbor, orb, iAtomStart, &
                  & img2CentCell)
            else
              call makeDensityMatrix(SSqrReal, HSqrReal(:,:,iSpin2), &
                  &filling(:,1,iSpin), eigen(:,1,iSpin))
            end if
            call packHS(ERhoPrim(:), SSqrReal, neighborList%iNeighbor, &
                &nNeighbor, orb%mOrb, iAtomStart, iPair, img2CentCell)
          end do
        else
          do iSpin = 1, nSpin
            if (tStoreEigvecs) then
              iSpin2 = 1
            else
              iSpin2 = iSpin
            end if

            do nK = 1, nKPoint
              !! Calculate eigenvectors, if necessary. Eigenvectors for the last
              !! spin in the last k-points are still there, so use those
              !! directly
              if (tStoreEigvecs) then
                iK2 = 1
                call get(storeEigvecsCplx(iSpin), HSqrCplx(:,:,iK2, iSpin2))
              else
                iK2 = nK
              end if

              if (tDensON) then
                call makeDensityMatrix(SSqrCplx, HSqrCplx(:,:,iK2,iSpin2), &
                    &filling(:,nK,iSpin), eigen(:,nK, iSpin), &
                    &neighborlist%iNeighbor, nNeighbor, orb, iAtomStart,&
                    &img2CentCell)
              else
                call makeDensityMatrix(SSqrCplx, HSqrCplx(:,:,iK2,iSpin2), &
                    &filling(:,nK,iSpin), eigen(:,nK, iSpin))
              end if
              call packHS(ERhoPrim(:), SSqrCplx, kPoint(:,nk), &
                  &kWeight(nk), neighborList%iNeighbor, nNeighbor, orb%mOrb, &
                  &iCellVec, cellVec, iAtomStart, iPair, img2CentCell)
            end do
          end do
        end if
      end if

      derivs(:,:) = 0.0_dp
      if (tExtChrg) then
        chrgForces(:,:) = 0.0_dp
      end if

      if (.not. tSCC) then
        if (tImHam) then
          call derivative_shift(derivs,rhoPrim,iRhoPrim,erhoPrim(:),skHamCont, &
              & skOverCont, coord, specie, neighborList%iNeighbor, nNeighbor, &
              & img2CentCell, iPair, orb, potential%intBlock, &
              & potential%iorbitalBlock)
        else
          call derivative_nonscc(derivs,rhoPrim(:,1), ERhoPrim(:), skHamCont, &
              &skOverCont, coord, specie, neighborList%iNeighbor, nNeighbor, &
              &img2CentCell, iPair, orb)
        end if        
      else
        potential%intBlock = potential%intBlock + potential%extBlock
        if (tDFTBU) then
          potential%intBlock = potential%orbitalBlock + potential%intBlock
        end if

        if (tImHam) then
          call derivative_shift(derivs,rhoPrim,iRhoPrim,erhoPrim(:),skHamCont, &
              & skOverCont, coord, specie, neighborList%iNeighbor, nNeighbor, &
              & img2CentCell, iPair, orb, potential%intBlock, &
              & potential%iorbitalBlock)
        else
          call derivative_shift(derivs,rhoPrim,erhoPrim(:),skHamCont, &
              &skOverCont, coord, specie, neighborList%iNeighbor, nNeighbor, &
              &img2CentCell, iPair, orb, potential%intBlock)
        end if

        ! add double counting terms in force :
        if (tExtChrg) then
          call addForceDCSCC(derivs, specie, neighborList%iNeighbor, &
              &img2CentCell, coord, chrgForces)
        else
          call addForceDCSCC(derivs, specie, neighborList%iNeighbor, &
              &img2CentCell, coord)
        end if
        if (tEField) then
          do ii = 1, 3
            derivs(ii,:) = derivs(ii,:) - &
                & sum(q0(:, :, 1)-qOutput(:,:,1),dim=1)*EField(ii)
          end do
        end if
        if (t3rdFull) then
          call addGradientDC(thirdOrd, neighborList, specie, coord, &
              &img2CentCell, derivs)
        end if
      end if

      call getERepDeriv(repulsiveDerivs, coord, nNeighbor, &
          &neighborList%iNeighbor, specie,pRepCont, img2CentCell)
      totalDeriv(:,:) = repulsiveDerivs(:,:) + derivs(:,:)

      if (tDispersion) then
        call addGradients(myDispersion, totalDeriv)
      end if
      
      if (tStress) then       
        call getRepulsiveStress(repulsiveStress, coord, nNeighbor, &
            & neighborList%iNeighbor, specie, img2CentCell, pRepCont, cellVol)
        if (tSCC) then
          if (tImHam) then
            call getBlockiStress(elecStress,rhoPrim,iRhoPrim,ERhoPrim, &
                & skHamCont, skOverCont, coord, specie, &
                & neighborList%iNeighbor, nNeighbor,&
                & img2CentCell, iPair, orb, potential%intBlock, &
                & potential%iorbitalBlock, cellVol)
          else
            call getBlockStress(elecStress, rhoPrim, ERhoPrim, skHamCont, &
                &skOverCont, coord, specie, neighborList%iNeighbor,nNeighbor,&
                &img2CentCell, iPair, orb, potential%intBlock, cellVol)
          end if
          
          call addStressDCSCC(elecStress,specie,neighborList%iNeighbor, &
              & img2CentCell,coord)
          
          if (tEField) then
            elecLatDeriv = 0.0_dp
            call cart2frac(coord0,latVec)
            do iAtom = 1, nAtom
              do ii = 1, 3
                do jj = 1, 3
                  elecLatDeriv(jj,ii) =  elecLatDeriv(jj,ii) - &
                      & sum(q0(:, iAtom, 1)-qOutput(:,iAtom,1),dim=1) &
                      & * EField(ii) * coord0(jj,iAtom)
                end do
              end do
            end do
            call frac2cart(coord0,latVec)
            elecStress = elecStress &
                & -matmul(elecLatDeriv,transpose(latVec))/cellVol
          end if
          
        else
          if (tImHam) then
            call getBlockiStress(elecStress,rhoPrim,iRhoPrim,ERhoPrim, &
                & skHamCont, skOverCont, coord, specie, &
                & neighborList%iNeighbor, nNeighbor,&
                & img2CentCell, iPair, orb, potential%intBlock, &
                & potential%iorbitalBlock, cellVol)
          else
            call getNonSCCStress(elecStress,rhoPrim(:,1), ERhoPrim, skHamCont, &
                &skOverCont, coord, specie, neighborList%iNeighbor, nNeighbor,&
                &img2CentCell, iPair, orb, cellVol)
          end if
        end if
        
        if (tDispersion) then
          call getStress(myDispersion,dispStress)
          dispLatDeriv = -cellVol * matmul(dispStress,invLatVec)
        end if
        repulsiveLatDeriv = -cellVol * matmul(repulsiveStress,invLatVec)
        elecLatDeriv = -cellVol * matmul(elecStress,invLatVec)
        
        totalStress = repulsiveStress + elecStress + dispStress 
        cellPressure = ( totalStress(1,1) + totalStress(2,2) &
            & + totalStress(3,3) )/3.0_dp
        totalLatDeriv = repulsiveLatDeriv + elecLatDeriv + dispLatDeriv
        
        write(*,format1Ue)'Volume',cellVol,'au^3'
        write(*,format2Ue)'Pressure',cellPressure,'au',&
            & cellPressure* au__pascal, 'Pa'
      end if
    end if

    if (tPeriodic.and.pressure/=0.0_dp) then
      write (*, format1U) "Gibbs free energy", &
          & energy%Etotal - sum(TS(:)) + cellVol*pressure,'au'
    end if

    !! Write out information at the end of the SCC loop

    if (tWriteDetailedOut) then
      if (tSCC) then
        if (tConverged) then
          write (fdUser,*) "SCC converged"
          write (fdUser,*) ""
        else
          write (fdUser,*)"SCC is NOT converged, maximal SCC&
              & iterations exceeded"
          write (fdUser,*) ""
          if (tConvrgForces) then
            call error("SCC is NOT converged, maximal SCC iterations&
                & exceeded")
          else
            call warning("SCC is NOT converged, maximal SCC iterations&
                & exceeded")
          end if
        end if
      else
        write (fdUser,*) "Non-SCC calculation"
        write (fdUser,*) ""
      end if

      if (tGeoOpt .or. tMD) then
        write(fdUser,*) "Full geometry written in ",trim(geoOutFile), &
            & ".{xyz|gen}"
        write(fdUser,*)''
      end if

      !! Write out forces
      if (tForces) then

        write(fdUser,*)'Total Forces'
        do ii = 1, nAtom
          write(fdUser,*) -totalDeriv(:,ii)
        end do
        write(fdUser,*)''

        if (tStress.and. .not.tMD) then
          write(fdUser,*)'Total stress tensor'
          do ii = 1, 3
            write(fdUser,"(3F20.12)")totalStress(:,ii)
          end do
          write(fdUser,*)
          write(fdUser,*)'Total lattice derivs'
          do ii = 1, 3
            write(fdUser,"(3F20.12)")totalLatDeriv(:,ii)
          end do
          write(fdUser,*)''
        end if

        write(fdUser,format1Ue) "Maximal derivative component", &
            & maxval(abs(totalderiv)),'au'
        if (nMovedAtom > 0) then
          write (fdUser,format1Ue) "Max force for moved atoms:", &
              &maxval(abs(totalDeriv(:,indMovedAtom))),'au'
        end if
        write(fdUser,*)''

        if (tExtChrg) then
          write (fdUser,*) "Forces on external charges"
          do ii = 1, nExtChrg
            write (fdUser, *) -chrgForces(:,ii)
          end do
          write(fdUser,*)''
        end if

        if (tPeriodic .and. .not. tMD) then
          write(fdUser,format1Ue)'Volume',cellVol,'au^3'
          if (tStress) then
            write(fdUser,format2Ue)'Pressure',cellPressure, 'au', &
                & cellPressure * au__pascal, 'Pa'
          end if
          write(fdUser,*)''
        end if

      end if
    end if

    if (tForces) then
      !! Set force components along constraint vectors zero
      do ii = 1, nGeoConstr
        iAtom = conAtom(ii)
        totalDeriv(:,iAtom) = totalDeriv(:,iAtom) &
            &- conVec(:,ii) * dot_product(conVec(:,ii), totalDeriv(:,iAtom))
      end do

      if (tCoordOpt) then
        tmpDerivs(1:nMovedCoord) = &
            & reshape(totalDeriv(:,indMovedAtom),(/ nMovedCoord /))
        write (*,"(' ',A,':',T30,E20.6)") "Maximal force component", &
            &maxval(abs(tmpDerivs))
      end if

      if (tLatOpt) then
        tmpLat3Vecs = totalLatDeriv + derivCellVol
        tmpLatVecs(1:9) = reshape(tmpLat3Vecs,(/ 9 /))

        if (tLatOptFixAng) then ! project forces to be along original lattice
          tmpLat3Vecs = tmpLat3Vecs * origLatVec
          tmpLatVecs(:) = 0.0_dp
          if (any(tLatOptFixLen)) then
            tmpLatVecs(:) = 0.0_dp
            do ii = 1, 3
              if (.not.tLatOptFixLen(ii)) then
                tmpLatVecs(ii) = sum(tmpLat3Vecs(:,ii))
              end if
            end do
          else
            tmpLatVecs(1:3) = sum(tmpLat3Vecs,dim=1)
          end if
        elseif (tLatOptIsotropic) then
          tmpLat3Vecs = tmpLat3Vecs * origLatVec
          tmpLatVecs(:) = 0.0_dp
          tmpLatVecs(1) = sum(tmpLat3Vecs)
        end if
        write (*,format1Ue) "Maximal Lattice force component", &
            &maxval(abs(tmpLatVecs)),'au'
      end if

      !! If geometry minimizer finished and the last calculated geometry is the
      !! minimal one (not necessary the case, depends on the optimizer!)
      !! -> we are finished.
      !! Otherwise we have to recalc everything in the converged geometry.

      if (tGeomEnd) then
        exit lpGeomOpt
      else
        if (tWriteRestart .and. tMulliken .and. tSCC) then
          if (tDFTBU) then
            if (tSpinOrbit) then
              call writeQToFile(qInput, fChargeIn, orb, qBlockIn, qiBlockIn)
            else
              call writeQToFile(qInput, fChargeIn, orb, qBlockIn)
            end if
          else
            call writeQToFile(qInput, fChargeIn, orb)
          end if
          print "('>> Charges saved for restart in ',A)", fChargeIn
        end if
        if (tDerivs) then
          call next(pDerivDriver,new3Coord,totalDeriv(:,indMovedAtom),tGeomEnd)
          coord0(:,indMovedAtom) = new3Coord(:,:)
          if (tGeomEnd) exit lpGeomOpt
        elseif (tGeoOpt) then
          if (tCoordStep) then
            call next(pGeoCoordOpt, energy%Etotal - sum(TS), tmpDerivs, &
                & tmpCoords,tCoordEnd)
            if (.not.tLatOpt) tGeomEnd = tCoordEnd
          else
            call next(pGeoLatOpt, energy%Etotal - sum(TS) + cellVol*pressure, &
                & tmpLatVecs, newLatVecs,tGeomEnd)
            if (tLatOptFixAng) then ! optimization uses scaling factor of
              !  lattice vectors
              if (any(tLatOptFixLen)) then
                do ii = 3, 1, -1
                  if (.not.tLatOptFixLen(ii)) then
                    newLatVecs(3*ii-2:3*ii) =  newLatVecs(ii)*origLatVec(:,ii)
                  else
                    newLatVecs(3*ii-2:3*ii) =  origLatVec(:,ii)
                  end if
                end do
              else
                newLatVecs(7:9) =  newLatVecs(3)*origLatVec(:,3)
                newLatVecs(4:6) =  newLatVecs(2)*origLatVec(:,2)
                newLatVecs(1:3) =  newLatVecs(1)*origLatVec(:,1)
              end if
            else if (tLatOptIsotropic) then ! optimization uses scaling factor
              !  unit cell
              do ii = 3,1, -1 ! loop downwards as reusing newLatVecs
                newLatVecs(3*ii-2:3*ii) =  newLatVecs(1)*origLatVec(:,ii)
              end do
            end if
            iLatGeoStep = iLatGeoStep + 1
          end if
        elseif(tMD) then
!         BXD here
          if (COMbxd .eqv. .true.) then
             COMreverse = .false.
             if (RealCOM .eqv. .true.) then
                call getCOM(pConnect, coord0,  movedMass, COM1, COM2, distance, constraint )
             else
                call getDistance(pConnect, coord0,  movedMass, COM1, COM2, distance, constraint)
                write(fdBXD,bxdFormat)"COMdist", distance, iGeoStep
                write(fdBXD,bxdFormat)"Energy", (energy%Etotal - sum(TS)),iGeoStep
             end if
             if ((distance .lt. minCOM) .and. (COMbxd .eqv. .true.)) then
                bxdWait = .false.
                if (restrain .eqv. .false.) then
                   COMbxd = .false.
                   call reinitialise(pConnect,coord0)
                   call writeXYZFormat('newMin.xyz', pCoord0Out, specie0, specieName)
                end if
             elseif (( distance .gt. oldDistance) .and. (distance .gt. minCOM) ) then
                COMreverse = .true.
                write(fdBXD,bxdFormat)"COMreversed", distance, iGeoStep
             end if
          end if 
          if (BXD .eqv. .true.) then
             if ((iGeoStep .ge. bxdEquil) .and. (bxdWait .eqv. .false.) .and. (bxdStarted .eqv. .false. )) then
                bxdStarted = .true.
                bxdBounded = .true.
                bxdLow = (energy%Etotal - sum(TS)) - bxdLower
                maxE = bxdLow + maxE
                if (adaptive .eqv. .true.) then
                    bxdUp = bxdLow + 11.0_dp
                else
                    bxdUp = bxdLow + bxdUpper
                end if
                write(fdBXD,bxdFormat)"bxdOn", bxdLow, iGeoStep
             end if
          end if 
             if (bxdBounded .eqv. .true.) then
                if ((energy%Etotal - sum(TS)) < bxdLow) then
                     DOSstuck = 0
                     write(fdBXD,bxdFormat4)"Lower", bxdLow, bxdUp, iGeoStep
                     write(fdBXD,bxdFormat2)"BoxNumb", currentBox
                     if (currentBox == DOSbox - 1) then
                        DOSreverse =  .false.
                     end if
                     if ((currentBox == DOSbox) .and. (DOSreverse .eqv. .false.)) then
                        upCount = 1
                        print *, DOSbox
                        hitLimit = .false.
                        rateDown = iGeoStep - ratePrev
                        ratePrev = iGeoStep
                        if (downCount .ne. 0) then
                            write (fdDOS,'(I10,I10,I10,I10)'), iGeoStep, rateDown, currentBox, DOSiterations 
                        end if
                        downCount = downCount + 1
                        if (downCount .ge. maxHits) then
                            !bxdUp = bxdUp - 11_dp
                            upCount = maxHits
                            DOSiterations = DOSiterations + 1
                            if ((DOSiterations == DOSiter) .and. (DOS .eqv. .true.)) then
                                tGeomEnd = .true.
                                close(fdDOS)
                                DEALLOCATE_(DOScounts)
                            end if
                        end if
                     end if
                     if (DOSreverse .eqv. .false.) then
                        iGeoStep = iGeoStep - 1
                        stuckCount = stuckCount + 1
                        if (COMreverse .eqv. .true.) then
                            call bxdInvert2(pMDIntegrator,movedMass,constraint,totalOldDeriv)
                        else
                            call bxdInvert(pMDIntegrator,movedMass,totalOldDeriv)
                        end if
                        movedAccel(:,:) = -totalOldDeriv(:,indMovedAtom) / movedMass
                     else
                        bxdUp = bxdLow
                        bxdLow = bxdLow - bxdUpper
                        currentBox = currentBox - 1
                        if (currentBox == DOSbox - 2) then
                            DOSreverse = .false.
                        end if
                        if (COMreverse .eqv. .true.) then
                            call bxdInvert(pMDIntegrator,movedMass,constraint)
                            iGeoStep = iGeoStep - 1
                            stuckCount = stuckCount + 1
                            movedAccel(:,:) = -totalOldDeriv(:,indMovedAtom) / movedMass
                        else
                            movedAccel(:,:) = -totalDeriv(:,indMovedAtom) / movedMass
                        end if
                     end if
                else if (((energy%Etotal - sum(TS)) > bxdUp) .and. (bxdCount .lt. bxdEvents)) then
                     bxdCount = bxdCount + 1
                     if (DOSreverse .eqv. .false.) then
                        if (COMreverse .eqv. .true.) then
                            call bxdInvert2(pMDIntegrator,movedMass,constraint,totalOldDeriv)
                        else
                            call bxdInvert(pMDIntegrator,movedMass,totalOldDeriv)
                        end if
                        iGeoStep = iGeoStep - 1
                        stuckCount = 0
                     end if
                     movedAccel(:,:) = -totalOldDeriv(:,indMovedAtom) / movedMass
                     write(fdBXD,bxdFormat)"Upper", bxdUp, iGeoStep
                     write(fdBXD,bxdFormat2)"BoxNumb", currentBox
                else if (((energy%Etotal - sum(TS)) > bxdUp) .and. (bxdCount >= bxdEvents)) then
                     DOSstuck = DOSstuck + 1
                     if (DOSstuck .ge. 10 .and. DOSreverse .eqv. .true.) then
                         write(fdBXD,bxdFormat2) "Rev", DOSiterations
                         DOSreverse = .false.
                         DOSstuck = 0
                     end if
                     bxdCount = 0
                     if ((currentBox == (bxdBoxes)) .and. (DOS .eqv. .true.)) then
                        DOSreverse = .true.
                        boxCount = 0
                     end if
                     if (currentBox == (DOSbox - 1) .and. (DOSreverse .eqv. .false.) ) then
                        rateUp = iGeoStep - rateUpPrev
                        rateUpPrev = iGeoStep
                        if (upCount .ne. 0) then
                            write (fdDOS2,'(I10,I10,I10,I10)'), iGeoStep, rateUp, currentBox, DOSiterations
                        end if
                        upCount = upCount + 1
                        write(fdBXD,bxdFormat2)"UpCount", upCount
                     end if
                     if ((DOSreverse .eqv. .true.) .or. ((upCount .lt. maxHits) .and. (upCount .gt. 0))) then
                     write(fdBXD,bxdFormat2)"UpCount", upCount
			movedAccel(:,:) = -totalOldDeriv(:,indMovedAtom) / movedMass
                        if (COMreverse .eqv. .true.) then
                            call bxdInvert2(pMDIntegrator,movedMass,constraint, totalOldDeriv)
                        else
                            call bxdInvert(pMDIntegrator,movedMass,totalOldDeriv)
                        end if
                     else
                         movedAccel(:,:) = -totalDeriv(:,indMovedAtom) / movedMass
                         upCount = 0
                         stuckCount = 0
                         bxdLow = bxdUp
                         boxCount = boxCount + 1
                         currentBox = currentBox + 1
                         adaptiveCount = 0
                         highE = -500000.0_dp
                         if (boxCount .ge. bxdBoxes) then
                            adaptive = .false.
                         end if
                         if (((boxCount .ge. bxdBoxes) .and. (DOS .eqv. .false.)) .or. (adaptive .eqv. .true.)) then
                            bxdUp = maxE
                         else if (boxCount .lt. bxdBoxes) then
                            bxdUp = bxdLow + bxdUpper
                         end if
                         if (currentBox == DOSbox) then
                            downCount = 0
                            ratePrev = iGeoStep
                            ! bxdUp = bxdUp + 11_dp
                         end if
                         if (currentBox == (DOSbox - 1)) then
                            rateUpPrev = iGeoStep
                            ! bxdLow = bxdLow - 11_dp
                        end if
                        write(fdBXD,bxdFormat2)"NewBox", iGeoStep
                        write(fdBXD,bxdFormat2)"BoxNumb", currentBox
                        write(fdBXD,bxdFormat5)"NewUpper", bxdUp
                        write(fdBXD, bxdFormat5)"NewLower", bxdLow
                      end if
                else
                     stuckCount = 0
                     if(adaptive .eqv. .true.) then
                        adaptiveCount = adaptiveCount + 1
                        highE = Max(highE,(energy%Etotal - sum(TS)))
                        if (adaptiveCount .gt. adaptiveLimit) then
                           bxdUp = highE
                        end if
                     end if
                     if (COMreverse .eqv. .true.) then
                        movedAccel(:,:) = -totalOldDeriv(:,indMovedAtom) / movedMass
                        call bxdInvert(pMDIntegrator,movedMass,constraint)
                        write(fdBXD,bxdFormat2)"REVERSED", iGeoStep
			write(*,*) 'constraint = ', constraint
                     else
                        movedAccel(:,:) = -totalDeriv(:,indMovedAtom) / movedMass
                     end if 
                end if
        else
            movedAccel(:,:) = -totalDeriv(:,indMovedAtom) / movedMass
        end if


!         Connectivity and reaction checks
!         If we are stuck at a boundary
             if (stuckCount > 100) then
                stuckCount = 0
                bxdLow = bxdLow - 0.001_dp
             end if

          if ( rxnCountDown == 1) then
                call writeXYZFormat('newMin.xyz', pCoord0Out, specie0, specieName)
                iGeoStep = 0
                bxdStarted = .false.
                bxdBounded = .false.
               if (endOnReact .eqv. .true.) then
                  tGeomEnd = .true.
                  close(fdNewMin)
                  close(fdTrans)
                  close(fdBXD)
                  close(fdConform)
               end if
          end if 

          print *, adaptive
          print *, bxdUp                
          call next(pMDIntegrator, movedAccel ,new3Coord, movedVelo, revIdx, rxnCountDown, reversal, movedMass)
          
         totalOldDeriv(:,:) = totalDeriv(:,:)
         oldDistance = distance
         reversal = .false.
         


!        Update connectivities to determine if new reaction has occured
         if (rxnCountDown == 0) then
            call update(pConnect,new3Coord,rxnCountDown,transition, revIdx)
         end if
         if (transition .eqv. .true.) then
              write(fdBXD,bxdFormat)"transition", bxdLow, iGeoStep
              write(fdBXD,bxdFormat)"Window", bxdLow, consitantWindow
              if (consitantWindow == 0) then
                consitantWindow = 100
              else if (consitantWindow == 1) then
                consitantWindow = consitantWindow
                transition = .false.
                rxnCountDown = rxnWindow
                write(fdBXD,bxdFormat)"transition", bxdLow, iGeoStep
                write(fdBXD,bxdFormat)"Window", bxdLow, consitantWindow
                write (fdTrans,"(A15)") "**************"
                write (fdTrans,"(A15,I10)") "Step:", iGeoStep
                write (fdTrans,"(A15)") "ReaxCoord"
                do ii = 1, nAtom
                    if (revIdx(ii) == 1) then
                        write (fdTrans,"(I10)") ii
                    end if
                end do
                call writeXYZFormat('transition.xyz', pCoord0Out, specie0, specieName)
               else
                    consitantWindow = consitantWindow - 1
              end if
          end if

          if (rxnCountDown .GT. 1) then
!          more printing options here?
           rxnCountDown = rxnCountDown - 1 
          end if

          if (associated(pTempProfile)) then
            call next(pTempProfile)
          end if
          call evalKE(KE, movedVelo, movedMass(1,:))
          call evalkT(pMDFrame, kT, movedVelo, movedMass(1,:))

          velocities(:, indMovedAtom) = movedVelo(:,:)
          if (tWriteRestart) then
            write(tmpStr,"('MD iter: ',i0)")iGeoStep
            print *, 'afterFirstWrite'
            print *, tMulliken
            if (tMulliken) then
               call writeXYZFormat(trim(lcTmp), pCoord0Out, specie0, &
                  &specieName, charges=sum(qOutput(:,:,1), dim=1),&
                  &velocities=velocities, comment=trim(tmpStr))
            else
              call writeXYZFormat(trim(lcTmp), pCoord0Out, specie0, &
                 &specieName, velocities=velocities, &
                  &comment=trim(tmpStr))
            end if
          end if
          if (tStress) then
            call getVirialStress(virialStress, mass,specie0,velocities, cellVol)
            virialLatDeriv = -cellVol * matmul(virialStress,invLatVec)

            totalStress(:,:) = repulsiveStress(:,:) + elecStress(:,:) &
                & + virialStress(:,:) + dispStress
            totalLatDeriv(:,:) = repulsiveLatDeriv(:,:) + elecLatDeriv(:,:) &
                & + virialLatDeriv(:,:) + dispLatDeriv
            cellPressure = ( totalStress(1,1) + totalStress(2,2) &
                & + totalStress(3,3) )/3.0_dp
          end if
          if (tMD.and.tStress .and. tWriteDetailedOut) then
            write(fdUser,*)'Total stress tensor'
            do ii = 1, 3
              write(fdUser,"(3F20.12)")totalStress(:,ii)
            end do
            write(fdUser,*)''
            write(fdUser,*)'Total lattice derivs'
            do ii = 1, 3
              write(fdUser,"(3F20.12)")totalLatDeriv(:,ii)
            end do
            write(fdUser,*)''
          end if
          if (tSetFillingTemp) then
            write(*,format2U)'Electronic Temperature:',tempElec,'H',&
                &tempElec/Boltzmann,'K'
          end if
          if (tEfield) then
            write(*,format1U1e)'External E field', absEField, 'au', &
                & absEField * au__V_m, 'V/m'
          end if
          write(*,format2U)"MD Temperature:",kT,"H",kT/Boltzmann,"K"          
          write(*,format1U)"MD Kinetic Energy", KE,"H"
          write(*,format1U)"Total MD Energy",KE + energy%Etotal - sum(TS),"H"

          if (tWriteDetailedOut) then
            if (tSetFillingTemp) then
              write (fdUser, format2U)"Electronic Temperature", tempElec,'au',&
                  &tempElec * Hartree__eV,'eV'
            end if
            write(fdUser,format1U)"MD Kinetic Energy", KE,"H"
            write(fdUser,format1U)"Total MD Energy", &
                & KE + energy%Etotal - sum(TS),"H"
            write(fdUser,format2U)"MD Temperature",kT,"H",kT/Boltzmann,"K"
          end if
        end if

        if (tGeomEnd.and.tGeoOpt) then
          diffGeo = 0.0_dp
          if (tLatOpt) then
            diffGeo = max( maxval(abs(reshape(latVec,(/9/))- newLatVecs)), &
                & diffGeo)
          end if

          if (tCoordOpt) then
            diffGeo = max(maxval(abs(reshape(coord0(:,indMovedAtom), &
                & (/nMovedCoord/))- tmpCoords)),diffGeo)
          end if
          if (diffGeo < tolSameDist) then
            tGeomEnd = .true.
            exit lpGeomOpt
          end if
        end if
        if (.not. tGeomEnd) then
          if (tGeoOpt) then
            if (tCoordStep) then
              if (tCoordEnd) then
                diffGeo = maxval(abs(reshape(coord0(:,indMovedAtom), &
                    & (/nMovedCoord/))- tmpCoords))
                if (diffGeo < tolSameDist) then
                  tCoordStep = .false.
                  if (tLatOpt) then
                    tCoordEnd = .false.
                  end if
                end if
              end if
              if (nMovedCoord > 0) then
                coord0(:,indMovedAtom) = reshape(tmpCoords(:), (/3, nMovedAtom/))
              end if
            else
              call cart2frac(coord0,latVec)
              latVec = reshape(newLatVecs, (/3,3/))
              call frac2cart(coord0,latVec)
              recVec2p = latVec(:,:)
              call matinv(recVec2p)
              recVec2p = reshape(recVec2p, (/3, 3/), order=(/2, 1/))
              recVec = 2.0_dp * pi * recVec2p
              volume = determinant33(latVec)
              recVolume = determinant33(recVec)
              call getCellTranslations(cellVec, rCellVec, latVec, recVec2p, &
                  & mCutoff)
              if (tSCC) then
                call updateLatVecs_SCC(latVec, recVec, volume)
                mCutoff = max(mCutoff, getSCCCutoff())
              end if
              if (tDispersion) then
                call updateLatVecs(myDispersion, latVec, recVec, volume, &
                    &specie0)
                mCutoff = max(mCutoff, getRCutoff(myDispersion))
              end if
              if (tCoordOpt) then
                tCoordStep = .true.
                tCoordEnd = .false.
                tmpCoords(1:nMovedCoord) = reshape(coord0(:, indMovedAtom), &
                    & (/ nMovedCoord /))
                call reset(pGeoCoordOpt, tmpCoords)
              end if
            end if
          elseif (tMD) then
            coord0(:,indMovedAtom) = new3Coord(:,:)

            if (tBarostat) then ! apply a Barostat
              call rescale(pMDIntegrator,coord0,latVec,totalStress)
              recVec2p = latVec(:,:)
              call matinv(recVec2p)
              recVec2p = reshape(recVec2p, (/3, 3/), order=(/2, 1/))
              recVec = 2.0_dp * pi * recVec2p
              call getCellTranslations(cellVec, rCellVec, latVec, recVec2p, &
                  & mCutoff)
            end if
            if (tWriteRestart) then
              if (tStress) then
                if (tBarostat) then
                  write(fdMD,*)'Lattice vectors (A)'
                  do ii = 1, 3
                    write(fdMD,*)latVec(:,ii)*Bohr__AA
                  end do
                  write(fdMD,format1Ue)'Volume',cellVol,'au^3'
                end if
                write(fdMD,format2Ue)'Pressure',cellPressure,'au',&
                    & cellPressure * au__pascal, 'Pa'
                if (pressure/=0.0_dp) then
                  write(fdMD,format2U) 'Gibbs free energy', energy%Etotal -&
                      & sum(TS(:))+cellVol*pressure,'H', Hartree__eV * &
                      & (energy%Etotal - sum(TS(:))+cellVol*pressure),'eV'
                end if
              end if
              write(fdMD,format2U)'Potential Energy', &
                  & energy%Etotal - sum(TS),'H', &
                  & (energy%Etotal - sum(TS))*Hartree__eV,'eV'
              write(fdMD,format2U)'MD Kinetic Energy',KE,'H',KE*Hartree__eV,'eV'
              write(fdMD,format2U)'Total MD Energy', &
                  & KE+energy%Etotal - sum(TS),'H', &
                  & (KE+energy%Etotal - sum(TS))*Hartree__eV,'eV'
              write(fdMD,format2U)'MD Temperature', kT,'au',kT/Boltzmann,'K'
              if (tEfield) then
                write(fdMD,format1U1e)'External E field', absEField, 'au', &
                    &absEField * au__V_m, 'V/m'
              end if
              if (tDipole) then
                write(fdMD,"(' Dipole moment  :',3f14.8,' au')")dipoleMoment
                write(fdMD,"(' Dipole moment  :',3f14.8,' Debye')") &
                    & dipoleMoment*au__Debye
              end if
            end if
          end if
        end if
      end if
      
      if (tWriteDetailedOut.and.tMD) then
        write(fdUser,format1U)"MD Kinetic Energy", KE,"H"
        write(fdUser,format1U)"Total MD Energy", &
            & KE + energy%Etotal - sum(TS),"H"
        write(fdUser,format2U)"MD Temperature",kT,"H",kT/Boltzmann,"K"
        write(fdUser,*) ""
      end if
    end if

    if (tWriteDetailedOut .and. tEfield) then
      write(fdUser,format1U1e)'External E field', absEField,'au', &
          & absEField * au__V_m, 'V/m'
    end if
    
    !! Stop reading of initial charges/block populations again
    tReadChrg = .false.

    !! Stop SCC if appropriate stop file is present
    if (.not. tStopSCC) then
      inquire(file=fStopDriver, exist=tStopDriver)
      if (tStopDriver) then
        write (*,*) "Stop file '" // fStopDriver // "' found."
      end if
    end if
    if (tStopSCC .or. tStopDriver) then
      nGeoSteps = iGeoStep
      write (*,*) "Setting max number of geometry steps to current step number."
    end if
    iGeoStep = iGeoStep + 1

  end do lpGeomOpt

  if (tWriteDetailedOut.and.tDipole) then
    write(fdUser,"(' Dipole moment  :',3f14.8,' au')")dipoleMoment
    write(fdUser,"(' Dipole moment  :',3f14.8,' Debye')") &
        & dipoleMoment*au__Debye
    write(fdUser,*)''
  end if
  
  tGeomEnd = tMD .or. tGeomEnd .or. tDerivs

  if (tWriteDetailedOut) then
    if (tGeoOpt) then
      if (tGeomEnd) then
        write (fdUser,*) "Geometry converged"
      else
        write (fdUser, *) "!!! Geometry did NOT converge!"
      end if
    elseif (tMD) then
      if (tGeomEnd) then
        write (fdUser,*) "Molecular dynamics completed"
      else
        write (fdUser, *) "!!! Molecular dynamics terminated abnormally!"
      end if
    elseif (tDerivs) then
      if (tGeomEnd) then
        write (fdUser,*) "Second derivatives completed"
      else
        write (fdUser, *) "!!! Second derivatives terminated abnormally!"
      end if
    end if
    write(fdUser,*)''
    close(fdUser)
  end if

  if (tGeoOpt) then
    if (tGeomEnd) then
      write (*,*)
      write (*,*) "Geometry converged"
    else
      call warning("!!! Geometry did NOT converge!")
    end if
  elseif (tMD) then
    if (tGeomEnd) then
      write (*,*)
      write (*,*) "Molecular dynamics completed"
    else
      call warning("!!! Molecular dynamics terminated abnormally!")
    end if
  elseif (tDerivs) then
    if (tGeomEnd) then
      write (*,*)
      write (*,*) "Second derivatives completed"
    else
      call warning("!!! Second derivatives terminated abnormally!")
    end if
  end if

  if (tMD) then
    write(*,*)'MD information accumulated in ',mdOut
    close(fdMD)
  end if

  if (tDerivs) then
    call getHessianMatrix(pDerivDriver,pDynMatrix)
    write(*,*)'Hessian matrix written to ',hessianOut
    do ii = 1, size(pDynMatrix,dim=2)
      write(fdHessian,formatHessian)pDynMatrix(:,ii)
    end do
    close(fdHessian)
  end if

  if (tWriteTagged) then
    !! Write out final results in tagged format
    open(fdTagged, file=taggedOut, position="append")
    call writeTagged(fdTagged, tag_eigenVal, eigen)
    call writeTagged(fdTagged, tag_filling, filling)
    call writeTagged(fdTagged, tag_egyBand, Eband(:mod(nSpin,3)))
    call writeTagged(fdTagged, tag_egyBandT0, E0(:mod(nSpin,3)))
    call writeTagged(fdTagged, tag_egyRep, energy%Erep)
    if (tAtomicEnergy) then
      call writeTagged(fdTagged, tag_egyRepAt, energy%atomRep)
    end if
    if (tMulliken) then
      call qm2ud(qOutput)
      call writeTagged(fdTagged, tag_qOutput,qOutput(:,:,1))
      call writeTagged(fdTagged, tag_qOutputAt, sum(qOutput,dim=1))
      call ud2qm(qOutput)
    end if
    if (tForces) then
      call writeTagged(fdTagged, tag_forceTot, &
          &-repulsiveDerivs - derivs)
      call writeTagged(fdTagged, tag_forceBand, -derivs)
      call writeTagged(fdTagged, tag_forceRep, -repulsiveDerivs)
      if (tExtChrg) then
        call writeTagged(fdTagged, tag_chrgForces, -chrgForces)
      end if
      if (tStress) then
        call writeTagged(fdTagged, tag_stressRep, repulsiveStress)
        call writeTagged(fdTagged, tag_stressElec, elecStress)
        if (pressure/=0.0_dp) then
          call writeTagged(fdTagged, tag_pV,derivCellVol)
        end if
        if (tMD) then
          call writeTagged(fdTagged, tag_stressVirial, cellVolStress)
        end if
        call writeTagged(fdTagged, tag_stressTot, totalStress)
      end if
    end if
    if (tDerivs) then
      call writeTagged(fdTagged, tag_HessianNum, pDynMatrix)
    end if
    if (tSCC) then
      call writeTagged(fdTagged, tag_egySCC, energy%ESCC)
      if (tAtomicEnergy) then
        call writeTagged(fdTagged, tag_egySCCAt, energy%atomSCC)
      end if
    end if
    if (tEfield) then
      call writeTagged(fdTagged, tag_egyExt, energy%Eext)
      if (tAtomicEnergy) then
        call writeTagged(fdTagged, tag_egyExtAt, energy%atomExt)
      end if
    end if
    if (tDispersion) then
      call writeTagged(fdTagged, tag_egyDispersn ,energy%eDisp )
      if (tAtomicEnergy) then
        call writeTagged(fdTagged, tag_egyDispAt ,energy%atomDisp)
      end if
    end if
    if (tSpin) then
      call writeTagged(fdTagged, tag_egySpin, energy%Espin)
      if (tAtomicEnergy) then
        call writeTagged(fdTagged, tag_egySpinAt, energy%atomSpin)
      end if
    end if
    if (tDFTBU) then
      call writeTagged(fdTagged, tag_egyDFTBU, energy%Edftbu)
      if (tAtomicEnergy) then
        call writeTagged(fdTagged, tag_egyDFTBUAt, energy%atomDftbu)
      end if
    end if
    if (tSpinOrbit) then
      call writeTagged(fdTagged, tag_egyLS, energy%ELS)
      if (tAtomicEnergy) then
        call writeTagged(fdTagged, tag_egyLSAt, energy%atomDftbu)
      end if
    end if
    call writeTagged(fdTagged, tag_egyTotElec, energy%Eelec)
    call writeTagged(fdTagged, tag_egyTotal, energy%Etotal)
    if (tAtomicEnergy) then
      call writeTagged(fdTagged, tag_egyTotElAt, energy%atomElec)
      call writeTagged(fdTagged, tag_egyTotalAt, energy%atomTotal)
    end if
    call writeTagged(fdTagged, tag_tempElec, tempElec)
    call writeTagged(fdTagged, tag_entropy, sum(TS))
    call writeTagged(fdTagged, tag_freeEgy, energy%Etotal - sum(TS))
    if (pressure/=0.0_dp) then
      call writeTagged(fdTagged, tag_Gibbsfree, &
          & energy%Etotal - sum(TS) + cellVol*pressure)
    end if
    call writeTagged(fdTagged, tag_endCoord, coord0)

    close(fdTagged)
  end if

  if (tWriteResultsTag) then
    fdResultsTag = getFileId()
    call initTaggedWriter()
    open(fdResultsTag, file=resultsTag, position="rewind", status="replace")
    call writeTagged(fdResultsTag, tag_egyTotal, energy%Etotal)
    if (tAtomicEnergy) then
      call writeTagged(fdResultsTag, tag_egyTotalAt, energy%atomTotal)
    end if
    call writeTagged(fdResultsTag, tag_forces, tForces)
    if (tForces) then
      totalDeriv(:,:) = repulsiveDerivs(:,:) + derivs(:,:)
      if (tDispersion) then
        call addGradients(myDispersion, totalDeriv)
      end if
      call writeTagged(fdResultsTag, tag_forceTot, &
          &-totalDeriv)
      if (tExtChrg) then
        call writeTagged(fdResultsTag, tag_chrgForces, -chrgForces)
      end if
    end if
    if (tStress) then
      call writeTagged(fdResultsTag, tag_stressTot, totalStress)
    end if
    if (tDerivs) then
      call writeTagged(fdResultsTag, tag_HessianNum, pDynMatrix)
    end if
    call writeTagged(fdResultsTag, tag_scc, tSCC)
    if (tSCC) then
      call writeTagged(fdResultsTag, tag_nSCC, iSCCIter)
      call writeTagged(fdResultsTag, tag_sccConv, tConverged)
    end if
    if (tMulliken) then
      call writeTagged(fdResultsTag, tag_qOutputAt, sum(qOutput, dim=1))
      call writeTagged(fdResultsTag, tag_qOutAtNet, &
          &sum(q0(:,:,1) - qOutput(:,:,1), dim=1))
    end if
    call writeTagged(fdResultsTag, tag_eigenVal, eigen)
    call writeTagged(fdResultsTag, tag_filling, filling)
    call writeTagged(fdResultsTag, tag_efermi, Ef)
    if (size(nEl) == 1) then
      call writeTagged(fdResultsTag, tag_nElUp, 0.5_dp*nEl(1))
      call writeTagged(fdResultsTag, tag_nElDown, 0.5_dp*nEl(1))
    else
      call writeTagged(fdResultsTag, tag_nElUp, nEl(1))
      call writeTagged(fdResultsTag, tag_nElDown, nEl(2))
    end if
    close(fdResultsTag)
  end if


  if (tWriteDetailedXML) then
    !! Ugly hack for printing out xml info, will be removed later
    call xml_OpenFile("detailed.xml", xf, indent=.true.)
    call xml_ADDXMLDeclaration(xf)
    call xml_NewElement(xf, "detailedout")
    call writeChildValue(xf, "identity", runId)
    call xml_NewElement(xf, "geometry")
    call writeChildValue(xf, "typenames", specieName)
    call writeChildValue(xf, "typesandcoordinates", &
        &reshape(specie0, (/ 1, size(specie0) /)), pCoord0Out)
    call writeChildValue(xf, "periodic", tPeriodic)
    if (tPeriodic) then
      call writeChildValue(xf, "latticevectors", latVec)
    end if
    call xml_EndElement(xf, "geometry")
    call writeChildValue(xf, "real", tRealHS)
    call writeChildValue(xf, "nrofkpoints", nKPoint)
    call writeChildValue(xf, "nrofspins", nSpin)
    call writeChildValue(xf, "nrofstates", size(eigen, dim=1))
    call writeChildValue(xf, "nroforbitals", nOrb)
    ALLOCATE_(bufferRealR2, (4, nKPoint))
    bufferRealR2(1:3, :) = kPoint(:,:)
    bufferRealR2(4, :) = kWeight(:)
    call writeChildValue(xf, "kpointsandweights", bufferRealR2)
    DEALLOCATE_(bufferRealR2)
    call xml_NewElement(xf, "occupations")
    do ii = 1, nSpin
      call xml_NewElement(xf, "spin" // i2c(ii))
      do jj = 1, nKpoint
        call writeChildValue(xf, "k" // i2c(jj), filling(:, jj, mod(ii,3)))
      end do
      call xml_EndElement(xf, "spin" // i2c(ii))
    end do
    call xml_EndElement(xf, "occupations")
    call xml_EndElement(xf, "detailedout")
    call xml_Close(xf)
  end if

  DEALLOCATE_(orbitalL)
  DEALLOCATE_(orbitalLPart)
  DEALLOCATE_(new3Coord)
  DEALLOCATE_(tmpDerivs)

  !! Deallocate arrays
  DEALLOCATE_(filling)

  if (tForces) then
    DEALLOCATE_(derivs)
    DEALLOCATE_(repulsiveDerivs)
    DEALLOCATE_(chrgForces)
  end if

  call destroy(energy)
  call destroy(potential)
  call destroy(pConnect)

  if (tMulliken) then
    DEALLOCATE_(qOutput)
    DEALLOCATE_(qInput)
    DEALLOCATE_(q0)
  end if
  DEALLOCATE_(rhoPrim)
  DEALLOCATE_(iRhoPrim)
  if (tForces) then
    DEALLOCATE_(ERhoPrim)
  end if

  if (tMD) then
    DEALLOCATE_(velocities)
    DEALLOCATE_(movedVelo)
    DEALLOCATE_(movedAccel)
    DEALLOCATE_(movedMass)
  end if

  DEALLOCATE_(HSqrCplx)
  DEALLOCATE_(SSqrCplx)
  DEALLOCATE_(HSqrReal)
  DEALLOCATE_(SSqrReal)
  DEALLOCATE_(eigen)

  if (tSCC) then
    call destruct_SCC()
    call destroy(pChrgMixer)
  end if

  call destroyProgramVariables()

contains

  ! Invokes the writing routines for the Hamiltonian and overlap matrices.
  subroutine writeHS(tWriteHS, tWriteRealHS, ham, over, iNeighbor, &
          &nNeighbor, iAtomStart, iPair, img2CentCell, kPoint, iCellVec, &
          &cellVec)
    logical, intent(in) :: tWriteHS, tWriteRealHS
    real(dp), intent(in) :: ham(:,:), over(:)
    integer, intent(in) :: iNeighbor(0:,:), nNeighbor(:)
    integer, intent(in) :: iAtomStart(:), iPair(0:,:), img2CentCell(:)
    real(dp), intent(in) :: kPoint(:,:)
    integer, intent(in) :: iCellVec(:)
    real(dp), intent(in) :: cellVec(:,:)

    integer :: iS, nSpin

    nSpin = size(ham, dim=2)

    if (tWriteRealHS) then
      do iS = 1, nSpin
        call writeSparse("hamreal" // i2c(iS) // ".dat", ham(:,iS), iNeighbor, &
            &nNeighbor, iAtomStart, iPair, img2CentCell, iCellVec, cellVec)
      end do
      call writeSparse("overreal.dat", over, iNeighbor, &
          &nNeighbor, iAtomStart, iPair, img2CentCell, iCellVec, cellVec)
    end if
    if (tWriteHS) then
      if (tRealHS) then
        do iS = 1, nSpin
          call writeSparseAsSquare("hamsqr" // i2c(iS) // ".dat", ham(:,iS), &
              &iNeighbor, nNeighbor, iAtomStart, iPair, img2CentCell)
        end do
        call writeSparseAsSquare("oversqr.dat", over, iNeighbor, nNeighbor, &
            &iAtomStart, iPair, img2CentCell)
      else
        do iS = 1, nSpin
          call writeSparseAsSquare("hamsqr" // i2c(iS) // ".dat", ham(:,iS), &
              &kPoint, iNeighbor, nNeighbor, iAtomStart, iPair, img2CentCell, &
              &iCellVec, cellVec)
        end do
        call writeSparseAsSquare("oversqr.dat", over, kPoint, iNeighbor, &
            &nNeighbor, iAtomStart, iPair, img2CentCell, iCellVec, cellVec)
      end if
    end if

  end subroutine writeHS


end program dftb
