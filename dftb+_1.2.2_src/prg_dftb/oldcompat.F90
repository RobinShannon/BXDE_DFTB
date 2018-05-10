!!* Contains routines to convert HSD input for old parser to the current format.
module OldCompat
  use Message
  use flib_dom
  use HSDUtils
  use HSDUtils2
  use CharManip
  use XMLUtils
  implicit None
  private

  public :: convertOldHSD


contains

  !!* Converts an HSD input for an older parser to the current format.
  !!* @param root Root tag of the HSD-tree
  !!* @param oldVersion Version number of the old parser
  !!* @param curVersion Version number of the current parser
  subroutine convertOldHSD(root, oldVersion, curVersion)
    type(fnode), pointer :: root
    integer, intent(in) :: oldVersion
    integer, intent(in) :: curVersion

    integer :: version

    version = oldVersion
    do while (version < curVersion)
      select case(version)
      case(1)
        call convert_1_2(root)
        version = 2
      case(2)
        call convert_2_3(root)
        version = 3
      case (3)
        call convert_3_4(root)
        version = 4
      end select
    end do
    
  end subroutine convertOldHSD

  

  !!* Converts input from version 1 to 2. (Version 2 introcuded in August 2006)
  !!* @param root Root tag.
  subroutine convert_1_2(root)
    type(fnode), pointer :: root

    type(fnode), pointer :: child1, child2

    call getChild(root, "Geometry", child1, requested=.false.)
    if (associated(child1)) then
      call setUnprocessed(child1)
      call getChild(child1, "SpecieNames", child2, requested=.false.)
      if (associated(child2)) then
        call setUnprocessed(child2)
        call setNodeName(child2, "TypeNames")
      end if
    end if
    
  end subroutine convert_1_2


  !!* Converts input from version 2 to 3. (Version 3 introduced in Nov. 2006)
  subroutine convert_2_3(root)
    type(fnode), pointer :: root

    type(fnode), pointer :: ch1, ch2, par
    logical :: tValue

    call getDescendant(root, &
        &"Driver/VelocityVerlet/Thermostat/Andersen/RescalingProbability", &
        &ch1)
    if (associated(ch1)) then
      call detailedWarning(ch1, "Keyword renamed to 'ReselectProbability'.")
      call setNodeName(ch1, "ReselectProbability")
    end if

    call getDescendant(root, &
        &"Driver/VelocityVerlet/Thermostat/Andersen/RescaleIndividually", &
        &ch1)
    if (associated(ch1)) then
      call detailedWarning(ch1, "Keyword renamed to 'ReselectIndividually'.")
      call setNodeName(ch1, "ReselectIndividually")
    end if

    call getDescendant(root, "Hamiltonian/DFTB/Variational", ch1)
    if (associated(ch1)) then
      call getChildValue(ch1, "", tValue)
      call setUnprocessed(ch1)
      if (.not. tValue) then
        call detailedError(ch1, "Sorry, non-variational energy calculation &
            &is not supported any more!")
      else
        call detailedWarning(ch1, "Energy calculation is made only &
            &variational, option removed.")
        call destroyNode(ch1)
      end if
    end if
    
    call getDescendant(root, "Hamiltonian/DFTB/SCC", ch1, parent=par)
    if (associated(ch1)) then
      call getChildValue(ch1, "", tValue)
      call setUnprocessed(ch1)
      if (tValue) then
        call setChildValue(par, "OrbitalResolvedSCC", .true., child=ch2)
        call setUnprocessed(ch2)
        call detailedWarning(ch2, "Calculations are not orbital resolved &
            &per default any more. Keyword 'OrbitalResolvedSCC' added.")
      end if
    end if

    call getDescendant(root, "Options/PrintEigenvectors", ch1)
    if (associated(ch1)) then
      call detailedWarning(ch1, "Keyword converted to 'WriteEigenvectors'")
      call setNodeName(ch1, "WriteEigenvectors")
    end if
    
    call getDescendant(root, "Options/WriteTaggedOut", ch1)
    if (associated(ch1)) then
      call detailedWarning(ch1, "Keyword converted to 'WriteAutotestTag'. &
          &Output file name changed to 'autotest.out'")
      call setNodeName(ch1, "WriteAutotestTag")
    end if
    
    call getDescendant(root, "Options/WriteBandDat", ch1)
    if (associated(ch1)) then
      call detailedWarning(ch1, "Keyword converted to 'WriteBandOut'. &
          &Output file name changed to 'band.out'")
      call setNodeName(ch1, "WriteBandOut")
    end if
    
  end subroutine convert_2_3



  !!* Converts input from version 3 to 4. (Version 3 introduced in Mar. 2010)
  subroutine convert_3_4(root)
    type(fnode), pointer :: root

    type(fnode),pointer :: node, node2, node3
    type(fnodeList), pointer :: children
    integer :: ii

    !! Replace range operator with short start:end syntax
    call getDescendant(root, "Driver/SteepestDescent/MovedAtoms", node)
    call replaceRange(node)
    call getDescendant(root, "Driver/ConjugateGradient/MovedAtoms", node)
    call replaceRange(node)
    call getDescendant(root, "Driver/SecondDerivatives/Atoms", node)
    call replaceRange(node)
    call getDescendant(root, "Driver/VelocityVerlet/MovedAtoms", node)
    call replaceRange(node)
    call getDescendant(root, "Hamiltonian/DFTB/SpinPolarisation/Colinear&
        &/InitialSpin", node)
    if (associated(node)) then
      call getChildren(node, "AtomSpin", children)
      do ii = 1, getLength(children)
        call getItem1(children, ii, node2)
        call getChild(node2, "Atoms", node3)
        call replaceRange(node3)
      end do
      call destroyNodeList(children)
    end if
    
    call getDescendant(root, "Hamiltonian/DFTB/SpinPolarisation/Colinear&
        &/InitialSpin", node)
    if (associated(node)) then
      call detailedWarning(node, "Keyword renamed to 'InitalSpins'.")
      call setNodeName(node, "InitialSpins")
    end if


  contains

    ! Helper function
    subroutine replaceRange(node)
      type(fnode), pointer :: node

      type(fnode), pointer :: node2
      integer :: bounds(2)

      if (associated(node)) then
        call getChild(node, "Range", node2, requested=.false.)
        if (associated(node2)) then
          call getChildValue(node2, "", bounds)
          call removeChildNodes(node)
          call setChildValue(node, "", &
              &i2c(bounds(1)) // ":" // i2c(bounds(2)), replace=.true.)
          call detailedWarning(node, "Specification 'Range { start end }' &
              &not supported any more, using 'start:end' instead")
        end if
      end if
      
    end subroutine replaceRange
    
  end subroutine convert_3_4


  
end module OldCompat
