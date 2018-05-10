# Build and test main dftb application.
############################################################################
#  Global variables
############################################################################

# Source root directory (containing Makefile.common, Makefile.objdir, etc.)
ROOT = ..

# target executables (have to be defined, if you want to use make (un)install)
TARGETS = dftb+

# Targets which do not need a source dependency tree
NODEPTARGETS = compare summary tar testclean


############################################################################
#  Put objects in separate dir
############################################################################

ifeq (,$(filter _obj%,$(notdir $(CURDIR))))
include $(ROOT)/Makefile.objdir
else

############################################################################
# Usual Makefile (executed in the object directory)
############################################################################

# Default target (declared *before* importing global makefiles).
.PHONY: all
all: $(TARGETS)

# Default autotest directory containing the test directories and the file
# describing which tests to carry out. (You can override them in Makefile.user) 
PRGDFTB_TESTDIR = $(ROOT)/../autotest
PRGDFTB_TESTFILE = $(PRGDFTB_TESTDIR)/tests

# Include common makefile
include $(ROOT)/Makefile.common


# Target definitions (declared *after* importing necessary makefiles).
dftb+: $(dftb+.o)
	$(link-target)


#---------------------------------------------------------------------
# Autotest
#---------------------------------------------------------------------
#
# NB:	- AUTOTEST_WORKDIR is now created by autotest2
#	- User may, but probably should not, override AUTOTEST_FLAGS.

# Program to use for running and comparing
AUTOTEST_SCRIPT= $(PRGDFTB_TESTDIR)/bin/autotest2
TAGDIFF_SCRIPT = $(PRGDFTB_TESTDIR)/bin/tagdiff

# Complete command with options
AUTOTEST_CMD = $(AUTOTEST_SCRIPT) \
		-r $(PRGDFTB_TESTDIR) \
		-w $(AUTOTEST_WORKDIR) \
		-d $(TAGDIFF_SCRIPT)

# Default autotest working directory (you can override it in Makefile.user)
# (if relative path: relative to the _obj_* directory)
AUTOTEST_WORKDIR = _autotest

# Variables for conveniently archiving (and presumably transferring)
# autotest results.
#
# Caveat:
#	  Careful when unpacking this after a transfer, as local
#	results from the same ARCH may be overwritten.
#	Appending HOSTNAME to OBJDIR is more hassle, since then both
#	ARCH and HOSTNAME would need overriding to do "make compare".

TAR_SRC = $(notdir $(CURDIR))/$(strip $(AUTOTEST_WORKDIR))
TAR_DST = $(HOME)/autotest-$(ARCH)

# Get list of tests to execute, but only if test is specified as goal
TESTS := $(if $(filter test summary compare,$(MAKECMDGOALS)),\
	$(shell $(AUTOTEST_CMD) -f $(PRGDFTB_TESTFILE) -l),)


.PHONY: $(TESTS) run test compare summary tar testclean

# Prepare and run individual tests, possibly in parallel.
$(TESTS): dftb+
	@$(AUTOTEST_CMD) -p ./dftb+ -s P,R $(AUTOTEST_FLAGS) $@

# Run calculations only (perhaps python is missing to run tagdiff)
run: AUTOTEST_FLAGS = -v
run: $(TESTS)

# Run calculations and immediately compare (possibly in parallel),
# then summarize.  The dependency ensures proper sequencing.
test: AUTOTEST_FLAGS = -s C -v
test: $(TESTS)
	@$(AUTOTEST_CMD) -s S $(TESTS)

# Targets to run analyses post facto.
compare:
	@$(AUTOTEST_CMD) -s C $(AUTOTEST_FLAGS) $(TESTS)

summary:
	-@$(AUTOTEST_CMD) -s S $(AUTOTEST_FLAGS) $(TESTS)

# NB:
#    - Initial cd goes to parent of OBJDIR, which is not necessarily SRCDIR.)
#    - Cannot assume tar == gtar.
#
# Idea:  For support, copy system and user config info (makefiles etc.)
#	into AUTOTEST_WORKDIR.
tar:
	@( cd ..; echo "# tarring in `pwd` ..."; \
	    tar cvf - $(TAR_SRC) | gzip > $(TAR_DST).tgz )
	@echo $(TAR_DST).tgz

testclean:
	rm -rf $(AUTOTEST_WORKDIR)

############################################################################
endif
