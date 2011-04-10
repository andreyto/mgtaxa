### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##

## Installation directories.

## We support a "home" type installation,
## where everything for our application is installed under
## application's own subdirectory, except for the 'datadir'.
## (or possibly two - one for arch dependent and another for arch independent files).
## This makes it easier to create test installs as opposed to
## "prefix" type installation (where executables for all applications go
## into a single common prefix/bin etc). The installation variables
## are still called 'prefix' and such, and are low case.

# Software will be installed here
prefix := $(INST)/mgtaxa

# Data files will be stored here
datadir := $(prefix)/var

# Temporary files will be stored here
tmpdir := $(prefix)/tmp

exec_prefix := $(INSTMACH)/mgtaxa

## Scripts go here
bindir := $(prefix)/bin
## Binary executables go here
exec_bindir := $(exec_prefix)/bin
## Binary libraries go here
libdir := $(exec_prefix)/lib
## profile to source goes here
sysconfdir := $(prefix)/etc
docdir := $(prefix)/doc

CC := $(shell which gcc)
CFLAGS := -O3 #-O0
CXX := $(shell which g++)
#To avoid "undefined symbols" message from gdb print command
#even when all optimization is off,
#use -gdwarf-2 or -gstabs compiler option that sets debug info
#format.
#-ggdb is essential for STL data structures to show in debugger
#and for STL gdb macros to work (from ~/.gdb/stl_views).
#It's a little odd that we need to tell GCC to produce gdb
#specific debugging data on a linux system...
CXXFLAGS := -O3 -fPIC -g #-g -gstabs -ggdb -O0 -fPIC #-O3 -fPIC 
#CXXFLAGS := -g -O0 -fPIC #-g -gstabs -ggdb -O0 -fPIC #-O3 -fPIC 
PYTHON := $(shell which python)
MAKE := make
AR := $(shell which ar)
ARFLAGS := -rvs
DOXYGEN := $(shell which doxygen)
#CPR := cp -dR
CPR := $(shell which rsync) -av
MKDIRP := mkdir -p
RMRF := rm -rf
MAKEDEPEND = $(MKDIRP) $(DEP_DIR); $(CXX) $(CXXFLAGS) $(EXTRA_CXXFLAGS) -MM -o $(DEP_DIR)/$*.d $<

PROGRAMS       := fastaFilter 
PROGRAMS_TEST  := test_kmers
LIBRARIES      := libMGT.a
#This should list only extensions that are to be installed, and omit testing ones
PYEXT          := MGTX/sample_boost.so MGTX/kmersx.so
PYEXT_TEST     := test_sample_boostx.so test_numpy_boostx.so

EXTRA_CXXFLAGS = -I$(PROJ_DIR)/include -I$(BOOST_INC_DIR) -I$(PY_INC_DIR) -I$(NUMPY_INC_DIR)

## Installation subdirectories

# Python extensions will be here. $(libdir) must be added to PYTHONPATH,
# so that the extensions can be 'import MGTX.<extmodule>'.
# Somehow having two identical subdirectories 'MGT' under 'prefix' and
# 'exec_prefix' does not work (extension modules are never found)

extdir := $(libdir)/MGTX

## How we generate autodependencies:
## Dependency autogeneration code is taken from:
## http://make.paulandlesley.org/autodep.html
## with the following change: We include here only those .P files
## that already exist (through wildcard search). This way, we do not
## need to define (or build through wildcard search) a list of all source files.
## What happens:
## Whenever a specific source file has to be compiled, an attempt is made to
## generate a dependency .P file for it. If preprocessor fails, the compilation
## of source file fails (as it should). If dependency file is generated, it will
## be included in this Makefile the next time make is run, in which case make
## will cause a recompilation of the source on any dependency update.
## If the dependency file does not exist (and not included in Makefile) then
## .o file does not exist either, in which case source compilation is done anyway.
## Compared to generating a list of all source files through wildcard search,
## this approach has the benefit that we can have any malformed source files in our
## source tree, as long as they are not involved in target compilation, and 'make'
## will still succeed.

PROJ_DIR := $(shell cd ..; pwd)
BUILD_DIR := $(shell pwd)
BUILD_EXT_DIR := $(BUILD_DIR)/MGTX
SRC_DIR := $(PROJ_DIR)/src
EXT_DIR := $(PROJ_DIR)/ext
PY_DIR := $(PROJ_DIR)/MGT
EXAMPLE_DATA_DIR := $(PROJ_DIR)/example_data
INC_DIR := $(PROJ_DIR)/include
TEST_DIR := $(PROJ_DIR)/test
TEST_SRC_DIR := $(TEST_DIR)/cpp
TEST_PY_DIR := $(TEST_DIR)/py
DOC_DIR := $(PROJ_DIR)/doc
BIN_PY_DIR := $(PROJ_DIR)/bin

DEP_DIR := $(BUILD_DIR)/.deps

SRC_DIRS := $(SRC_DIR) $(EXT_DIR) $(TEST_SRC_DIR)
PY_DIRS := $(PY_DIR) $(TEST_PY_DIR)

vpath %.h $(INC_DIR)
vpath %.hpp $(INC_DIR)
vpath %.cpp $(SRC_DIRS)
vpath %.py $(PY_DIRS)

####################### Support for building Python extensions ############

PY_INC_DIR := $(shell $(PYTHON) -c 'from distutils.sysconfig import *; print get_python_inc()')
## As per distutils documentation,
## get_config_var() might not be portable.
PY_LIB_DIR := $(shell $(PYTHON) -c 'from distutils.sysconfig import *; print get_config_var("LIBPL")')
PY_LIB := $(shell $(PYTHON) -c 'from distutils.sysconfig import *; print get_config_var("LIBRARY")')

NUMPY_INC_DIR := $(shell $(PYTHON) -c 'import numpy; print numpy.get_include()')

#We build a single Python extension with Boost interface library, and therefore
#use the static liboost_python.a. This saves a lot of grief with
#finicky -Wl<xx> linker options, especially if those are already used
#by custom gcc installation. As far as I remember, cross-module support in
#Boost.python would require linking with a shared boost library.

ifeq ($(BOOST_OS),YES)
ifneq (,$(filter x86_32%,$(MACH)))
BOOST_LIB_DIR := /usr/lib
else 
ifneq (,$(filter x86_64%,$(MACH)))
BOOST_LIB_DIR := /usr/lib64
else
$(error "Unknown MACH variable value: $(MACH))
endif
BOOST_INC_DIR := /usr/include
BOOST_PY_ST_LIB := libboost_python.a
BOOST_PY_SH_LIB := boost_python
endif
else
BOOST_INC_DIR := $(INST)/include/boost-1_37
BOOST_LIB_DIR := $(INSTMACH)/lib
BOOST_PY_ST_LIB := libboost_python-gcc41-mt-1_37.a
BOOST_PY_SH_LIB := boost_python-gcc41-mt
endif

ifdef ($BOOST_STATIC)
BOOST_PY_LINK := $(BOOST_LIB_DIR)/$(BOOST_PY_ST_LIB)
else
BOOST_PY_LINK := -L$(BOOST_LIB_DIR) -l$(BOOST_PY_SH_LIB)
endif

#Example of debugging the make process.
#Debugging using 'echo' outside of target definition works fine as 
#long as we redirect all output to a file.
#$(shell echo $(PY_INC_DIR) &> make.log)
#$(info $(PY_INC_DIR))

####################### .PHONY Target Definitions #########################

.PHONY: all
all: build doc

.PHONY: build
build: $(PROGRAMS) $(PROGRAMS_TEST) $(LIBRARIES) $(PYEXT) $(PYEXT_TEST) mgtaxa.shrc mgtaxa.insrc.shrc


.PHONY: doc
doc: $(DOC_DIR)/html

$(DOC_DIR)/html: $(SRC) $(PY) $(DOC_DIR)/Doxyfile
	@echo && echo "Running Doxygen" && echo
	cd $(PROJ_DIR) && $(DOXYGEN) $(DOC_DIR)/Doxyfile > /dev/null

## Testing is done "in place", before the "install" target is made.
## Still, we can run Python tests as though they call the installed
## application modules due to our directory structure. In particular,
## Python extensions are built into a subdirectory MGTX, and so can be
## "import MGTX.<extension_name>"

test: build
	rm -rf .testrun; mkdir .testrun
	echo "#!/bin/sh"
	echo "export PYTHONPATH=$(TEST_PY_DIR):$(BUILD_DIR):$(PROJ_DIR):$(PYTHONPATH)" >> .testrun/run
	echo "../test_kmers" >> .testrun/run
	echo "$(PYTHON) $(TEST_PY_DIR)/main.py" >> .testrun/run
	chmod +x .testrun/run
	cd .testrun && ./run

.PHONY: install
install: all
	install -d $(bindir)
	install -d $(exec_bindir)
	install -d $(extdir)
	install -d $(docdir)
	install -d $(sysconfdir)
	$(CPR) $(PY_DIR) $(prefix)
	$(CPR) $(EXAMPLE_DATA_DIR) $(prefix)
	install -t $(extdir) $(PYEXT) MGTX/__init__.py
	$(CPR) $(DOC_DIR)/html $(docdir)
	install -t $(sysconfdir) mgtaxa.shrc $(PROJ_DIR)/etc/gt.sketch.default.style
	install -t $(exec_bindir) $(PROGRAMS)
	$(CPR) $(BIN_PY_DIR)/* $(bindir)

.PHONY: clean
clean:		
	$(RMRF) $(PROGRAMS) $(PROGRAMS_TEST) $(LIBRARIES) $(BUILD_EXT_DIR) *.o *.so *.pyc *.pyo $(DEP_DIR)/*.P
	$(RMRF) $(DOC_DIR)/html $(DOC_DIR)/tex


mgtaxa.shrc: $(PROJ_DIR)/etc/mgtaxa.shrc.in
	sed -e 's|__MGT_HOME__|$(prefix)|' \
	    -e 's|__MGT_BIN__|$(bindir)|' \
	    -e 's|__MGT_DATA__|$(datadir)|' \
	    -e 's|__MGT_TMP__|$(tmpdir)|' \
	    -e 's|__MGT_EXEC_BIN__|$(exec_bindir)|' \
	    -e 's|__MGT_PY_PATH__|$(prefix):$(libdir)|' \
	    -e 's|__MGT_RC__|$(sysconfdir)/mgtaxa.shrc|' \
	    $(PROJ_DIR)/etc/mgtaxa.shrc.in > $@ || rm $@

mgtaxa.insrc.shrc: $(PROJ_DIR)/etc/mgtaxa.shrc.in
	sed -e 's|__MGT_HOME__|$(PROJ_DIR)|' \
	    -e 's|__MGT_BIN__|$(BUILD_DIR)|' \
	    -e 's|__MGT_DATA__|$(datadir)|' \
	    -e 's|__MGT_TMP__|$(tmpdir)|' \
	    -e 's|__MGT_EXEC_BIN__|$(BUILD_DIR)|' \
	    -e 's|__MGT_PY_PATH__|$(PROJ_DIR):$(BUILD_DIR)|' \
	    -e 's|__MGT_RC__|$(BUILD_DIR)/mgtaxa.insrc.shrc|' \
	    $(PROJ_DIR)/etc/mgtaxa.shrc.in > $@ || rm $@

############################ Compilation Rules #############################

##### A couple of macro definitions for autodependency ######
##### generation in compile rules:                     ######

# Dynamically expands into dependency file name:

define df 
$(DEP_DIR)/$(*F)
endef

# Create .P file from .d file:

define df_pf
cp $(df).d $(df).P && \
            sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
                -e '/^$$/ d' -e 's/$$/ :/' < $(df).d >> $(df).P && \
            rm $(df).d
endef

####### Compile here ########

%.o: %.cpp
	$(MAKEDEPEND) && $(df_pf)
	$(CXX) $(CXXFLAGS) $(EXTRA_CXXFLAGS) -c -o $*.o $<

%.o: %.c
	@$(MAKEDEPEND); $(df_pf)
	$(CC) $(CFLAGS) $(EXTRA_CFLAGS) -c -o $*.o $<


######################### Real Target Definitions ##########################

# Define macro for executable linking
define LINK_EXE
$(CXX) $(CXXFLAGS) $(EXTRA_CXXFLAGS) $(LDFLAGS) $(EXTRA_LIBS) -lm -o $@ $^
endef

# Define macro for static library building
define BUILD_LIB
$(AR) $(ARFLAGS) $@ $?
endef

# Define macro for shared library linking
define LINK_SO
$(CXX) $(CXXFLAGS) $(EXTRA_CXXFLAGS) $(LDFLAGS) $(EXTRA_LIBS) -shared -lm -o $@ $^
endef

# Define macro for Python extension linking with Boost.Python
define LINK_EXT
install -d $(BUILD_EXT_DIR)
touch $(BUILD_EXT_DIR)/__init__.py
$(LINK_SO) $(BOOST_PY_LINK)
endef


#####################################################################################
#### Define each target, adding custom libs and options after LINK_EXE as needed ####

libMGT.a: kmers.o py_num_util.o
	$(BUILD_LIB)

fastaFilter: fastaFilter.o
	$(LINK_EXE)

test_kmers: test_kmers.o libMGT.a
	$(LINK_EXE)

MGTX/sample_boost.so: sample_boost.o
	$(LINK_EXT)

MGTX/kmersx.so: kmersx.o libMGT.a
	$(LINK_EXT) -lz

test_sample_boostx.so: test_sample_boostx.o
	$(LINK_EXT)

test_numpy_boostx.so: test_numpy_boostx.o
	$(LINK_EXT)


########################## End target definitions ###################################
#####################################################################################

#################### Include generated dependency files as Makefiles #######

#We always create an empty dependency file (it needs to have at least a new line in it)
#so that -include would not print a warning after make is run after previous 'make clean'

$(shell echo '' > $(DEP_DIR)/dummy_make_empty.P)
-include $(wildcard $(DEP_DIR)/*.P)
