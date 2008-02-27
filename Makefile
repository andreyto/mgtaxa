CC := $(shell which gcc)
CFLAGS := -O0
CXX := $(shell which g++)
#To avoid "undefined symbols" message from gdb print command
#even when all optimization is off,
#use -gdwarf-2 or -gstabs compiler option that sets debug info
#format.
#-ggdb is essential for STL data structures to show in debugger
#and for STL gdb macros to work (from ~/.gdb/stl_views).
#It's a little odd that we need to tell GCC to produce gdb
#specific debugging data on a linux system...
CXXFLAGS := -g -gstabs -ggdb -O0
MAKE := make
AR := $(shell which ar)
ARFLAGS := -rvs
DOXYGEN := $(shell which doxygen)
MAKEDEPEND = mkdir -p $(DEP_DIR); $(CXX) $(CXXFLAGS) $(EXTRA_CXXFLAGS) -MM -o $(DEP_DIR)/$*.d $<

PROGRAMS  := test_kmers
LIBRARIES := libMGT.a
EXTRA_CXXFLAGS = -I$(PROJ_DIR)/include

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
## source tree, as long as they are not involved in target compilation, and make
## will still succeed.

PROJ_DIR := $(HOME)/work/mgtaxa
BUILD_DIR := $(PROJ_DIR)/build/$(MACH)
SRC_DIR := $(PROJ_DIR)/src
PY_DIR := $(PROJ_DIR)/MGT
INC_DIR := $(PROJ_DIR)/include
TEST_DIR := $(PROJ_DIR)/test
TEST_SRC_DIR := $(TEST_DIR)/cpp
TEST_PY_DIR := $(TEST_DIR)/py
DOC_DIR := $(PROJ_DIR)/doc

DEP_DIR := $(BUILD_DIR)/.deps

SRC_DIRS := $(SRC_DIR) $(TEST_SRC_DIR)
PY_DIRS := $(PY_DIR) $(TEST_PY_DIR)

vpath %.h $(INC_DIR)
vpath %.hpp $(INC_DIR)
vpath %.cpp $(SRC_DIRS)
vpath %.py $(PY_DIRS)

####################### .PHONY Target Definitions #########################

.PHONY: all
all: $(PROGRAMS) $(LIBRARIES) doc

.PHONY: doc
doc: $(DOC_DIR)/html

$(DOC_DIR)/html: $(SRC) $(PY) $(DOC_DIR)/Doxyfile
	@echo && echo "Running Doxygen" && echo
	cd $(PROJ_DIR) && $(DOXYGEN) $(DOC_DIR)/Doxyfile > /dev/null

.PHONY: clean
clean:		
	rm -f $(PROGRAMS) $(LIBRARIES) *.o *.pyc *.pyo $(DEP_DIR)/*.P
	rm -rf $(DOC_DIR)/html $(DOC_DIR)/tex


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

# Define macro for library building
define BUILD_LIB
$(AR) $(ARFLAGS) $@ $?
endef

#####################################################################################
#### Define each target, adding custom libs and options after LINK_EXE as needed ####

libMGT.a: kmers.o
	$(BUILD_LIB)

test_kmers: test_kmers.o libMGT.a
	$(LINK_EXE)

########################## End target definitions ###################################
#####################################################################################

#################### Include generated dependency files as Makefiles #######

#We always create an empty dependency file (it needs to have at least a new line in it)
#so that -include would not print a warning after make is run after previous 'make clean'

$(shell echo '' > $(DEP_DIR)/dummy_make_empty.P)
-include $(wildcard $(DEP_DIR)/*.P)
