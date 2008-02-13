CC := $(shell which gcc)
CFLAGS = -O3
CXX := $(shell which g++)
CXXFLAGS=
MAKE=make
AR=/usr/bin/ar
ARFLAGS=-rvs
DOXYGEN=$(shell which doxygen)
MAKEDEPEND = mkdir -p $(DEP_DIR); $(CXX) $(CXXFLAGS) $(EXTRA_CXXFLAGS) -MM -o $(DEP_DIR)/$*.d $<
DEP_DIR = .deps
df = $(DEP_DIR)/$(*F)

## Dependency autogeneration code is taken from:
##http://make.paulandlesley.org/autodep.html

PROJ_DIR = $(HOME)/work/mgtaxa
SRC_DIR = $(PROJ_DIR)/src
PY_DIR = $(PROJ_DIR)/MGT
INC_DIR = $(PROJ_DIR)/include
TEST_DIR = $(PROJ_DIR)/test
TEST_SRC_DIR = $(TEST_DIR)/cpp
TEST_PY_DIR = $(TEST_DIR)/py
DOC_DIR = $(PROJ_DIR)/doc

PROGRAMS   = test_kmers

EXTRA_CXXFLAGS = -I$(PROJ_DIR)/include

#EXTRA_LIBS = -L/usr/X11R6/lib -L/usr/local/lib -lX11 -lgzstream -lz -lm

#EXTRA_LIBS = -L/usr/X11R6/lib -L/home/syooseph/gzstream_bioinfo/gzstream -lX11 -lgzstream -lz -lm

#VPATH = $(SRC_DIR):$(TEST_DIR)
vpath %.h $(INC_DIR)
vpath %.hpp $(INC_DIR)
vpath %.cpp $(SRC_DIR) $(TEST_SRC_DIR)
vpath $.py $(PY_DIR) $(TEST_PY_DIR)

# Source files
SRC_CPP = $(wildcard *.cpp)
SRC_C = $(wildcard *.c)
SRC = $(SRC_CPP) $(SRC_C)
PY = $(wildcard *.py)
#DEP = $(SRC_CPP:.cpp=.d) $(SRC_C:.c=.d)
OBJ = $(SRC:.cpp=.o) $(SRC:.c=.o)

######################### Target Definitions ##########################

# Targets that you shouldn't need to change

all: $(PROGRAMS) $(LIBRARIES) doc

doc: $(SRC) $(PY) $(DOC_DIR)/Doxyfile
	cd $(PROJ_DIR); $(DOXYGEN) $(DOC_DIR)/Doxyfile

clean:		
	rm -f $(PROGRAMS) $(LIBRARIES) *.o *.pyc *.pyo $(DEP_DIR)/*.P
	rm -rf $(DOC_DIR)/html $(DOC_DIR)/tex

install: all
ifdef PROGRAMS
	for PR in $(PROGRAMS) ; do cp $$PR $(PREFIX)/bin ; done
endif
ifdef LIBRARIES
	for LB in $(LIBRARIES) ; do cp $$LB $(PREFIX)/lib ; done
endif

%.o: %.cpp
	@$(MAKEDEPEND); \
            cp $(df).d $(df).P; \
            sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
                -e '/^$$/ d' -e 's/$$/ :/' < $(df).d >> $(df).P; \
            rm -f $(df).d
	$(CXX) $(CXXFLAGS) $(EXTRA_CXXFLAGS) -c -o $*.o $<

%.o: %.c
	@$(MAKEDEPEND); \
            cp $(df).d $(df).P; \
            sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
                -e '/^$$/ d' -e 's/$$/ :/' < $(df).d >> $(df).P; \
            rm -f $(df).d
	$(CC) $(CFLAGS) $(EXTRA_CFLAGS) -c -o $*.o $<

#%.d: %.c
#	$(CC) -M $(EXTRA_FLAGS) $< > $@

#%.d: %.cpp
#	$(CXX) -M $(EXTRA_FLAGS) $< > $@

# Add .d to Make's recognized suffixes.
#SUFFIXES += .d

# Generate file.d from file.cpp.
#%.d: %.cpp
#	$(CXX) $(CXXFLAGS) $(EXTRA_CXXFLAGS) -MM $< > $@

#%.d: %.cpp
#	$(SHELL) -ec '$(CXX) -MM $(CXXFLAGS) $(EXTRA_CXXFLAGS) $(CPPFLAGS) $< \
#                      | sed '\''s/\($*\)\.o[ :]*/\1.o $@ : /g'\'' > $@; \
#                      [ -s $@ ] || rm -f $@'

# In order to support multiple targets, you need to place a target
# for each program and library down here.

# Sample program target
#
# myprog: $(MYPROG_OBJ) $(EXTRA_LIBS)
#	$(CC) -o myprog $(CFLAGS) $(EXTRA_FLAGS) $(MYPROG_OBJ) $(EXTRA_LIBS)

# Sample library target
#
# mylib.a: $(MYLIB_OBJ)
#       $(AR) $(ARFLAGS) mylib.a $(MYLIB_OBJ)

test_kmers: test_kmers.o
	$(CXX) -o test_kmers $(CXXFLAGS) $(EXTRA_CXXFLAGS) test_kmers.o \
	$(EXTRA_LIBS) -lm


##
## Dependencies
##

#C_DEPS   = $(patsubst %.c, %.d, $(filter %.c, $(SRC)))
#CXX_DEPS = $(patsubst %.cpp, %.d, $(filter %.cpp, $(SRC)))

-include $(SRC:%.c=$(DEP_DIR)/%.P)
-include $(SRC:%.cpp=$(DEP_DIR)/%.P)

#sinclude $(CXX_DEPS) $(C_DEPS)
#ifeq (0, $(words $(findstring $(MAKECMDGOALS), clean distclean)))
#include $(DEP)
#endif
