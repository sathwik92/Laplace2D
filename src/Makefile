######################################################################
#Include the variables from PETSC directory
######################################################################
include ${PETSC_DIR}/lib/petsc/conf/variables
######################################################################
# MODIFY PARAMETERS IN THIS SECTION
######################################################################
# Add the names of your .cpp files below
######################################################################
FILES = input_reader main utiltrig make_global apply_bc petsc_diagonalizer solution


######################################################################
# The name of the binary to be output
######################################################################
BINARY_NAME= Laplace2d
######################################################################

######################################################################
# The MPI enabled compiler to use
######################################################################
MPICXX= mpicxx
######################################################################
# GENERATE NAMES OF FILES
######################################################################
OBJ_FILES = $(addsuffix .o, $(FILES))

######################################################################
# COMPILATION FLAGS
######################################################################
CPPFLAGS = -I $(PETSC_DIR)/include/ -I $(PETSC_DIR)/$(PETSC_ARCH)/include/

######################################################################
# LINKER FLAGS
######################################################################
LDFLAGS = -L $(PETSC_DIR)/$(PETSC_ARCH)/lib/\
          -L $(SLEPC_DIR)/$(PETSC_ARCH)/lib/ \
          -Wl,-rpath,$(PETSC_DIR)/$(PETSC_ARCH)/lib\
          -lpetsc -llapack ${PETSC_SYS_LIB}

######################################################################
# TARGETS BEGIN HERE
######################################################################
all: $(BINARY_NAME)

%.o: %.cpp
	$(MPICXX) $(CPPFLAGS) -c $< -o $@

$(BINARY_NAME): $(OBJ_FILES)
	$(MPICXX) $^ $(LDFLAGS) -o $@

.PHONY = clear remake

clear:
	$(RM) $(RMFL) $(OBJ_FILES) $(BINARY_NAME) *~
remake: clear all
