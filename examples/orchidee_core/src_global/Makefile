#-
#- $Id: AA_make.ldef 1311 2013-06-04 08:13:54Z josefine.ghattas $
#-
#---------------------------------------------------------------------
#- Creation des elements relatifs a ORGLOB
#---------------------------------------------------------------------
SHELL = /bin/sh
#---------------------------------------------------------------------
LIBDIR = ../../../lib
BINDIR = ../../../bin
MODDIR = $(LIBDIR)
#---------------------------------------------------------------------
ORCHIDEE_LIB = $(LIBDIR)/liborchidee.a
MODEL_LIB = $(LIBDIR)/liborglob.a
SXMODEL_LIB = $(MODEL_LIB)
#-
#- $Id: AA_make.gdef 2634 2015-08-31 11:44:24Z mafoipsl $
#-
#- Validate the correlation between the target and the environment
#-
UTIL_DIR = ../../../util
#-
#-Q-
######-Q- ada      F_O = -DCPP_PARA -p -g -traceback -fp-stack-check -ftrapuv -check bounds $(F_D) $(F_P) -I$(MODDIR) -module $(MODDIR)
######-Q- curie  F_O = -DCPP_PARA -p -g -traceback -fp-stack-check -ftrapuv -check bounds $(F_D) $(F_P) -I$(MODDIR) -module $(MODDIR)
######-Q- cur_mono  F_O = -DCPP_PARA -p -g -traceback -fp-stack-check -ftrapuv -check bounds $(F_D) $(F_P) -I$(MODDIR) -module $(MODDIR)
#-
#- Global definitions for ciclad, ifort
LIB_MPI = MPI1
LIB_MPI_BIS = MPI1
FCM_ARCH = ifort_CICLAD
M_K = make
P_C = cpp
P_O = -P -C -traditional $(P_P)
F_C = mpiifort -c -cpp
F_D =
F_P = -i4 -r8
F_O = -DCPP_PARA -O3 $(F_D) $(F_P) -I$(MODDIR) -module $(MODDIR) -fp-model precise
F_L = mpiifort
M_M = 0
L_X = 0
L_O =
A_C = ar -rs
A_G = ar -x
C_C = cc -c
C_O =
C_L = cc
#-
NCDF_INC = /usr/sw-cluster/apps/netcdf-with-parallel/include
NCDF_LIB = -L/usr/sw-cluster/apps/netcdf-with-parallel/lib -lnetcdff -lnetcdf
####-Q- lxiv8    F_O = -DCPP_PARA -p -g -traceback -fp-stack-check -ftrapuv -check bounds $(F_D) $(F_P) -I$(MODDIR) -module $(MODDIR)
######-Q- occigen  F_O = -DCPP_PARA -p -g -traceback -fp-stack-check -ftrapuv -check bounds $(F_D) $(F_P) -I$(MODDIR) -module $(MODDIR)
RM      = rm -f
STRIP   = strip
SIZE    = size

#-
#- $Id: AA_make 4063 2017-01-27 14:13:56Z albert.jornet $
#-
PARAM_LIB = $(LIBDIR)/libparameters.a
SXPARAM_LIB = $(PARAM_LIB)
#-
MODS1 = grid.f90           \
        solar.f90          \
	interpol_help.f90  \
	gauss_jordan_method.f90 \
	interpweight.f90 

OBJSMODS1 = $(MODS1:.f90=.o)
#-
.PRECIOUS : $(MODEL_LIB)
#-
all:
	$(M_K) libparallel
	$(M_K) libparameters
	$(M_K) m_all
	@echo orglob is OK

m_all: $(MODEL_LIB)($(OBJSMODS1))

libparallel:
	(cd ../src_parallel; $(M_K) -f Makefile)

libparameters:
	(cd ../src_parameters; $(M_K) -f Makefile)

$(MODEL_LIB)(%.o) : %.f90
	$(F_C) $(F_O) -I$(NCDF_INC) $*.f90
	$(A_C) $(MODEL_LIB) $*.o
	$(A_C) $(ORCHIDEE_LIB) $*.o
	$(RM) $*.o

config : 
	$(BINDIR)/Fparser -name ORGLOB $(MODS1)
	echo 'Configuration of ORGLOB done'

clean:
	$(RM) $(MODEL_LIB)

$(MODEL_LIB)(grid.o): \
  $(PARAM_LIB)(constantes.o)

$(MODEL_LIB)(interpol_help.o): \
  $(PARAM_LIB)(constantes.o)

$(MODEL_LIB)(gauss_jordan_method.o): \
  $(PARAM_LIB)(constantes.o)

$(MODEL_LIB)(interpweight.o): \
  $(PARAM_LIB)(constantes.o)
