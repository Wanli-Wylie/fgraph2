#-
#- $Id: AA_make.ldef 1311 2013-06-04 08:13:54Z josefine.ghattas $
#-
#---------------------------------------------------------------------
#- Creation des elements relatifs a STOMATE
#---------------------------------------------------------------------
SHELL = /bin/sh
#---------------------------------------------------------------------
LIBDIR = ../../../lib
BINDIR = ../../../bin
MODDIR = $(LIBDIR)
#---------------------------------------------------------------------
ORCHIDEE_LIB = $(LIBDIR)/liborchidee.a
MODEL_LIB = $(LIBDIR)/libstomate.a
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
#- $Id: AA_make 4511 2017-07-18 14:51:56Z albert.jornet $
#-
PARALLEL_LIB = $(LIBDIR)/libparallel.a
SXPARALLEL_LIB = $(PARALLEL_LIB)
#-
PARAM_LIB = $(LIBDIR)/libparameters.a
SXPARAM_LIB = $(PARAM_LIB)
#-
ORGLOB_LIB = $(LIBDIR)/liborglob.a
SXORGLOB_LIB = $(ORGLOB_LIB)
#-
MODS1 = stomate_data.f90      \
        lpj_constraints.f90        \
        lpj_cover.f90              \
        lpj_crown.f90              \
        lpj_establish.f90          \
        lpj_fire.f90               \
        lpj_gap.f90                \
        lpj_kill.f90               \
        lpj_light.f90              \
        lpj_pftinout.f90           \
        stomate_growth_fun_all.f90          \
        stomate_data.f90           \
        stomate_io.f90             \
        stomate_litter.f90         \
        stomate_phenology.f90      \
        stomate_prescribe.f90      \
        stomate_season.f90         \
        stomate_soilcarbon.f90     \
        stomate_turnover.f90       \
        stomate_vmax.f90           \
        stomate_lcchange.f90       \
        stomate_lpj.f90            \
        stomate_resp.f90           \
        grassland_constantes.f90         \
        grassland_fonctions.f90          \
        grassland_cutting.f90          \
        grassland_fertilisation.f90   \
        grassland_grazing.f90         \
        grassland_management.f90     \
        stomate.f90

OBJSMODS1 = $(MODS1:.f90=.o)
#-
.PRECIOUS : $(MODEL_LIB)
#-
all:
	$(M_K) libparallel
	$(M_K) libparameters
	$(M_K) liborglob
	$(M_K) m_all
	@echo stomate is OK

m_all: $(MODEL_LIB)($(OBJSMODS1))

memory:

libparallel:
	(cd ../src_parallel; $(M_K) -f Makefile)

libparameters:
	(cd ../src_parameters; $(M_K) -f Makefile)

liborglob:
	(cd ../src_global; $(M_K) -f Makefile)

$(MODEL_LIB)(%.o): %.f90
	$(F_C) $(F_O) -I$(NCDF_INC) $*.f90
	$(A_C) $(MODEL_LIB) $*.o
	$(A_C) $(ORCHIDEE_LIB) $*.o
	$(RM) $*.o

config: 
	$(BINDIR)/Fparser -name STOMATE $(MODS1)
	echo 'Configuration of STOMATE done'

clean:
	$(RM) $(MODEL_LIB)

$(MODEL_LIB)(stomate.o):          \
  $(PARAM_LIB)(pft_parameters.o)  \
  $(MODEL_LIB)(stomate_io.o)      \
  $(MODEL_LIB)(stomate_lpj.o)

$(MODEL_LIB)(stomate_data.o): \
  $(PARAM_LIB)(pft_parameters.o)	\
  $(PARAM_LIB)(constantes.o)  

$(MODEL_LIB)(lpj_constraints.o): \
  $(MODEL_LIB)(stomate_data.o)

$(MODEL_LIB)(lpj_cover.o): \
  $(MODEL_LIB)(stomate_data.o)

$(MODEL_LIB)(lpj_crown.o): \
  $(MODEL_LIB)(stomate_data.o)

$(MODEL_LIB)(lpj_establish.o): \
  $(MODEL_LIB)(stomate_data.o)

$(MODEL_LIB)(lpj_fire.o): \
  $(MODEL_LIB)(stomate_data.o)

$(MODEL_LIB)(lpj_gap.o): \
  $(MODEL_LIB)(stomate_data.o)

$(MODEL_LIB)(lpj_kill.o): \
  $(MODEL_LIB)(stomate_data.o)

$(MODEL_LIB)(lpj_light.o): \
  $(MODEL_LIB)(stomate_data.o)

$(MODEL_LIB)(lpj_pftinout.o): \
  $(MODEL_LIB)(stomate_data.o)

$(MODEL_LIB)(stomate_growth_fun_all.o): \
  $(MODEL_LIB)(stomate_data.o) \
  $(PARAM_LIB)(pft_parameters.o) \
  $(PARAM_LIB)(function_library.o) 
  

$(MODEL_LIB)(stomate_io.o): \
  $(MODEL_LIB)(stomate_data.o) \
  $(PARALLEL_LIB)(mod_orchidee_para.o)


$(MODEL_LIB)(stomate_litter.o): \
  $(MODEL_LIB)(stomate_data.o)

$(MODEL_LIB)(stomate_phenology.o): \
  $(MODEL_LIB)(stomate_data.o)

$(MODEL_LIB)(stomate_prescribe.o): \
  $(MODEL_LIB)(stomate_data.o)

$(MODEL_LIB)(stomate_resp.o): \
  $(MODEL_LIB)(stomate_data.o)

$(MODEL_LIB)(stomate_season.o): \
  $(PARAM_LIB)(pft_parameters.o)	\
  $(MODEL_LIB)(stomate_data.o)

$(MODEL_LIB)(stomate_soilcarbon.o): \
  $(MODEL_LIB)(stomate_data.o) \
  $(ORGLOB_LIB)(grid.o)\
  $(ORGLOB_LIB)(solar.o)

$(MODEL_LIB)(stomate_turnover.o): \
  $(MODEL_LIB)(stomate_data.o)

$(MODEL_LIB)(stomate_vmax.o): \
  $(MODEL_LIB)(stomate_data.o)

$(MODEL_LIB)(stomate_lcchange.o): \
  $(MODEL_LIB)(stomate_data.o)

$(MODEL_LIB)(grassland_constantes.o): \
  $(PARAM_LIB)(pft_parameters.o)        \
  $(PARAM_LIB)(constantes.o)

$(MODEL_LIB)(grassland_fonctions.o): \
  $(PARAM_LIB)(pft_parameters.o)        \
  $(PARAM_LIB)(constantes.o)            \
  $(MODEL_LIB)(grassland_constantes.o)

$(MODEL_LIB)(grassland_fertilisation.o): \
  $(PARAM_LIB)(pft_parameters.o)        \
  $(PARAM_LIB)(constantes.o)            \
  $(MODEL_LIB)(grassland_constantes.o)      \
  $(MODEL_LIB)(grassland_fonctions.o)

$(MODEL_LIB)(grassland_cutting.o): \
  $(PARAM_LIB)(pft_parameters.o)        \
  $(PARAM_LIB)(constantes.o)            \
  $(MODEL_LIB)(grassland_constantes.o)      \
  $(MODEL_LIB)(grassland_fonctions.o)

$(MODEL_LIB)(grassland_grazing.o): \
  $(PARAM_LIB)(pft_parameters.o)        \
  $(PARAM_LIB)(constantes.o)            \
  $(MODEL_LIB)(grassland_constantes.o)      \
  $(MODEL_LIB)(grassland_fonctions.o)

$(MODEL_LIB)(grassland_management.o): \
  $(PARAM_LIB)(pft_parameters.o)        \
  $(PARAM_LIB)(constantes.o)            \
  $(MODEL_LIB)(grassland_constantes.o)      \
  $(MODEL_LIB)(grassland_fonctions.o)       \
  $(MODEL_LIB)(grassland_grazing.o)      \
  $(MODEL_LIB)(grassland_cutting.o)       \
  $(MODEL_LIB)(grassland_fertilisation.o)

$(MODEL_LIB)(stomate_lpj.o):            \
  $(MODEL_LIB)(lpj_constraints.o)       \
  $(MODEL_LIB)(lpj_cover.o)             \
  $(MODEL_LIB)(lpj_crown.o)             \
  $(MODEL_LIB)(lpj_establish.o)         \
  $(MODEL_LIB)(lpj_fire.o)              \
  $(MODEL_LIB)(lpj_gap.o)               \
  $(MODEL_LIB)(lpj_kill.o)              \
  $(MODEL_LIB)(lpj_light.o)             \
  $(MODEL_LIB)(lpj_pftinout.o)          \
  $(MODEL_LIB)(stomate_data.o)          \
  $(MODEL_LIB)(stomate_litter.o)        \
  $(MODEL_LIB)(stomate_growth_fun_all.o)           \
  $(MODEL_LIB)(stomate_phenology.o)     \
  $(MODEL_LIB)(stomate_prescribe.o)     \
  $(MODEL_LIB)(stomate_season.o)        \
  $(MODEL_LIB)(stomate_soilcarbon.o)    \
  $(MODEL_LIB)(stomate_turnover.o)      \
  $(MODEL_LIB)(stomate_vmax.o)          \
  $(MODEL_LIB)(stomate_resp.o)          \
  $(MODEL_LIB)(stomate_lcchange.o)      \
  $(MODEL_LIB)(grassland_management.o)


