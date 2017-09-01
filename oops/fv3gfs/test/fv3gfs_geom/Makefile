FC = mpif90 -f90=ifort

FCFLAGS = -fno-alias -stack_temps -safe_cray_ptr -ftz -assume byterecl -i4 -r8 -nowarn -g -sox -traceback -qopenmp -O2

RM = /bin/rm -f

NETCDF_INC = $(shell $(NETCDF)/bin/nf-config --fflags)
NETCDF_LIB = $(shell $(NETCDF)/bin/nf-config --flibs)

# Cheyenne
FMS = /glade/p/work/rahulm/GFDL/FMS
# Theia
FMS = /scratch4/NCEPDEV/global/save/Rahul.Mahajan/git/GFDL/FMS
FMS_LIB = $(FMS)/build/libfms.a
FMS_INC = -I$(FMS)/build -I$(FMS)/include
FMSCPP = -Duse_libMPI -Duse_netCDF

CPPFLAGS = $(FMS_CPP)

SRCS_geom = module_fv3gfs.f90 fv3gfs_geom.f90
SRCS_fields = module_fv3gfs.f90 fv3gfs_fields.f90
OBJS_geom = $(SRCS_geom:.f90=.o)
OBJS_fields = $(SRCS_fields:.f90=.o)

.SUFFIXES: .f90 .F90 .o

.F90.o:
	$(FC) $(CPPFLAGS) $(FCFLAGS) $(FMS_INC) $(NETCDF_INC) -c $*.F90

.f90.o:
	$(FC) $(CPPFLAGS) $(FCFLAGS) $(FMS_INC) $(NETCDF_INC) -c $*.f90

.PHONY: clean

all: fv3gfs_geom.x fv3gfs_fields.x

fv3gfs_geom.x: $(OBJS_geom)
	$(FC) $(OBJS_geom) $(CPPFLAGS) $(FCFLAGS) -I. $(FMS_INC) $(NETCDF_INC) -I. -o $@ $(FMS_LIB) $(NETCDF_LIB)

fv3gfs_fields.x: $(OBJS_fields)
	$(FC) $(OBJS_fields) $(CPPFLAGS) $(FCFLAGS) -I. $(FMS_INC) $(NETCDF_INC) -I. -o $@ $(FMS_LIB) $(NETCDF_LIB)

clean:
	$(RM) *.o *.x *.mod
