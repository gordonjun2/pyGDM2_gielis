# Start of the makefile
# Defining variables
OBJDIR=obj

objects = $(addprefix $(OBJDIR)/, \
	precision.o \
	propagator_generalized.o \
	propagator_elec_elec_123.o \
	propagator_elec_mag_freespace.o \
	routines_decayrate.o \
	routines_incidentfields.o \
	routines_linear.o \
	routines_other.o)


	
### PRECISION
PRECISION=precision_single.f90
# PRECISION=precision_double.f90




### GFORTRAN
f90comp = gfortran
SWITCH = -O3 -mcmodel=medium -fopenmp -fPIC -fbounds-check -finit-real=snan
F2PYFLAGS=--fcompiler=gnu95 --f90flags="-O3 -mcmodel=medium -fopenmp -fbounds-check" --opt="-O3" -lgomp -DF2PY_REPORT_ON_ARRAY_COPY=1

### IFORT
# f90comp = ifort
# SWITCH = -xAVX -openmp -O3 -mcmodel=medium -xhost -fp-model fast -fPIC
# F2PYFLAGS=--fcompiler=intelem --f90flags='-openmp -O3 -mcmodel=medium -xhost -fp-model fast' -liomp5 --opt='-fast' -DF2PY_REPORT_ON_ARRAY_COPY=1



### F2PY - Python module building
f2pycomp = f2py2.7
f2pycomp3 = f2py3

F2PYFILE=propagator_generalized
MODULE=pyGDMfor
FORTFILE=$(F2PYFILE).f90
PYF=$(F2PYFILE).pyf
# PYFMAN=pyGDM_manual_single.pyf
# PYFMAN=pyGDM_manual_double.pyf




all: gf pyf3 module3
	cp $(MODULE)*.so ../pyGDM2/$(MODULE).so

python2: gf pyf2 module2
	cp $(MODULE)*.so ../pyGDM2/$(MODULE).so

# Main compilation
gf: $(objects)
	$(f90comp)  $(SWITCH) -c $(objects)

## Precision Definition
$(OBJDIR)/precision.o: $(PRECISION)
	$(f90comp) -c $(SWITCH) $(PRECISION) -o $@


	## Main Python-Interface Routines
	$(OBJDIR)/propagator_generalized.o: propagator_generalized.f90
		$(f90comp) -c $(SWITCH) propagator_generalized.f90 -o $@
	$(OBJDIR)/propagator_elec_elec_123.o: propagator_elec_elec_123.f90
		$(f90comp) -c $(SWITCH) propagator_elec_elec_123.f90 -o $@
	$(OBJDIR)/propagator_elec_mag_freespace.o: propagator_elec_mag_freespace.f90
		$(f90comp) -c $(SWITCH) propagator_elec_mag_freespace.f90 -o $@


	## Routines
	$(OBJDIR)/routines_decayrate.o: routines_decayrate.f90
		$(f90comp) -c $(SWITCH) routines_decayrate.f90 -o $@
	$(OBJDIR)/routines_incidentfields.o: routines_incidentfields.f90
		$(f90comp) -c $(SWITCH) routines_incidentfields.f90 -o $@
	$(OBJDIR)/routines_linear.o: routines_linear.f90
		$(f90comp) -c $(SWITCH) routines_linear.f90 -o $@









pyf2: $(FORTFILE)
	$(f2pycomp2) -m $(MODULE) -h $(PYF)  \
		--overwrite-signature \
		*.f*

pyf3: $(FORTFILE)
	$(f2pycomp3) -m $(MODULE) -h $(PYF)  \
		--overwrite-signature \
		*.f*


module2: $(objects)
	$(f2pycomp2) -c $(PYF) $(objects) \
	        ${F2PYFLAGS}

module3: $(objects)
	$(f2pycomp3) -c $(PYF) $(objects) \
	        ${F2PYFLAGS}


manconf2: $(objects)
	$(f2pycomp2) -c $(PYFMAN) $(objects) \
	        ${F2PYFLAGS}
	cp $(MODULE).so ../pyGDM2/

manconf3: $(objects)
	$(f2pycomp3) -c $(PYFMAN) $(objects) \
	        ${F2PYFLAGS}
	cp $(MODULE).so ../pyGDM2/






# Cleaning everything
clean:
	rm $(objects)
	rm $(MODULE)*.so
# End of the makefile
