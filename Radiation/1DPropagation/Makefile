##########################################
# Makefile for Turbulence Studies
##########################################

##########################################
# Programs
##########################################

exem1 = m1.x
exesn = sn.x
exefld = fld.x

######################	
# complile options
######################
# intel
#fc=ifort -diag-disable=10448
#fopt=-g -traceback -O2
#fopt=-g -traceback -check all -fpe0

# GNU 
fc=gfortran
fopt=-g -fbacktrace -O2
#fopt= -O0 -g -fbacktrace -fcheck=all


##########################################
# Timeseqencial data and figures
##########################################

######################
# directory
######################

all: ${exem1} ${exesn} ${exefld}

.PHONY: all clean allclean

#################
# simulation
#################

${exem1}: main_m1.f90
	mkdir -p m1mod
	${fc} ${fopt} $< -o $@
	rm *.mod

${exesn}: main_sn.f90
	${fc} ${fopt} $< -o $@
	rm *.mod

${exefld}: main_fld.f90
	${fc} ${fopt} $< -o $@
	rm *.mod

#################
# clean up
#################
clean:
	rm -f ${exem1} ${exesn} ${exefld} *.o *.mod *~ 

allclean:
	rm -fr ${dird} 
