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
fc=ifort -diag-disable=10448 -extend-source
fopt=-g -traceback -O2
#fopt=-g -traceback -check all -fpe0

##########################################
# Timeseqencial data and figures
##########################################

######################
# directory
######################

all: ${exem1} ${exesn}

.PHONY: all clean allclean

#################
# simulation
#################

${exem1}: main_m1.f90
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
	rm -f ${exe} *.o *.mod *~

allclean:
	rm -fr ${dird} 
