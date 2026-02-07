# Radiation / 1DPropagation

このディレクトリには、輻射輸送の **1D 伝播（Propagation）** のサンプル計算を走らせて、出力を可視化するための一式が入っています。  

## compile
解析サーバの場合、以下のような `Makefile`でコンパイルできます。
<pre>
##########################################
# Makefile for 1Dpropataion
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
</pre>
