dird := snapshots

exe  := bw.x

all: ${exe}

.PHONY: all clean allclean

fc= mpif90
#foptopenacc = -Minfo=accel -acc
foptopenacc = -Minfo=accel -acc -gpu=mem:separate
fopt = -g -traceback -O2
#fopt=-g -traceback -check all -fpe0


config.o: config.f90
	${fc} ${fopt} ${foptopenacc} $< -c $@

mpi_setup.o: mpi_setup.f90 config.o
	${fc} ${fopt} ${foptopenacc} $< -c $@

mhd.o: mhd.f90 mpi_setup.o config.o
	${fc} ${fopt} ${foptopenacc} $< -c $@

boundary.o: boundary.f90 mhd.o mpi_setup.o config.o
	${fc} ${fopt} ${foptopenacc} $< -c $@

output.o: output.f90 mhd.o mpi_setup.o config.o
	${fc} ${fopt} ${foptopenacc} $< -c $@

main.o: main.f90 output.o boundary.o config.o
	${fc} ${fopt} ${foptopenacc} $< -c $@

objs = config.o mpi_setup.o mhd.o boundary.o output.o main.o

${exe}: ${objs}
	${fc} ${fopt} ${foptopenacc} ${objs} -o $@

#################
# clean up
#################
clean:
	rm -f ${exe} *.o *.mod *~

allclean:
	rm -f ${dirm}/*.wmv ${dirf}/*.png ${dird}/*.xss
