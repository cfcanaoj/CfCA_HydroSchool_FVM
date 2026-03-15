dird := snapshots

exe  := bw.x

all: ${exe}

.PHONY: all clean allclean

fc= nvfortran
#foptopenacc = -Minfo=accel -acc
foptopenacc = -Minfo=accel -acc -gpu=managed
fopt = -g -traceback -O2
#fopt=-g -traceback -check all -fpe0

${exe}: main_HDC.f90
	${fc} ${fopt} ${foptopenacc} $< -o ${exe}

#################
# clean up
#################
clean:
	rm -f ${exe} *.o *.mod *~

allclean:
	rm -f ${dirm}/*.wmv ${dirf}/*.png ${dird}/*.xss
