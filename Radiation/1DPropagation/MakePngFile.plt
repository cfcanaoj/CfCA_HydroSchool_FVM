# Usage: gnuplot -e num  MakePngFile.plt
set terminal push
set terminal pngcairo enhanced

if (!exists("num")){
   print "Usage: gnuplot -e num=?? MakePngFile.plt"
   print "Now set num=1"
   num=1
}

set xlabel "x"

set yrange [*:*]
set xrange [0:1]

# set dir name
model="m1"
print "Plot data in ".model

inpfile=model.sprintf("/snap%05d.dat",num)
outfile=model.sprintf("/snap%05d.png",num)

# extract time from the first line
command=sprintf("awk 'NR==1{print($3)}' %s", inpfile)
time=system(command)
set output outfile
print "The plot is saved as ".outfile
set title sprintf("time = %s",time) 
plot inpfile u 1:2 ti model w l lw 6


set terminal pop
