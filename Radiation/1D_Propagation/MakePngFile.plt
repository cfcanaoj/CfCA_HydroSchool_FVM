# Usage: gnuplot -c MakePngFile.plt model num
set terminal push
set terminal pngcairo enhanced

# 引数チェック
if (ARGC < 2) {
   print "Usage: gnuplot -c MakePngFile.plt model num" 
   print "Now set model=m1 and num=1"
   model="m1"
   num=1
}else{
model = (ARG1)
num = int(ARG2)

}

set xlabel "x"

set yrange [*:*]
set xrange [0:1]

# set dir name
model="m1"

if (exists("dir")) model=dir

print "Plot data in ".model

inpfile=model.sprintf("/snap%05d.dat",num)
# below we plot E
outfile=model.sprintf("/E%05d.png",num)

# extract time from the first line
command=sprintf("awk 'NR==1{print($3)}' %s", inpfile)
time=system(command)
set output outfile
print "The plot is saved as ".outfile
set title sprintf("time = %s",time) 

plot inpfile u 1:2 ti model w l lw 6


set terminal pop
