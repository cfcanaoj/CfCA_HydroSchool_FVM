set terminal x11
set xlabel "x"


set yrange [*:*]
set xrange [*:*]

model=(ARG1)

print "Plot data in ".model

do for [i = 1:40 ] {
  inpfile=model.sprintf("/snap%05d.dat",i)

  # extract time from the first line
  command=sprintf("awk 'NR==1{print($3)}' %s", inpfile)
  time=system(command)

  set title sprintf("time = %s",time) 
  plot inpfile u 1:2 ti model w l 
  pause 0.1
}
